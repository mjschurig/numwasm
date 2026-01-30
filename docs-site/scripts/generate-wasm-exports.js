#!/usr/bin/env node
/**
 * Generate WASM exports JSON file
 *
 * Extracts WASM function exports from:
 * - TypeScript interface definitions (wasm-types.ts files)
 * - Build scripts (EXPORTED_FUNCTIONS arrays)
 * - TypeScript source files (to find TS→WASM function mappings)
 *
 * Output: public/wasm-exports.json
 */

import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';

/**
 * Recursively find all files matching a pattern in a directory
 */
function findFiles(dir, extension, files = []) {
  if (!fs.existsSync(dir)) return files;
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      findFiles(fullPath, extension, files);
    } else if (entry.name.endsWith(extension)) {
      files.push(fullPath);
    }
  }
  return files;
}

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const DOCS_SITE_DIR = path.resolve(__dirname, '..');
const PACKAGES_DIR = path.resolve(DOCS_SITE_DIR, '..', 'packages');

/**
 * Find all WASM function calls in a code block
 */
function findWasmCalls(code) {
  // Patterns: Module._func, wasm._func, ._wasmModule._func, wasm[wasmFn] (dynamic calls)
  const wasmCallRegex = /(?:Module|wasm|_wasmModule)(?:\._(\w+)|\[(\w+)\])\s*\(/g;
  const wasmCalls = new Set();
  let match;

  while ((match = wasmCallRegex.exec(code)) !== null) {
    const funcName = match[1] || match[2];
    if (funcName) {
      const wasmFuncName = funcName.startsWith('_') ? funcName : '_' + funcName;
      // Skip utility functions
      if (wasmFuncName !== '_free' && wasmFuncName !== '_malloc') {
        wasmCalls.add(wasmFuncName);
      }
    }
  }

  return wasmCalls;
}

/**
 * Find the end of a brace-delimited block
 */
function findBlockEnd(content, start) {
  let braceCount = 1;
  let end = start;
  for (let i = start; i < content.length && braceCount > 0; i++) {
    if (content[i] === '{') braceCount++;
    else if (content[i] === '}') braceCount--;
    end = i;
  }
  return end;
}

/**
 * Extract TS function to WASM function mappings from TypeScript source files
 * Parses exported functions, classes, and their methods to find WASM calls
 */
function extractTsToWasmMappings(packageDir) {
  const mappings = {};
  const tsDir = path.join(packageDir, 'src', 'ts');

  if (!fs.existsSync(tsDir)) {
    return mappings;
  }

  // Find all TypeScript files
  const tsFiles = findFiles(tsDir, '.ts');

  for (const filePath of tsFiles) {
    // Skip type definition files
    if (filePath.includes('types.ts') || filePath.includes('wasm-types.ts')) {
      continue;
    }

    const content = fs.readFileSync(filePath, 'utf-8');

    // 1. Find exported functions: export [async] function name(...) { ... }
    const funcRegex = /export\s+(?:async\s+)?function\s+(\w+)\s*\([^)]*\)[^{]*\{/g;
    let funcMatch;

    while ((funcMatch = funcRegex.exec(content)) !== null) {
      const funcName = funcMatch[1];
      const funcStart = funcMatch.index + funcMatch[0].length;
      const funcEnd = findBlockEnd(content, funcStart);
      const funcBody = content.substring(funcStart, funcEnd);
      const wasmCalls = findWasmCalls(funcBody);

      if (wasmCalls.size > 0) {
        mappings[funcName] = [...wasmCalls];
      }
    }

    // 2. Find exported const arrow functions: export const name = [async] (...) => { ... }
    const arrowFuncRegex = /export\s+const\s+(\w+)\s*=\s*(?:async\s*)?\([^)]*\)\s*(?::\s*[^=]+)?\s*=>\s*\{/g;
    let arrowMatch;

    while ((arrowMatch = arrowFuncRegex.exec(content)) !== null) {
      const funcName = arrowMatch[1];
      const funcStart = arrowMatch.index + arrowMatch[0].length;
      const funcEnd = findBlockEnd(content, funcStart);
      const funcBody = content.substring(funcStart, funcEnd);
      const wasmCalls = findWasmCalls(funcBody);

      if (wasmCalls.size > 0) {
        mappings[funcName] = [...wasmCalls];
      }
    }

    // 3. Find exported classes and their methods
    const classRegex = /export\s+(?:abstract\s+)?class\s+(\w+)(?:\s+extends\s+\w+)?(?:\s+implements\s+[\w,\s]+)?\s*\{/g;
    let classMatch;

    while ((classMatch = classRegex.exec(content)) !== null) {
      const className = classMatch[1];
      const classStart = classMatch.index + classMatch[0].length;
      const classEnd = findBlockEnd(content, classStart);
      const classBody = content.substring(classStart, classEnd);

      // Find all methods in the class
      // Match: methodName(...) { ... } or async methodName(...) { ... }
      // Also match: protected/private/public methodName(...)
      const methodRegex = /(?:(?:public|private|protected)\s+)?(?:async\s+)?(\w+)\s*\([^)]*\)(?:\s*:\s*[^{]+)?\s*\{/g;
      let methodMatch;

      const classWasmCalls = new Set();

      while ((methodMatch = methodRegex.exec(classBody)) !== null) {
        const methodName = methodMatch[1];
        // Skip constructor and private methods starting with _
        if (methodName === 'constructor' || methodName.startsWith('_')) {
          continue;
        }

        const methodStart = methodMatch.index + methodMatch[0].length;
        const methodEnd = findBlockEnd(classBody, methodStart);
        const methodBody = classBody.substring(methodStart, methodEnd);
        const wasmCalls = findWasmCalls(methodBody);

        // Add method-level mappings as ClassName.methodName
        if (wasmCalls.size > 0) {
          const qualifiedName = `${className}.${methodName}`;
          mappings[qualifiedName] = [...wasmCalls];
          // Also add to class-level aggregate
          for (const call of wasmCalls) {
            classWasmCalls.add(call);
          }
        }
      }

      // Add class-level mapping (all WASM calls from all methods)
      if (classWasmCalls.size > 0) {
        mappings[className] = [...classWasmCalls];
      }
    }
  }

  return mappings;
}

/**
 * Extract WASM function names from a TypeScript interface definition
 * Matches lines like: _functionName(params): returnType;
 */
function extractFromTypeScriptInterface(content) {
  const exports = [];
  // Match method definitions starting with underscore
  // Handles both single-line and multi-line definitions
  const regex = /^\s+(_\w+)\s*\(/gm;
  let match;
  while ((match = regex.exec(content)) !== null) {
    exports.push(match[1]);
  }
  return [...new Set(exports)]; // Deduplicate
}

/**
 * Extract WASM function names from a build script's EXPORTED_FUNCTIONS
 * Matches strings like: "_functionName_", "_functionName"
 */
function extractFromBuildScript(content) {
  const exports = [];
  // Find EXPORTED_FUNCTIONS block
  const exportedMatch = content.match(/EXPORTED_FUNCTIONS\s*=\s*'\[([^\]]+)\]'/s);
  if (exportedMatch) {
    const functionsStr = exportedMatch[1];
    // Match quoted function names
    const regex = /"(_\w+_?)"/g;
    let match;
    while ((match = regex.exec(functionsStr)) !== null) {
      // Skip _malloc and _free as they're utility functions
      if (match[1] !== '_malloc' && match[1] !== '_free') {
        exports.push(match[1]);
      }
    }
  }
  return [...new Set(exports)];
}

/**
 * Extract TypeScript exports from API JSON file
 */
function extractTypeScriptExports(apiJsonPath) {
  if (!fs.existsSync(apiJsonPath)) {
    console.warn(`  Warning: API JSON not found at ${apiJsonPath}`);
    return [];
  }

  const data = JSON.parse(fs.readFileSync(apiJsonPath, 'utf-8'));
  const exports = [];

  // ReflectionKind values from TypeDoc
  const FUNCTION = 64;
  const CLASS = 128;
  const VARIABLE = 32;
  const TYPE_ALIAS = 2097152;

  function processChildren(children) {
    if (!children) return;
    for (const child of children) {
      // Include functions, classes, variables, and type aliases
      if ([FUNCTION, CLASS, VARIABLE, TYPE_ALIAS].includes(child.kind)) {
        // Skip internal/private exports
        if (!child.name.startsWith('_')) {
          exports.push({
            name: child.name,
            kind: child.kind === FUNCTION ? 'function' :
                  child.kind === CLASS ? 'class' :
                  child.kind === VARIABLE ? 'variable' : 'type'
          });
        }
      }
      // Recurse into namespaces
      if (child.children) {
        processChildren(child.children);
      }
    }
  }

  processChildren(data.children);
  return exports;
}

// Main extraction
const result = {
  packages: {},
  typescript: {},
  tsToWasm: {}
};

console.log('Extracting WASM exports...\n');

// numwasm
console.log('Processing numwasm...');
const numwasmTypesPath = path.join(PACKAGES_DIR, 'numwasm', 'src', 'ts', 'types.ts');
if (fs.existsSync(numwasmTypesPath)) {
  const content = fs.readFileSync(numwasmTypesPath, 'utf-8');
  const exports = extractFromTypeScriptInterface(content);
  result.packages.numwasm = {
    numjs: exports
  };
  console.log(`  numjs.wasm: ${exports.length} functions`);
} else {
  console.warn(`  Warning: types.ts not found at ${numwasmTypesPath}`);
}

// numwasm TypeScript exports
const numwasmApiPath = path.join(DOCS_SITE_DIR, 'public', 'api-numwasm.json');
result.typescript.numwasm = extractTypeScriptExports(numwasmApiPath);
console.log(`  TypeScript exports: ${result.typescript.numwasm.length}`);

// numwasm TS→WASM mappings
const numwasmDir = path.join(PACKAGES_DIR, 'numwasm');
result.tsToWasm.numwasm = extractTsToWasmMappings(numwasmDir);
console.log(`  TS→WASM mappings: ${Object.keys(result.tsToWasm.numwasm).length} functions`);

// sciwasm
console.log('\nProcessing sciwasm...');
result.packages.sciwasm = {};

// Main sciwasm module
const sciwasmTypesPath = path.join(PACKAGES_DIR, 'sciwasm', 'src', 'ts', 'wasm-types.ts');
if (fs.existsSync(sciwasmTypesPath)) {
  const content = fs.readFileSync(sciwasmTypesPath, 'utf-8');
  const exports = extractFromTypeScriptInterface(content);
  result.packages.sciwasm.sciwasm = exports;
  console.log(`  sciwasm.wasm: ${exports.length} functions`);
}

// ARPACK module
const arpackBuildPath = path.join(PACKAGES_DIR, 'sciwasm', 'scripts', 'build-arpack.sh');
if (fs.existsSync(arpackBuildPath)) {
  const content = fs.readFileSync(arpackBuildPath, 'utf-8');
  const exports = extractFromBuildScript(content);
  result.packages.sciwasm.arpack = exports;
  console.log(`  arpack.wasm: ${exports.length} functions`);
}

// QUADPACK module
const quadpackBuildPath = path.join(PACKAGES_DIR, 'sciwasm', 'scripts', 'build-quadpack.sh');
if (fs.existsSync(quadpackBuildPath)) {
  const content = fs.readFileSync(quadpackBuildPath, 'utf-8');
  const exports = extractFromBuildScript(content);
  result.packages.sciwasm.quadpack = exports;
  console.log(`  quadpack.wasm: ${exports.length} functions`);
}

// LINPACK module
const linpackBuildPath = path.join(PACKAGES_DIR, 'sciwasm', 'scripts', 'build-linpack.sh');
if (fs.existsSync(linpackBuildPath)) {
  const content = fs.readFileSync(linpackBuildPath, 'utf-8');
  const exports = extractFromBuildScript(content);
  result.packages.sciwasm.linpack = exports;
  console.log(`  linpack.wasm: ${exports.length} functions`);
}

// SuperLU module
const superluBuildPath = path.join(PACKAGES_DIR, 'sciwasm', 'scripts', 'build-superlu.sh');
if (fs.existsSync(superluBuildPath)) {
  const content = fs.readFileSync(superluBuildPath, 'utf-8');
  const exports = extractFromBuildScript(content);
  result.packages.sciwasm.superlu = exports;
  console.log(`  superlu.wasm: ${exports.length} functions`);
}

// sciwasm TypeScript exports
const sciwasmApiPath = path.join(DOCS_SITE_DIR, 'public', 'api-sciwasm.json');
result.typescript.sciwasm = extractTypeScriptExports(sciwasmApiPath);
console.log(`  TypeScript exports: ${result.typescript.sciwasm.length}`);

// sciwasm TS→WASM mappings
const sciwasmDir = path.join(PACKAGES_DIR, 'sciwasm');
result.tsToWasm.sciwasm = extractTsToWasmMappings(sciwasmDir);
console.log(`  TS→WASM mappings: ${Object.keys(result.tsToWasm.sciwasm).length} functions`);

// symwasm
console.log('\nProcessing symwasm...');
const symwasmTypesPath = path.join(PACKAGES_DIR, 'symwasm', 'src', 'ts', 'wasm-types.ts');
if (fs.existsSync(symwasmTypesPath)) {
  const content = fs.readFileSync(symwasmTypesPath, 'utf-8');
  const exports = extractFromTypeScriptInterface(content);
  result.packages.symwasm = {
    symwasm: exports
  };
  console.log(`  symwasm.wasm: ${exports.length} functions`);
}

// symwasm TypeScript exports
const symwasmApiPath = path.join(DOCS_SITE_DIR, 'public', 'api-symwasm.json');
result.typescript.symwasm = extractTypeScriptExports(symwasmApiPath);
console.log(`  TypeScript exports: ${result.typescript.symwasm.length}`);

// symwasm TS→WASM mappings
const symwasmDir = path.join(PACKAGES_DIR, 'symwasm');
result.tsToWasm.symwasm = extractTsToWasmMappings(symwasmDir);
console.log(`  TS→WASM mappings: ${Object.keys(result.tsToWasm.symwasm).length} functions`);

// Write output
const outputPath = path.join(DOCS_SITE_DIR, 'public', 'wasm-exports.json');
fs.writeFileSync(outputPath, JSON.stringify(result, null, 2));

console.log(`\nOutput written to: ${outputPath}`);

// Summary
console.log('\n=== Summary ===');
let totalWasm = 0;
let totalTs = 0;
let totalMappings = 0;
for (const [, modules] of Object.entries(result.packages)) {
  for (const [, funcs] of Object.entries(modules)) {
    totalWasm += funcs.length;
  }
}
for (const [, exports] of Object.entries(result.typescript)) {
  totalTs += exports.length;
}
for (const [, mappings] of Object.entries(result.tsToWasm)) {
  totalMappings += Object.keys(mappings).length;
}
console.log(`Total WASM functions: ${totalWasm}`);
console.log(`Total TypeScript exports: ${totalTs}`);
console.log(`Total TS→WASM mappings: ${totalMappings}`);
