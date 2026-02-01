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
  // Patterns to match:
  // - Module._func(  - Emscripten Module
  // - module._func(  - lowercase module variable
  // - wasm._func(    - wasm variable
  // - this._module._func(  - instance module reference
  // - arr._module._func(   - any variable._module reference
  // - _wasmModule._func(   - underscore prefixed
  // - wasm[wasmFn]( or module[funcName]( - dynamic calls
  // - xsf._wasm_func( - xsfwasm pattern
  // - arpack._func( - arwasm pattern
  // - anyVar._wasm_func( - generic WASM function call pattern
  const wasmCallRegex = /(?:(?:this\.|[\w]+\.)?_?[mM]odule|wasm|_wasmModule|xsf|arpack|lapack|linpack|quadpack|superlu|ode)(?:\._(\w+)|\[['"]?(\w+)['"]?\])\s*\(/g;
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

  // Also match any variable calling _wasm_* or _* functions (common pattern)
  // This catches patterns like: someVar._wasm_gamma(x)
  const genericWasmRegex = /\w+\.(_wasm_\w+)\s*\(/g;
  while ((match = genericWasmRegex.exec(code)) !== null) {
    const funcName = match[1];
    if (funcName && funcName !== '_free' && funcName !== '_malloc') {
      wasmCalls.add(funcName);
    }
  }

  return wasmCalls;
}

/**
 * Find the end of a brace-delimited block
 * Properly handles strings, template literals, and comments
 */
function findBlockEnd(content, start) {
  let braceCount = 1;
  let end = start;
  let i = start;

  while (i < content.length && braceCount > 0) {
    const ch = content[i];
    const next = content[i + 1];

    // Skip single-line comments
    if (ch === '/' && next === '/') {
      while (i < content.length && content[i] !== '\n') i++;
      i++;
      continue;
    }

    // Skip multi-line comments
    if (ch === '/' && next === '*') {
      i += 2;
      while (i < content.length - 1 && !(content[i] === '*' && content[i + 1] === '/')) i++;
      i += 2;
      continue;
    }

    // Skip string literals (single quote)
    if (ch === "'") {
      i++;
      while (i < content.length && content[i] !== "'") {
        if (content[i] === '\\') i++; // Skip escaped character
        i++;
      }
      i++;
      continue;
    }

    // Skip string literals (double quote)
    if (ch === '"') {
      i++;
      while (i < content.length && content[i] !== '"') {
        if (content[i] === '\\') i++; // Skip escaped character
        i++;
      }
      i++;
      continue;
    }

    // Skip template literals (backtick) - these can contain ${} expressions
    if (ch === '`') {
      i++;
      while (i < content.length && content[i] !== '`') {
        if (content[i] === '\\') {
          i++; // Skip escaped character
        } else if (content[i] === '$' && content[i + 1] === '{') {
          // Handle template expression - find matching }
          i += 2;
          let templateBraceCount = 1;
          while (i < content.length && templateBraceCount > 0) {
            if (content[i] === '{') templateBraceCount++;
            else if (content[i] === '}') templateBraceCount--;
            if (templateBraceCount > 0) i++;
          }
        }
        i++;
      }
      i++;
      continue;
    }

    // Count braces
    if (ch === '{') braceCount++;
    else if (ch === '}') braceCount--;

    end = i;
    i++;
  }

  return end;
}

/**
 * Extract TS function to WASM function mappings from TypeScript source files
 * Parses exported functions, classes, and their methods to find WASM calls.
 *
 * For files where exports call internal helper functions that in turn call WASM,
 * we associate all WASM calls in the file with all exports from that file.
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

    // Find all WASM calls in the entire file (for fallback association)
    const fileWasmCalls = findWasmCalls(content);

    // Collect exported function names from this file
    const exportedFunctions = [];

    // 1. Find exported functions: export [async] function name(...)
    // Note: We find the function declaration, then scan forward to find the opening brace
    // that starts the function body (skipping braces in return type annotations)
    const funcRegex = /export\s+(?:async\s+)?function\s+(\w+)\s*\(/g;
    let funcMatch;

    while ((funcMatch = funcRegex.exec(content)) !== null) {
      const funcName = funcMatch[1];
      exportedFunctions.push(funcName);

      // Find the opening brace of the function body by tracking paren/angle bracket depth
      let pos = funcMatch.index + funcMatch[0].length;
      let parenDepth = 1; // We're inside the params already
      let angleDepth = 0;

      // Skip past parameters
      while (pos < content.length && parenDepth > 0) {
        const ch = content[pos];
        if (ch === '(') parenDepth++;
        else if (ch === ')') parenDepth--;
        pos++;
      }

      // Now find the opening brace, tracking angle brackets for generic types
      while (pos < content.length) {
        const ch = content[pos];
        if (ch === '<') angleDepth++;
        else if (ch === '>') angleDepth--;
        else if (ch === '{' && angleDepth === 0) {
          // Found the function body opening brace
          break;
        }
        pos++;
      }

      if (pos >= content.length) continue;

      const funcStart = pos + 1;
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
      exportedFunctions.push(funcName);

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
      exportedFunctions.push(className);

      const classStart = classMatch.index + classMatch[0].length;
      const classEnd = findBlockEnd(content, classStart);
      const classBody = content.substring(classStart, classEnd);

      // Find all methods in the class
      // Match: methodName(...) { ... } or async methodName(...) { ... }
      // Also match: protected/private/public/static methodName(...)
      const methodRegex = /(?:(?:public|private|protected)\s+)?(?:static\s+)?(?:async\s+)?(\w+)\s*\([^)]*\)(?:\s*:\s*[^{]+)?\s*\{/g;
      let methodMatch;

      const classWasmCalls = new Set();

      while ((methodMatch = methodRegex.exec(classBody)) !== null) {
        const methodName = methodMatch[1];
        // Skip constructor, private methods starting with _, and JavaScript keywords
        const jsKeywords = ['if', 'for', 'while', 'switch', 'catch', 'with', 'return', 'throw', 'new', 'delete', 'typeof', 'void', 'function'];
        if (methodName === 'constructor' || methodName.startsWith('_') || jsKeywords.includes(methodName)) {
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

    // 4. Associate file-level WASM calls with exports that don't have direct WASM calls
    // This handles the case where exports call internal helper functions that call WASM
    if (fileWasmCalls.size > 0 && exportedFunctions.length > 0) {
      for (const funcName of exportedFunctions) {
        // If this export has no direct WASM calls, associate it with file-level calls
        if (!mappings[funcName]) {
          mappings[funcName] = [...fileWasmCalls];
        }
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
const numwasmTypesPath = path.join(PACKAGES_DIR, 'numwasm', 'src', 'ts', '_core', 'types.ts');
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

// Note: ARPACK, QUADPACK, LINPACK, SuperLU modules have been moved to standalone packages
// (arwasm, quadwasm, linwasm, superluwasm) and are handled in the standalonePackages section below

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

// New standalone WASM packages (arwasm, lawasm, linwasm, quadwasm, superluwasm, xsfwasm, odewasm)
const standalonePackages = [
  { name: 'arwasm', wasmName: 'arpack' },
  { name: 'lawasm', wasmName: 'lapack' },
  { name: 'linwasm', wasmName: 'linpack' },
  { name: 'quadwasm', wasmName: 'quadpack' },
  { name: 'superluwasm', wasmName: 'superlu' },
  { name: 'xsfwasm', wasmName: 'xsf' },
  { name: 'odewasm', wasmName: 'ode' },
];

for (const pkg of standalonePackages) {
  console.log(`\nProcessing ${pkg.name}...`);

  // Try multiple possible locations for types.ts
  const possiblePaths = [
    path.join(PACKAGES_DIR, pkg.name, 'src', 'types.ts'),
    path.join(PACKAGES_DIR, pkg.name, 'src', 'ts', 'types.ts'),
    path.join(PACKAGES_DIR, pkg.name, 'src', 'ts', 'core', 'types.ts'),
  ];

  let typesPath = possiblePaths.find(p => fs.existsSync(p));

  if (typesPath) {
    const content = fs.readFileSync(typesPath, 'utf-8');
    const exports = extractFromTypeScriptInterface(content);
    result.packages[pkg.name] = {
      [pkg.wasmName]: exports
    };
    console.log(`  ${pkg.wasmName}.wasm: ${exports.length} functions`);
  } else {
    console.warn(`  Warning: types.ts not found in ${pkg.name}`);
    result.packages[pkg.name] = {};
  }

  // TypeScript exports (if api-*.json exists)
  const apiPath = path.join(DOCS_SITE_DIR, 'public', `api-${pkg.name}.json`);
  result.typescript[pkg.name] = extractTypeScriptExports(apiPath);
  console.log(`  TypeScript exports: ${result.typescript[pkg.name].length}`);

  // TS→WASM mappings
  const pkgDir = path.join(PACKAGES_DIR, pkg.name);
  result.tsToWasm[pkg.name] = extractTsToWasmMappings(pkgDir);
  console.log(`  TS→WASM mappings: ${Object.keys(result.tsToWasm[pkg.name]).length} functions`);
}

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
