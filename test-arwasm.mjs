import * as fs from 'fs';
import * as path from 'path';

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

function findWasmCalls(code) {
  const wasmCallRegex = /(?:(?:this\.|[\w]+\.)?_?[mM]odule|wasm|_wasmModule|xsf|arpack|lapack|linpack|quadpack|superlu|ode)(?:\._(\w+)|\[['"]?(\w+)['"]?\])\s*\(/g;
  const wasmCalls = new Set();
  let match;

  while ((match = wasmCallRegex.exec(code)) !== null) {
    const funcName = match[1] || match[2];
    if (funcName) {
      const wasmFuncName = funcName.startsWith('_') ? funcName : '_' + funcName;
      if (wasmFuncName !== '_free' && wasmFuncName !== '_malloc') {
        wasmCalls.add(wasmFuncName);
      }
    }
  }

  const genericWasmRegex = /\w+\.(_wasm_\w+)\s*\(/g;
  while ((match = genericWasmRegex.exec(code)) !== null) {
    const funcName = match[1];
    if (funcName && funcName !== '_free' && funcName !== '_malloc') {
      wasmCalls.add(funcName);
    }
  }

  return wasmCalls;
}

const content = fs.readFileSync('/workspace/packages/arwasm/src/ts/eigs.ts', 'utf-8');
console.log('Full file WASM calls:', findWasmCalls(content));

// Now test the full extractTsToWasmMappings logic
function extractTsToWasmMappings(packageDir) {
  const mappings = {};
  const tsDir = path.join(packageDir, 'src', 'ts');

  console.log('Looking for TS files in:', tsDir);
  console.log('Exists:', fs.existsSync(tsDir));

  if (!fs.existsSync(tsDir)) {
    return mappings;
  }

  const tsFiles = findFiles(tsDir, '.ts');
  console.log('Found TS files:', tsFiles);

  for (const filePath of tsFiles) {
    if (filePath.includes('types.ts') || filePath.includes('wasm-types.ts')) {
      console.log('Skipping:', filePath);
      continue;
    }
    console.log('Processing:', filePath);

    const content = fs.readFileSync(filePath, 'utf-8');

    // Find exported functions
    const funcRegex = /export\s+(?:async\s+)?function\s+(\w+)\s*\(/g;
    let funcMatch;

    while ((funcMatch = funcRegex.exec(content)) !== null) {
      const funcName = funcMatch[1];
      console.log('  Found exported function:', funcName);

      // Find the opening brace of the function body
      let pos = funcMatch.index + funcMatch[0].length;
      let parenDepth = 1;
      let angleDepth = 0;

      // Skip past parameters
      while (pos < content.length && parenDepth > 0) {
        const ch = content[pos];
        if (ch === '(') parenDepth++;
        else if (ch === ')') parenDepth--;
        pos++;
      }

      // Find the opening brace
      while (pos < content.length) {
        const ch = content[pos];
        if (ch === '<') angleDepth++;
        else if (ch === '>') angleDepth--;
        else if (ch === '{' && angleDepth === 0) {
          break;
        }
        pos++;
      }

      if (pos >= content.length) {
        console.log('  Could not find function body');
        continue;
      }

      const funcStart = pos + 1;
      const funcEnd = findBlockEnd(content, funcStart);
      const funcBody = content.substring(funcStart, funcEnd);
      console.log('  Function body length:', funcBody.length);

      const wasmCalls = findWasmCalls(funcBody);
      console.log('  WASM calls found:', wasmCalls);

      if (wasmCalls.size > 0) {
        mappings[funcName] = [...wasmCalls];
      }
    }
  }

  return mappings;
}

const arwasmDir = '/workspace/packages/arwasm';
const mappings = extractTsToWasmMappings(arwasmDir);
console.log('Final mappings:', mappings);
