import * as fs from 'fs';

const content = fs.readFileSync('/workspace/packages/arwasm/src/ts/eigs.ts', 'utf-8');

function findBlockEnd(content, start) {
  let braceCount = 1;
  let end = start;
  let i = start;
  let debugEvents = [];

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
    if (ch === '\`') {
      i++;
      while (i < content.length && content[i] !== '\`') {
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
    if (ch === '{') {
      braceCount++;
      if (braceCount <= 3) debugEvents.push({i, ch, braceCount, ctx: content.substring(i-10, i+20)});
    }
    else if (ch === '}') {
      braceCount--;
      if (braceCount <= 2) debugEvents.push({i, ch, braceCount, ctx: content.substring(i-10, i+20)});
    }

    end = i;
    i++;
  }

  console.log('Debug events (first 10):', debugEvents.slice(0, 10));
  return end;
}

// Find the opening brace of eigs function
const funcMatch = content.match(/export\s+async\s+function\s+eigs\s*\(/);
if (!funcMatch) {
  console.log('Function not found!');
  process.exit(1);
}

let pos = funcMatch.index + funcMatch[0].length;
let parenDepth = 1;

// Skip past parameters
while (pos < content.length && parenDepth > 0) {
  const ch = content[pos];
  if (ch === '(') parenDepth++;
  else if (ch === ')') parenDepth--;
  pos++;
}

// Skip whitespace and return type to find opening brace
let angleDepth = 0;
while (pos < content.length) {
  const ch = content[pos];
  if (ch === '<') angleDepth++;
  else if (ch === '>') angleDepth--;
  else if (ch === '{' && angleDepth === 0) {
    break;
  }
  pos++;
}

console.log('Opening brace at:', pos);
console.log('Context:', JSON.stringify(content.substring(pos - 20, pos + 30)));

const funcStart = pos + 1;
const funcEnd = findBlockEnd(content, funcStart);
const funcBody = content.substring(funcStart, funcEnd);

console.log('Function body length:', funcBody.length);
console.log('First 200 chars of body:', funcBody.substring(0, 200));
console.log('Last 200 chars of body:', funcBody.substring(funcBody.length - 200));

// Check if Module._dsaupd_ is in the body
if (funcBody.includes('Module._dsaupd_')) {
  console.log('YES: Module._dsaupd_ found in body!');
} else {
  console.log('NO: Module._dsaupd_ NOT found in body!');
  // Search the full file
  const dsaupdIdx = content.indexOf('Module._dsaupd_');
  console.log('Module._dsaupd_ in full file at:', dsaupdIdx);
  console.log('funcEnd is at:', funcStart + funcEnd - funcStart);
}
