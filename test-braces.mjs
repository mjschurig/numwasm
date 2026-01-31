// Debug script to trace findBlockEnd

const code = `export async function eigs(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  options?: EigsOptions
): Promise<EigsResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n) {
    throw new Error(\`nev must satisfy 0 < nev < n (got nev=\${nev}, n=\${n})\`);
  }

  // More code here that should be included...
  console.log("test");
}`;

function findBlockEnd(content, start) {
  let braceCount = 1;
  let end = start;
  let i = start;

  while (i < content.length && braceCount > 0) {
    const ch = content[i];
    const next = content[i + 1];

    // Debug output
    if (ch === '{' || ch === '}') {
      console.log(`  At ${i}: char='${ch}', braceCount before=${braceCount}`);
    }

    // Skip template literals (backtick) - these can contain ${} expressions
    if (ch === '\`') {
      console.log(`  Found backtick at ${i}, starting template skip`);
      i++;
      while (i < content.length && content[i] !== '\`') {
        if (content[i] === '\\\\') {
          i++; // Skip escaped character
        } else if (content[i] === '$' && content[i + 1] === '{') {
          console.log(`    Found \${} at ${i}, skipping template expression`);
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
      console.log(`  Template literal ended at ${i}`);
      i++;
      continue;
    }

    // Count braces
    if (ch === '{') braceCount++;
    else if (ch === '}') braceCount--;

    if (ch === '{' || ch === '}') {
      console.log(`    braceCount after=${braceCount}`);
    }

    end = i;
    i++;
  }

  return end;
}

// Find the opening brace
const funcStart = code.indexOf(': Promise<EigsResult> {') + ': Promise<EigsResult> {'.length;
console.log('Starting at:', funcStart);
console.log('Code around start:', code.substring(funcStart - 5, funcStart + 50));

const funcEnd = findBlockEnd(code, funcStart);
console.log('Function end at:', funcEnd);
console.log('Function body:', code.substring(funcStart, funcEnd));
