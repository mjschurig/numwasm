#!/usr/bin/env node
/**
 * Generate function documentation JSON file
 *
 * Extracts JSDoc descriptions from:
 * 1. TypeDoc API JSON files (TypeScript documentation)
 * 2. C header files (WASM function documentation)
 * 3. Reverse mappings from TS→WASM (for WASM functions without direct docs)
 *
 * Output: public/function-docs.json
 */

import * as fs from 'fs';
import * as path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const DOCS_SITE_DIR = path.resolve(__dirname, '..');
const PUBLIC_DIR = path.join(DOCS_SITE_DIR, 'public');
const PACKAGES_DIR = path.resolve(DOCS_SITE_DIR, '..', 'packages');

/**
 * Extract first sentence/line from a description
 */
function getFirstLine(text) {
  if (!text) return '';
  // Get first sentence (ending with . or first line)
  const firstSentence = text.split(/[.\n]/)[0];
  return firstSentence.trim() + (firstSentence.endsWith('.') ? '' : '.');
}

/**
 * Extract description from TypeDoc comment object
 */
function extractDescription(comment) {
  if (!comment?.summary) return null;

  const parts = comment.summary
    .filter((part) => part.kind === 'text')
    .map((part) => part.text);

  if (parts.length === 0) return null;

  const fullText = parts.join('').trim();
  return getFirstLine(fullText);
}

/**
 * Recursively extract function/class names and descriptions
 */
function extractDocs(children, docs = {}) {
  if (!children) return docs;

  for (const child of children) {
    // Skip private/internal items
    if (child.name?.startsWith('_')) continue;

    // Try to get description from the item itself first
    let description = extractDescription(child.comment);

    // For functions (kind 64), check signatures for comments
    if (!description && child.signatures?.length > 0) {
      for (const sig of child.signatures) {
        description = extractDescription(sig.comment);
        if (description) break;
      }
    }

    // For call signatures on classes/interfaces
    if (!description && child.kind === 4096 && child.comment) {
      description = extractDescription(child.comment);
    }

    if (description && child.name) {
      docs[child.name] = description;
    }

    // Recurse into namespaces, modules, and classes
    if (child.children) {
      extractDocs(child.children, docs);
    }
  }

  return docs;
}

/**
 * Process a single package's API JSON
 */
function processPackage(packageName) {
  const apiPath = path.join(PUBLIC_DIR, `api-${packageName}.json`);

  if (!fs.existsSync(apiPath)) {
    console.warn(`  Warning: API JSON not found at ${apiPath}`);
    return {};
  }

  const data = JSON.parse(fs.readFileSync(apiPath, 'utf-8'));
  return extractDocs(data.children);
}

/**
 * Extract class-level JSDoc comments directly from TypeScript source files
 * This handles cases where TypeDoc doesn't include the class comment in the API JSON
 */
function extractClassDocsFromSource(packageDir) {
  const docs = {};
  const tsDir = path.join(packageDir, 'src', 'ts');

  if (!fs.existsSync(tsDir)) return docs;

  // Find all TypeScript files
  const tsFiles = [];
  function findTsFiles(dir) {
    const entries = fs.readdirSync(dir, { withFileTypes: true });
    for (const entry of entries) {
      const fullPath = path.join(dir, entry.name);
      if (entry.isDirectory()) {
        findTsFiles(fullPath);
      } else if (entry.name.endsWith('.ts') && !entry.name.endsWith('.d.ts')) {
        tsFiles.push(fullPath);
      }
    }
  }
  findTsFiles(tsDir);

  for (const filePath of tsFiles) {
    const content = fs.readFileSync(filePath, 'utf-8');

    // Match JSDoc comment followed by export class declaration
    // Pattern: /** ... */ export [abstract] class ClassName
    const classDocRegex = /\/\*\*\s*([\s\S]*?)\s*\*\/\s*export\s+(?:abstract\s+)?class\s+(\w+)/g;

    let match;
    while ((match = classDocRegex.exec(content)) !== null) {
      const rawComment = match[1];
      const className = match[2];

      // Clean up the comment
      const description = rawComment
        .replace(/\s*\*\s*/g, ' ')  // Remove comment asterisks
        .replace(/@\w+[\s\S]*$/, '') // Remove @tags and everything after
        .replace(/\s+/g, ' ')        // Normalize whitespace
        .trim();

      if (description && description.length > 5) {
        docs[className] = getFirstLine(description);
      }
    }
  }

  return docs;
}

// Main
console.log('Extracting function documentation...\n');

const result = {};

// All packages to process
const ALL_PACKAGES = [
  'numwasm', 'sciwasm', 'symwasm',
  'arwasm', 'lawasm', 'linwasm', 'quadwasm', 'superluwasm', 'xsfwasm', 'odewasm'
];

// Process each package
for (const pkg of ALL_PACKAGES) {
  console.log(`Processing ${pkg}...`);
  result[pkg] = processPackage(pkg);
  console.log(`  Extracted ${Object.keys(result[pkg]).length} descriptions`);
}

// Extract class docs from source files (to fill gaps in TypeDoc output)
console.log('\nExtracting class documentation from source...');
for (const pkg of ALL_PACKAGES) {
  const pkgDir = path.join(PACKAGES_DIR, pkg);
  const classDocs = extractClassDocsFromSource(pkgDir);
  let added = 0;
  for (const [name, desc] of Object.entries(classDocs)) {
    if (!result[pkg][name]) {
      result[pkg][name] = desc;
      added++;
    }
  }
  console.log(`  ${pkg}: ${added} class docs added`);
}

// Fallback descriptions for important classes/constants without docs
const fallbackDocs = {
  numwasm: {
    NDArray: 'NumPy-inspired N-dimensional array for TypeScript/WebAssembly.',
    Ellipsis: 'Placeholder for full slicing along an axis (equivalent to NumPy\'s ...).',
    Newaxis: 'Placeholder for adding a new axis (equivalent to NumPy\'s np.newaxis).',
    linalg: 'Linear algebra submodule with matrix operations.',
    fft: 'Fast Fourier Transform submodule.',
    random: 'Random number generation submodule.',
  },
  sciwasm: {
    // Add sciwasm fallbacks as needed
  },
  symwasm: {
    Expr: 'Base class for all symbolic expressions.',
    Symbol: 'Symbolic variable for algebraic expressions.',
    Integer: 'Symbolic integer type.',
    Rational: 'Symbolic rational number type.',
    Float: 'Symbolic floating-point number type.',
  },
  arwasm: {
    // ARPACK - eigenvalue solvers
  },
  lawasm: {
    // LAPACK - linear algebra routines
  },
  linwasm: {
    // LINPACK - legacy linear algebra
  },
  quadwasm: {
    // QUADPACK - numerical integration
  },
  superluwasm: {
    // SuperLU - sparse matrix solvers
  },
  xsfwasm: {
    // Special functions (Bessel, Gamma, etc.)
  },
};

console.log('\nAdding fallback documentation...');
for (const [pkg, fallbacks] of Object.entries(fallbackDocs)) {
  let added = 0;
  for (const [name, desc] of Object.entries(fallbacks)) {
    if (!result[pkg][name]) {
      result[pkg][name] = desc;
      added++;
    }
  }
  if (added > 0) {
    console.log(`  ${pkg}: ${added} fallback docs added`);
  }
}

/**
 * Find all .h files recursively in a directory
 */
function findHeaderFiles(dir, files = []) {
  if (!fs.existsSync(dir)) return files;
  const entries = fs.readdirSync(dir, { withFileTypes: true });
  for (const entry of entries) {
    const fullPath = path.join(dir, entry.name);
    if (entry.isDirectory()) {
      findHeaderFiles(fullPath, files);
    } else if (entry.name.endsWith('.h')) {
      files.push(fullPath);
    }
  }
  return files;
}

/**
 * Extract documentation from C header files
 * Parses patterns like: [JSDoc comment] function_name(...)
 */
function extractCHeaderDocs(packageDir) {
  const wasmDir = path.join(packageDir, 'src', 'wasm');
  const docs = {};

  if (!fs.existsSync(wasmDir)) {
    return docs;
  }

  const headerFiles = findHeaderFiles(wasmDir);

  for (const filePath of headerFiles) {
    const content = fs.readFileSync(filePath, 'utf-8');

    // Match: /** single-line description */ followed by function declaration
    // Only match comments that don't contain newlines (to skip file-level headers)
    // Pattern: /** [^*\n]+ */ [return_type] function_name(
    const regex = /\/\*\*\s*([^*\n][^\n]*?)\s*\*\/\s*\n?\s*(?:NDArray\*?|double|int|void|int32_t|float|size_t|bool|uint32_t|int64_t|uint64_t)\s*\*?\s*(\w+)\s*\(/g;

    let match;
    while ((match = regex.exec(content)) !== null) {
      const description = match[1].trim();
      const funcName = match[2];

      // Skip internal functions and very short descriptions
      if (funcName.startsWith('_') || description.length < 5) continue;

      // Get first sentence
      const firstSentence = getFirstLine(description);
      if (firstSentence) {
        // Store with underscore prefix (WASM export name)
        docs[`_${funcName}`] = firstSentence;
      }
    }
  }

  return docs;
}

/**
 * Build reverse mapping: WASM function -> TS function description
 */
function buildReverseWasmDocs(tsDocs, tsToWasm) {
  const wasmDocs = {};

  for (const [tsFuncName, wasmFuncs] of Object.entries(tsToWasm)) {
    // Get the TS function's description
    // Handle both "funcName" and "ClassName.methodName" patterns
    const baseName = tsFuncName.includes('.') ? tsFuncName.split('.')[0] : tsFuncName;
    const description = tsDocs[tsFuncName] || tsDocs[baseName];

    if (description && wasmFuncs) {
      for (const wasmFunc of wasmFuncs) {
        // Only set if not already documented
        if (!wasmDocs[wasmFunc]) {
          wasmDocs[wasmFunc] = description;
        }
      }
    }
  }

  return wasmDocs;
}

// Load wasm-exports.json for reverse mappings
let wasmExports = { tsToWasm: {} };
const wasmExportsPath = path.join(PUBLIC_DIR, 'wasm-exports.json');
if (fs.existsSync(wasmExportsPath)) {
  wasmExports = JSON.parse(fs.readFileSync(wasmExportsPath, 'utf-8'));
}

// Extract C header documentation for each package
console.log('\nExtracting C header documentation...');
const cHeaderDocs = {};
for (const pkg of ALL_PACKAGES) {
  const pkgDir = path.join(PACKAGES_DIR, pkg);
  cHeaderDocs[pkg] = extractCHeaderDocs(pkgDir);
  console.log(`  ${pkg}: ${Object.keys(cHeaderDocs[pkg]).length} C function docs`);
}

// Build reverse WASM mappings
console.log('\nBuilding reverse WASM mappings...');
const reverseWasmDocs = {};
for (const pkg of ALL_PACKAGES) {
  const tsToWasm = wasmExports.tsToWasm?.[pkg] || {};
  reverseWasmDocs[pkg] = buildReverseWasmDocs(result[pkg], tsToWasm);
  console.log(`  ${pkg}: ${Object.keys(reverseWasmDocs[pkg]).length} reverse mappings`);
}

// Merge all documentation sources
// Priority: C headers > reverse mappings (for WASM functions)
console.log('\nMerging documentation sources...');
for (const pkg of ALL_PACKAGES) {
  // Add C header docs (these are WASM function names like _ufunc_sqrt)
  for (const [name, desc] of Object.entries(cHeaderDocs[pkg])) {
    if (!result[pkg][name]) {
      result[pkg][name] = desc;
    }
  }

  // Add reverse mapping docs as fallback
  for (const [name, desc] of Object.entries(reverseWasmDocs[pkg])) {
    if (!result[pkg][name]) {
      result[pkg][name] = desc;
    }
  }
}

// Write output
const outputPath = path.join(PUBLIC_DIR, 'function-docs.json');
fs.writeFileSync(outputPath, JSON.stringify(result, null, 2));

console.log(`\nOutput written to: ${outputPath}`);

// Summary
const total = Object.values(result).reduce(
  (sum, pkg) => sum + Object.keys(pkg).length,
  0
);
console.log(`Total descriptions: ${total}`);

// Show file size
const stats = fs.statSync(outputPath);
console.log(`File size: ${(stats.size / 1024).toFixed(1)} KB`);
