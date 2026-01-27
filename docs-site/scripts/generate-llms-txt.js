/**
 * Generate llms.txt and llms-full.txt for AI agent access.
 * Reads api.json and produces structured markdown files following the llmstxt.org spec.
 *
 * Usage: node scripts/generate-llms-txt.js
 * Run after: vite build (so dist/ exists)
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const distDir = path.resolve(__dirname, '../dist');
const apiJsonPath = path.resolve(__dirname, '../public/api.json');

// ---------------------------------------------------------------------------
// ReflectionKind constants (mirrored from typedoc.ts)
// ---------------------------------------------------------------------------

const ReflectionKind = {
  Namespace: 4,
  Enum: 8,
  Variable: 32,
  Function: 64,
  Class: 128,
  Interface: 256,
  Property: 1024,
  Method: 2048,
  TypeAlias: 2097152,
  Reference: 4194304,
};

// ---------------------------------------------------------------------------
// Module & Category definitions (mirrored from apiData.ts)
// ---------------------------------------------------------------------------

const MODULE_DEFS = [
  { slug: 'linalg', displayName: 'Linear Algebra', varName: 'linalg' },
  { slug: 'fft', displayName: 'FFT', varName: 'fftModule' },
  {
    slug: 'random',
    displayName: 'Random',
    varName: null,
    names: [
      'Generator', 'BitGenerator', 'MT19937', 'PCG64', 'Philox', 'SFC64', 'SeedSequence',
      'default_rng', 'random', 'randint', 'randn', 'seed', 'initRandom',
      'getBitGenerator', 'listBitGenerators',
      'MT19937State', 'PCG64State', 'PhiloxState', 'SFC64State',
    ],
  },
  { slug: 'ma', displayName: 'Masked Arrays', varName: 'ma' },
  { slug: 'rec', displayName: 'Record Arrays', varName: 'rec' },
  { slug: 'strings', displayName: 'Strings', varName: 'strings' },
  {
    slug: 'polynomial',
    displayName: 'Polynomials',
    varName: null,
    namePatterns: [/^poly/, /^cheb/, /^herm/, /^lag/, /^leg/],
    names: [
      'ABCPolyBase', 'Chebyshev', 'Hermite', 'HermiteE', 'Laguerre', 'Legendre', 'Polynomial',
      'PolyDomainWarning', 'PolyError',
      'as_series', 'getdomain', 'mapdomain', 'mapparms', 'maxpower', 'trimcoef', 'trimseq',
    ],
  },
  { slug: 'testing', displayName: 'Testing', varName: null, isNamespace: true },
];

const CATEGORY_DEFS = [
  {
    title: 'Array Creation',
    names: [
      'NDArray', 'asarray', 'frombuffer', 'fromfile', 'fromregex', 'genfromtxt',
      'loadtxt', 'load', 'meshgrid', 'indices', 'ix_',
      'diag_indices', 'tril_indices', 'triu_indices',
      'broadcastArrays', 'broadcastShapes', 'broadcastShapesMulti', 'broadcastTo',
    ],
  },
  {
    title: 'Array Manipulation',
    names: [
      'concatenate', 'stack', 'hstack', 'vstack', 'dstack', 'column_stack', 'row_stack', 'block',
      'split', 'array_split', 'hsplit', 'vsplit', 'dsplit', 'unstack',
      'repeat', 'tile', 'resize', 'append', 'insert', 'deleteArr',
      'flip', 'fliplr', 'flipud', 'roll', 'rot90',
      'atleast_1d', 'atleast_2d', 'atleast_3d',
      'pad', 'trim_zeros', 'diagonal',
      'copyto', 'put', 'putmask', 'place', 'take', 'takeFlat', 'take_along_axis', 'put_along_axis',
      'compress', 'extract', 'select', 'choose',
      'ediff1d', 'unique', 'uniqueAll', 'uniqueCounts', 'uniqueIndex', 'uniqueInverse', 'uniqueValues',
    ],
  },
  {
    title: 'Indexing',
    names: [
      'ndenumerate', 'ndindex', 'nditer', 'flatnonzero', 'nonzero', 'argwhere',
      'unravelIndex', 'ravelMultiIndex',
      'slice', 'Slice', 'buildIndexSpecs', 'expandEllipsis',
    ],
  },
  {
    title: 'Math',
    names: [
      'abs', 'absolute', 'add', 'subtract', 'multiply', 'divide', 'true_divide',
      'floor_divide', 'negative', 'positive', 'power', 'mod', 'fmod', 'remainder', 'divmod',
      'reciprocal', 'square', 'sqrt', 'cbrt',
      'exp', 'exp2', 'expm1', 'log', 'log2', 'log10', 'log1p', 'logaddexp', 'logaddexp2',
      'sin', 'cos', 'tan', 'arcsin', 'arccos', 'arctan', 'arctan2',
      'sinh', 'cosh', 'tanh', 'arcsinh', 'arccosh', 'arctanh',
      'hypot', 'degrees', 'radians', 'deg2rad', 'rad2deg',
      'ceil', 'floor', 'trunc', 'round', 'rint', 'spacing',
      'sign', 'signbit', 'copysign', 'heaviside', 'ldexp', 'frexp', 'modf',
      'fmax', 'fmin', 'maximum', 'minimum',
      'gcd', 'lcm', 'nextafter',
      'nan_to_num', 'sinc', 'i0', 'clip',
    ],
  },
  {
    title: 'Logic',
    names: [
      'all', 'any', 'logical_and', 'logical_or', 'logical_not', 'logical_xor',
      'where', 'piecewise',
    ],
  },
  {
    title: 'Comparison',
    names: [
      'equal', 'not_equal', 'greater', 'less', 'greater_equal', 'less_equal',
      'allclose', 'isclose', 'array_equal', 'array_equiv',
      'isfinite', 'isinf', 'isnan', 'isneginf', 'isposinf',
      'isreal', 'isrealobj', 'iscomplex', 'iscomplexobj', 'isscalar', 'isfortran',
    ],
  },
  {
    title: 'Statistics',
    names: [
      'mean', 'median', 'std', 'var_', 'variance',
      'min', 'max', 'sum', 'prod',
      'cumsum', 'cumprod',
      'nanmean', 'nanmedian', 'nanstd', 'nanvar',
      'nanmin', 'nanmax', 'nansum', 'nanprod',
      'nancumsum', 'nancumprod',
      'nanargmin', 'nanargmax', 'nanpercentile', 'nanquantile',
      'histogram', 'histogram2d', 'histogramdd', 'histogram_bin_edges',
      'bincount', 'digitize',
    ],
  },
  {
    title: 'Sorting & Searching',
    names: [
      'sort', 'argsort', 'argmax', 'argmin',
      'searchsorted', 'partition', 'argpartition',
      'countNonzero',
    ],
  },
  {
    title: 'Set Operations',
    names: ['union1d', 'intersect1d', 'setdiff1d', 'setxor1d', 'in1d', 'isin'],
  },
  {
    title: 'Bitwise',
    names: [
      'bitwise_and', 'bitwise_or', 'bitwise_xor', 'bitwise_not', 'bitwise_count',
      'left_shift', 'right_shift', 'invert',
    ],
  },
  {
    title: 'I/O',
    names: [
      'save', 'savetxt', 'openMemmap',
      'array2string', 'arrayRepr', 'arrayStr',
      'formatFloatPositional', 'formatFloatScientific', 'formatValue',
      'getPrintoptions', 'setPrintoptions', 'resetPrintoptions', 'withPrintoptions',
      'baseRepr', 'baseReprArray', 'binaryRepr', 'binaryReprArray',
      'hexRepr', 'octalRepr', 'fromBaseRepr', 'fromBinaryRepr',
    ],
  },
  {
    title: 'Window Functions',
    names: ['bartlett', 'blackman', 'hamming', 'hanning', 'kaiser'],
  },
];

// ---------------------------------------------------------------------------
// Load and process api.json
// ---------------------------------------------------------------------------

const apiData = JSON.parse(fs.readFileSync(apiJsonPath, 'utf-8'));
const allChildren = apiData.children || [];

// Build module children map (same logic as apiData.ts)
const moduleClaimedNames = new Set();
const moduleChildrenMap = new Map();

for (const mod of MODULE_DEFS) {
  let children = [];

  if (mod.isNamespace) {
    const ns = allChildren.find(c => c.kind === ReflectionKind.Namespace && c.name === mod.slug);
    if (ns?.children) children = ns.children;
    moduleClaimedNames.add(mod.slug);
  } else if (mod.varName) {
    const varItem = allChildren.find(c => c.kind === ReflectionKind.Variable && c.name === mod.varName);
    if (varItem) {
      const decl = varItem.type?.declaration;
      if (decl?.children) children = decl.children;
      moduleClaimedNames.add(mod.varName);
    }
  }

  if (mod.names) {
    const nameSet = new Set(mod.names);
    for (const child of allChildren) {
      if (nameSet.has(child.name) && child.kind !== ReflectionKind.Reference) {
        children.push(child);
        moduleClaimedNames.add(child.name);
      }
    }
  }
  if (mod.namePatterns) {
    for (const child of allChildren) {
      if (
        child.kind !== ReflectionKind.Reference &&
        mod.namePatterns.some(p => p.test(child.name)) &&
        !moduleClaimedNames.has(child.name)
      ) {
        children.push(child);
        moduleClaimedNames.add(child.name);
      }
    }
  }

  moduleChildrenMap.set(mod.slug, children);
}

// Build categorised top-level items
const categoryClaimedNames = new Set();
for (const cat of CATEGORY_DEFS) {
  for (const n of cat.names) categoryClaimedNames.add(n);
}

function getCategorisedTopLevelItems() {
  const result = [];

  for (const cat of CATEGORY_DEFS) {
    const nameSet = new Set(cat.names);
    const items = allChildren.filter(
      c => nameSet.has(c.name) && c.kind !== ReflectionKind.Reference && !moduleClaimedNames.has(c.name)
    );
    if (items.length > 0) result.push({ title: cat.title, items });
  }

  // Types
  const typeItems = allChildren.filter(
    c =>
      (c.kind === ReflectionKind.Interface || c.kind === ReflectionKind.TypeAlias || c.kind === ReflectionKind.Enum) &&
      !moduleClaimedNames.has(c.name)
  );
  if (typeItems.length > 0) result.push({ title: 'Types & Interfaces', items: typeItems });

  // Constants
  const constantNames = new Set([
    'e', 'pi', 'inf', 'nan', 'euler_gamma', 'newaxis', 'ellipsis',
    'NAN', 'NINF', 'NZERO', 'PINF', 'PZERO',
  ]);
  const constants = allChildren.filter(
    c => c.kind === ReflectionKind.Variable && constantNames.has(c.name) && !moduleClaimedNames.has(c.name)
  );
  if (constants.length > 0) result.push({ title: 'Constants', items: constants });

  return result;
}

// ---------------------------------------------------------------------------
// Helpers for extracting documentation from reflections
// ---------------------------------------------------------------------------

function getSignature(item) {
  return item.signatures?.[0] ?? item.type?.declaration?.signatures?.[0] ?? null;
}

function getDescription(item) {
  const summary =
    item.comment?.summary ||
    item.signatures?.[0]?.comment?.summary ||
    item.type?.declaration?.signatures?.[0]?.comment?.summary;
  if (!summary) return '';
  return summary.map(p => p.text).join('').replace(/\n/g, ' ').trim();
}

function formatType(typeInfo) {
  if (!typeInfo) return 'unknown';

  switch (typeInfo.type) {
    case 'intrinsic':
      return typeInfo.name || 'unknown';
    case 'reference':
      if (typeInfo.typeArguments?.length) {
        return `${typeInfo.name}<${typeInfo.typeArguments.map(formatType).join(', ')}>`;
      }
      return typeInfo.name || 'unknown';
    case 'array':
      return `${formatType(typeInfo.elementType)}[]`;
    case 'union':
      return typeInfo.types?.map(formatType).join(' | ') || 'unknown';
    case 'intersection':
      return typeInfo.types?.map(formatType).join(' & ') || 'unknown';
    case 'literal':
      if (typeInfo.value === null) return 'null';
      if (typeof typeInfo.value === 'string') return `"${typeInfo.value}"`;
      return String(typeInfo.value);
    case 'tuple':
      return `[${(typeInfo.types || typeInfo.elements || []).map(formatType).join(', ')}]`;
    case 'reflection':
      return 'object';
    case 'typeOperator':
      return `${typeInfo.operator} ${formatType(typeInfo.target)}`;
    case 'predicate':
      return `${typeInfo.name} is ${formatType(typeInfo.targetType)}`;
    case 'conditional':
      return 'conditional';
    case 'mapped':
      return 'mapped';
    case 'indexedAccess':
      return `${formatType(typeInfo.objectType)}[${formatType(typeInfo.indexType)}]`;
    case 'query':
      return `typeof ${formatType(typeInfo.queryType)}`;
    default:
      return typeInfo.name || 'unknown';
  }
}

function formatSignatureString(name, sig) {
  if (!sig) return name;

  const params = sig.parameters || [];
  const typeParams = sig.typeParameter || [];

  const typeParamStr = typeParams.length > 0
    ? `<${typeParams.map(tp => tp.name).join(', ')}>`
    : '';

  const paramStr = params.map(p => {
    const optional = p.flags?.isOptional ? '?' : '';
    const rest = p.flags?.isRest ? '...' : '';
    const type = formatType(p.type);
    return `${rest}${p.name}${optional}: ${type}`;
  }).join(', ');

  const returnType = formatType(sig.type);
  return `${name}${typeParamStr}(${paramStr}): ${returnType}`;
}

function getKindLabel(kind) {
  switch (kind) {
    case ReflectionKind.Function: return 'function';
    case ReflectionKind.Class: return 'class';
    case ReflectionKind.Interface: return 'interface';
    case ReflectionKind.Variable: return 'variable';
    case ReflectionKind.Enum: return 'enum';
    case ReflectionKind.TypeAlias: return 'type';
    case ReflectionKind.Property: return 'function';
    case ReflectionKind.Method: return 'function';
    default: return 'function';
  }
}

// ---------------------------------------------------------------------------
// Generate full entry for an item (used in llms-full.txt)
// ---------------------------------------------------------------------------

function generateFullEntry(item, prefix) {
  const lines = [];
  const name = prefix ? `${prefix}.${item.name}` : item.name;
  const sig = getSignature(item);
  const desc = getDescription(item);
  const kind = getKindLabel(item.kind);

  lines.push(`### ${item.name}\n`);

  if (kind === 'class' || kind === 'interface') {
    lines.push(`\`${kind} ${item.name}\`\n`);
  } else if (sig) {
    lines.push(`\`${formatSignatureString(name, sig)}\`\n`);
  } else {
    lines.push(`\`${name}\`\n`);
  }

  if (desc) {
    lines.push(`${desc}\n`);
  }

  // Parameters (for functions)
  if (sig?.parameters?.length) {
    lines.push('**Parameters:**');
    for (const param of sig.parameters) {
      const type = formatType(param.type);
      const paramDesc = param.comment?.summary?.map(p => p.text).join('').trim() || '';
      const optional = param.flags?.isOptional ? ' (optional)' : '';
      const rest = param.flags?.isRest ? ' (rest)' : '';
      const descPart = paramDesc ? ` — ${paramDesc}` : '';
      lines.push(`- \`${param.name}\` (${type})${optional}${rest}${descPart}`);
    }
    lines.push('');
  }

  // Return type
  if (sig?.type) {
    lines.push(`**Returns:** \`${formatType(sig.type)}\`\n`);
  }

  // Class members
  if (item.kind === ReflectionKind.Class && item.children?.length) {
    const methods = item.children.filter(
      c => c.kind === ReflectionKind.Method && !c.flags?.isPrivate
    );
    const properties = item.children.filter(
      c => c.kind === ReflectionKind.Property && !c.flags?.isPrivate
    );

    if (properties.length > 0) {
      lines.push('**Properties:**');
      for (const prop of properties) {
        const type = formatType(prop.type);
        const propDesc = getDescription(prop);
        const descPart = propDesc ? ` — ${propDesc}` : '';
        lines.push(`- \`${prop.name}\`: \`${type}\`${descPart}`);
      }
      lines.push('');
    }

    if (methods.length > 0) {
      lines.push('**Methods:**');
      for (const method of methods) {
        const mSig = getSignature(method);
        if (mSig) {
          const mDesc = getDescription(method);
          const descPart = mDesc ? ` — ${mDesc}` : '';
          lines.push(`- \`${formatSignatureString(method.name, mSig)}\`${descPart}`);
        }
      }
      lines.push('');
    }
  }

  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// Generate llms.txt (concise overview)
// ---------------------------------------------------------------------------

function generateLlmsTxt() {
  const baseUrl = 'https://numwasm.dev';
  const lines = [];

  lines.push('# numwasm');
  lines.push('');
  lines.push('> NumPy-inspired n-dimensional array operations in TypeScript with WebAssembly acceleration.');
  lines.push('');
  lines.push('numwasm provides a comprehensive NumPy-compatible API for n-dimensional arrays');
  lines.push('in TypeScript, with performance-critical operations compiled to WebAssembly.');
  lines.push('');

  // Modules
  lines.push('## Modules');
  lines.push('');
  for (const mod of MODULE_DEFS) {
    const children = moduleChildrenMap.get(mod.slug) || [];
    const names = children.slice(0, 8).map(c => c.name);
    const suffix = children.length > 8 ? ', ...' : '';
    lines.push(`- [${mod.displayName}](${baseUrl}/docs/${mod.slug}): ${names.join(', ')}${suffix}`);
  }
  lines.push('');

  // Core API categories
  lines.push('## Core API');
  lines.push('');
  const categories = getCategorisedTopLevelItems();
  for (const cat of categories) {
    const names = cat.items.slice(0, 8).map(c => c.name);
    const suffix = cat.items.length > 8 ? ', ...' : '';
    lines.push(`- [${cat.title}](${baseUrl}/docs): ${names.join(', ')}${suffix}`);
  }
  lines.push('');

  // Links
  lines.push('## Optional');
  lines.push('');
  lines.push(`- [Full API Reference](${baseUrl}/llms-full.txt)`);
  lines.push(`- [Documentation](${baseUrl}/docs)`);
  lines.push(`- [Benchmarks](${baseUrl}/benchmarks)`);
  lines.push('');

  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// Generate llms-full.txt (complete API reference)
// ---------------------------------------------------------------------------

function generateLlmsFullTxt() {
  const lines = [];

  lines.push('# numwasm API Reference');
  lines.push('');
  lines.push('> Complete API reference for numwasm — NumPy-inspired n-dimensional array operations in TypeScript with WebAssembly acceleration.');
  lines.push('');

  // Modules
  for (const mod of MODULE_DEFS) {
    const children = moduleChildrenMap.get(mod.slug) || [];
    if (children.length === 0) continue;

    lines.push(`## ${mod.displayName} (nw.${mod.slug})`);
    lines.push('');

    for (const child of children) {
      lines.push(generateFullEntry(child, mod.slug));
    }
  }

  // Core API categories
  const categories = getCategorisedTopLevelItems();
  for (const cat of categories) {
    lines.push(`## ${cat.title}`);
    lines.push('');

    for (const item of cat.items) {
      lines.push(generateFullEntry(item, null));
    }
  }

  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

function main() {
  if (!fs.existsSync(apiJsonPath)) {
    console.error('Error: api.json not found at', apiJsonPath);
    console.error('Run typedoc first to generate api.json.');
    process.exit(1);
  }

  fs.mkdirSync(distDir, { recursive: true });

  const llmsTxt = generateLlmsTxt();
  const llmsFullTxt = generateLlmsFullTxt();

  fs.writeFileSync(path.join(distDir, 'llms.txt'), llmsTxt);
  fs.writeFileSync(path.join(distDir, 'llms-full.txt'), llmsFullTxt);

  const llmsSize = (Buffer.byteLength(llmsTxt) / 1024).toFixed(1);
  const fullSize = (Buffer.byteLength(llmsFullTxt) / 1024).toFixed(1);

  console.log(`Generated llms.txt (${llmsSize} KB)`);
  console.log(`Generated llms-full.txt (${fullSize} KB)`);
}

main();
