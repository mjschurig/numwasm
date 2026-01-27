import apiJson from '../../public/api.json';
import type { ProjectReflection, DeclarationReflection } from '../types/typedoc';
import { ReflectionKind } from '../types/typedoc';

export const apiData: ProjectReflection = apiJson as unknown as ProjectReflection;

// ---------------------------------------------------------------------------
// Module definitions
// ---------------------------------------------------------------------------

export interface ModuleDef {
  slug: string;
  displayName: string;
  /** Name of the Variable (kind=32) in api.json that holds the module children,
   *  or null if items are identified by name patterns at the top level. */
  varName: string | null;
  /** For namespace kind=4 entries instead of variable kind=32 */
  isNamespace?: boolean;
  /** Function/class name patterns to match at top level (when varName is null) */
  namePatterns?: RegExp[];
  /** Explicit names to include */
  names?: string[];
}

export const MODULE_DEFS: ModuleDef[] = [
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

export const MODULE_SLUGS = new Set(MODULE_DEFS.map(m => m.slug));

// ---------------------------------------------------------------------------
// Category definitions for top-level (non-module) items
// ---------------------------------------------------------------------------

export interface CategoryDef {
  title: string;
  names: string[];
}

export const CATEGORY_DEFS: CategoryDef[] = [
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
    names: [
      'union1d', 'intersect1d', 'setdiff1d', 'setxor1d', 'in1d', 'isin',
    ],
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
// Derived lookups (built once at import time)
// ---------------------------------------------------------------------------

const allChildren: DeclarationReflection[] = apiData.children || [];

/** Map from item id to item for quick lookup */
const idMap = new Map<number, DeclarationReflection>();
function buildIdMap(items: DeclarationReflection[]) {
  for (const item of items) {
    idMap.set(item.id, item);
    if (item.children) buildIdMap(item.children);
  }
}
buildIdMap(allChildren);

/** Set of names claimed by any module */
const moduleClaimedNames = new Set<string>();

/** Map from module slug to its children */
const moduleChildrenMap = new Map<string, DeclarationReflection[]>();

for (const mod of MODULE_DEFS) {
  let children: DeclarationReflection[] = [];

  if (mod.isNamespace) {
    // Find the Namespace (kind=4) entry
    const ns = allChildren.find(
      c => c.kind === ReflectionKind.Namespace && c.name === mod.slug
    );
    if (ns?.children) {
      children = ns.children;
    }
    moduleClaimedNames.add(mod.slug);
  } else if (mod.varName) {
    // Find the Variable (kind=32) whose type.declaration.children holds the module
    const varItem = allChildren.find(
      c => c.kind === ReflectionKind.Variable && c.name === mod.varName
    );
    if (varItem) {
      const decl = (varItem.type as any)?.declaration;
      if (decl?.children) {
        children = decl.children;
      }
      moduleClaimedNames.add(mod.varName);
    }
  }

  // Also match top-level items by explicit names or patterns
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

/** Set of names claimed by a category */
const categoryClaimedNames = new Set<string>();
for (const cat of CATEGORY_DEFS) {
  for (const n of cat.names) categoryClaimedNames.add(n);
}

// Reference kind items (re-exports) – these should not appear in listings
const referenceNames = new Set(
  allChildren.filter(c => c.kind === ReflectionKind.Reference).map(c => c.name)
);

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

export function getModuleDef(slug: string): ModuleDef | undefined {
  return MODULE_DEFS.find(m => m.slug === slug);
}

export function getModuleChildren(slug: string): DeclarationReflection[] {
  return moduleChildrenMap.get(slug) || [];
}

/** Get a module child by slug and item name */
export function getModuleChild(
  moduleSlug: string,
  itemName: string
): DeclarationReflection | undefined {
  const children = getModuleChildren(moduleSlug);
  return children.find(c => c.name === itemName);
}

/** Get a top-level item (not inside a module) by name */
export function getTopLevelItem(name: string): DeclarationReflection | undefined {
  return allChildren.find(
    c => c.name === name && c.kind !== ReflectionKind.Reference
  );
}

/** Resolve item by path: "linalg/matmul" or "abs" */
export function resolveItemByPath(path: string): DeclarationReflection | undefined {
  const parts = path.split('/');
  if (parts.length === 2) {
    return getModuleChild(parts[0], parts[1]);
  }
  return getTopLevelItem(parts[0]);
}

/**
 * Get categorised top-level items (not claimed by any module).
 * Returns the CATEGORY_DEFS groups plus an "Other" group for uncategorised items.
 */
export function getCategorisedTopLevelItems(): Array<{
  title: string;
  items: DeclarationReflection[];
}> {
  const result: Array<{ title: string; items: DeclarationReflection[] }> = [];

  for (const cat of CATEGORY_DEFS) {
    const nameSet = new Set(cat.names);
    const items = allChildren.filter(
      c =>
        nameSet.has(c.name) &&
        c.kind !== ReflectionKind.Reference &&
        !moduleClaimedNames.has(c.name)
    );
    if (items.length > 0) {
      result.push({ title: cat.title, items });
    }
  }

  // Types group
  const typeItems = allChildren.filter(
    c =>
      (c.kind === ReflectionKind.Interface ||
        c.kind === ReflectionKind.TypeAlias ||
        c.kind === ReflectionKind.Enum) &&
      !moduleClaimedNames.has(c.name)
  );
  if (typeItems.length > 0) {
    result.push({ title: 'Types & Interfaces', items: typeItems });
  }

  // Constants group
  const constantNames = new Set([
    'e', 'pi', 'inf', 'nan', 'euler_gamma', 'newaxis', 'ellipsis',
    'NAN', 'NINF', 'NZERO', 'PINF', 'PZERO',
  ]);
  const constants = allChildren.filter(
    c =>
      c.kind === ReflectionKind.Variable &&
      constantNames.has(c.name) &&
      !moduleClaimedNames.has(c.name)
  );
  if (constants.length > 0) {
    result.push({ title: 'Constants', items: constants });
  }

  // Other (uncategorised)
  const allClaimedNames = new Set([
    ...moduleClaimedNames,
    ...categoryClaimedNames,
    ...constantNames,
    ...Array.from(referenceNames),
  ]);
  // Also claim type items
  for (const t of typeItems) allClaimedNames.add(t.name);

  const other = allChildren.filter(
    c =>
      !allClaimedNames.has(c.name) &&
      c.kind !== ReflectionKind.Reference &&
      c.kind !== ReflectionKind.Namespace
  );
  if (other.length > 0) {
    result.push({ title: 'Other', items: other });
  }

  return result;
}

/** Extract a short description from a declaration's comment/signature */
export function getItemDescription(item: DeclarationReflection): string {
  const summary =
    item.comment?.summary
    || item.signatures?.[0]?.comment?.summary
    || (item.type as any)?.declaration?.signatures?.[0]?.comment?.summary;
  if (!summary) return `numwasm ${item.name}`;
  return summary
    .map((p: any) => p.text)
    .join('')
    .replace(/\n/g, ' ')
    .slice(0, 160);
}

/** Generate all doc route paths (for sitemap + SSG) */
export function getAllRoutes(): string[] {
  const routes: string[] = ['/docs'];

  // Module overview pages
  for (const mod of MODULE_DEFS) {
    routes.push(`/docs/${mod.slug}`);
    // Module children
    const children = getModuleChildren(mod.slug);
    for (const child of children) {
      routes.push(`/docs/${mod.slug}/${child.name}`);
    }
  }

  // Top-level items (non-module, non-reference)
  for (const child of allChildren) {
    if (
      child.kind !== ReflectionKind.Reference &&
      child.kind !== ReflectionKind.Namespace &&
      !moduleClaimedNames.has(child.name)
    ) {
      routes.push(`/docs/${child.name}`);
    }
  }

  return routes;
}
