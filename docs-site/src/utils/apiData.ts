import type { ProjectReflection, DeclarationReflection } from '../types/typedoc';
import { ReflectionKind } from '../types/typedoc';

// ---------------------------------------------------------------------------
// Package definitions
// ---------------------------------------------------------------------------

export type PackageId =
  | 'numwasm'
  | 'sciwasm'
  | 'symwasm'
  | 'arwasm'
  | 'lawasm'
  | 'linwasm'
  | 'quadwasm'
  | 'superluwasm'
  | 'xsfwasm'
  | 'odewasm';

export const PACKAGES: { id: PackageId; displayName: string; description: string }[] = [
  { id: 'numwasm', displayName: 'numwasm', description: 'NumPy-compatible n-dimensional arrays' },
  { id: 'sciwasm', displayName: 'sciwasm', description: 'SciPy-compatible scientific computing' },
  { id: 'symwasm', displayName: 'symwasm', description: 'SymPy-compatible symbolic math' },
  { id: 'arwasm', displayName: 'arwasm', description: 'ARPACK eigenvalue solvers' },
  { id: 'lawasm', displayName: 'lawasm', description: 'LAPACK linear algebra routines' },
  { id: 'linwasm', displayName: 'linwasm', description: 'LINPACK legacy linear algebra' },
  { id: 'quadwasm', displayName: 'quadwasm', description: 'QUADPACK numerical integration' },
  { id: 'superluwasm', displayName: 'superluwasm', description: 'SuperLU sparse matrix solvers' },
  { id: 'xsfwasm', displayName: 'xsfwasm', description: 'Special functions (Bessel, Gamma, etc.)' },
  { id: 'odewasm', displayName: 'odewasm', description: 'ODE solvers' },
];

export const PACKAGE_IDS = new Set<string>(PACKAGES.map(p => p.id));

// ---------------------------------------------------------------------------
// Module definitions (per package)
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

const NUMWASM_MODULE_DEFS: ModuleDef[] = [
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

const SCIWASM_MODULE_DEFS: ModuleDef[] = [
  { slug: 'optimize', displayName: 'Optimization', varName: 'optimize', isNamespace: true },
  { slug: 'integrate', displayName: 'Integration', varName: 'integrate', isNamespace: true },
  { slug: 'interpolate', displayName: 'Interpolation', varName: 'interpolate', isNamespace: true },
  { slug: 'stats', displayName: 'Statistics', varName: 'stats', isNamespace: true },
  { slug: 'signal', displayName: 'Signal Processing', varName: 'signal', isNamespace: true },
  { slug: 'spatial', displayName: 'Spatial', varName: 'spatial', isNamespace: true },
  { slug: 'special', displayName: 'Special Functions', varName: 'special', isNamespace: true },
  { slug: 'sparse', displayName: 'Sparse Matrices', varName: 'sparse', isNamespace: true },
  { slug: 'cluster', displayName: 'Clustering', varName: 'cluster', isNamespace: true },
  { slug: 'io', displayName: 'I/O', varName: 'io', isNamespace: true },
  { slug: 'ndimage', displayName: 'N-D Image', varName: 'ndimage', isNamespace: true },
  { slug: 'constants', displayName: 'Constants', varName: 'constants', isNamespace: true },
];

const SYMWASM_MODULE_DEFS: ModuleDef[] = [
  { slug: 'core', displayName: 'Core', varName: 'core', isNamespace: true },
  { slug: 'simplify', displayName: 'Simplify', varName: 'simplify', isNamespace: true },
  { slug: 'solvers', displayName: 'Solvers', varName: 'solvers', isNamespace: true },
  { slug: 'calculus', displayName: 'Calculus', varName: 'calculus', isNamespace: true },
  { slug: 'matrices', displayName: 'Matrices', varName: 'matrices', isNamespace: true },
  { slug: 'printing', displayName: 'Printing', varName: 'printing', isNamespace: true },
];

// Standalone WASM packages (lower-level, fewer TS wrappers)
const ARWASM_MODULE_DEFS: ModuleDef[] = [];
const LAWASM_MODULE_DEFS: ModuleDef[] = [];
const LINWASM_MODULE_DEFS: ModuleDef[] = [];
const QUADWASM_MODULE_DEFS: ModuleDef[] = [];
const SUPERLUWASM_MODULE_DEFS: ModuleDef[] = [];
const XSFWASM_MODULE_DEFS: ModuleDef[] = [];
const ODEWASM_MODULE_DEFS: ModuleDef[] = [];

const PACKAGE_MODULE_DEFS: Record<PackageId, ModuleDef[]> = {
  numwasm: NUMWASM_MODULE_DEFS,
  sciwasm: SCIWASM_MODULE_DEFS,
  symwasm: SYMWASM_MODULE_DEFS,
  arwasm: ARWASM_MODULE_DEFS,
  lawasm: LAWASM_MODULE_DEFS,
  linwasm: LINWASM_MODULE_DEFS,
  quadwasm: QUADWASM_MODULE_DEFS,
  superluwasm: SUPERLUWASM_MODULE_DEFS,
  xsfwasm: XSFWASM_MODULE_DEFS,
  odewasm: ODEWASM_MODULE_DEFS,
};

// ---------------------------------------------------------------------------
// Category definitions for top-level (non-module) items — numwasm only
// ---------------------------------------------------------------------------

export interface CategoryDef {
  title: string;
  names: string[];
}

const NUMWASM_CATEGORY_DEFS: CategoryDef[] = [
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

const PACKAGE_CATEGORY_DEFS: Record<PackageId, CategoryDef[]> = {
  numwasm: NUMWASM_CATEGORY_DEFS,
  sciwasm: [],
  symwasm: [],
  arwasm: [],
  lawasm: [],
  linwasm: [],
  quadwasm: [],
  superluwasm: [],
  xsfwasm: [],
  odewasm: [],
};

// ---------------------------------------------------------------------------
// Per-package data store
// ---------------------------------------------------------------------------

interface PackageData {
  apiData: ProjectReflection | null;
  allChildren: DeclarationReflection[];
  moduleClaimedNames: Set<string>;
  moduleChildrenMap: Map<string, DeclarationReflection[]>;
  categoryClaimedNames: Set<string>;
  referenceNames: Set<string>;
}

function loadPackageJson(packageId: PackageId): ProjectReflection | null {
  try {
    // Vite glob import with eager loading for all api JSON files
    const modules = import.meta.glob('../../public/api-*.json', { eager: true });
    const key = `../../public/api-${packageId}.json`;
    const mod = modules[key] as { default: unknown } | undefined;
    return mod ? (mod.default as unknown as ProjectReflection) : null;
  } catch {
    return null;
  }
}

function buildPackageData(packageId: PackageId): PackageData {
  const apiData = loadPackageJson(packageId);
  const allChildren: DeclarationReflection[] = apiData?.children || [];
  const moduleDefs = PACKAGE_MODULE_DEFS[packageId];
  const categoryDefs = PACKAGE_CATEGORY_DEFS[packageId];

  const moduleClaimedNames = new Set<string>();
  const moduleChildrenMap = new Map<string, DeclarationReflection[]>();

  for (const mod of moduleDefs) {
    let children: DeclarationReflection[] = [];

    if (mod.isNamespace) {
      const ns = allChildren.find(
        c => c.kind === ReflectionKind.Namespace && c.name === mod.slug
      );
      if (ns?.children) {
        children = ns.children;
      }
      moduleClaimedNames.add(mod.slug);
    } else if (mod.varName) {
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

  const categoryClaimedNames = new Set<string>();
  for (const cat of categoryDefs) {
    for (const n of cat.names) categoryClaimedNames.add(n);
  }

  const referenceNames = new Set(
    allChildren.filter(c => c.kind === ReflectionKind.Reference).map(c => c.name)
  );

  return { apiData, allChildren, moduleClaimedNames, moduleChildrenMap, categoryClaimedNames, referenceNames };
}

// Build data for all packages at import time
const packageDataMap = new Map<PackageId, PackageData>();
for (const pkg of PACKAGES) {
  packageDataMap.set(pkg.id, buildPackageData(pkg.id));
}

// ---------------------------------------------------------------------------
// Backward-compatible exports (numwasm is the default)
// ---------------------------------------------------------------------------

export const apiData: ProjectReflection = packageDataMap.get('numwasm')!.apiData || ({ children: [] } as unknown as ProjectReflection);
export const MODULE_DEFS = NUMWASM_MODULE_DEFS;
export const MODULE_SLUGS = new Set(NUMWASM_MODULE_DEFS.map(m => m.slug));
export const CATEGORY_DEFS = NUMWASM_CATEGORY_DEFS;

// ---------------------------------------------------------------------------
// Per-package public API
// ---------------------------------------------------------------------------

export function getPackageModuleDefs(packageId: PackageId): ModuleDef[] {
  return PACKAGE_MODULE_DEFS[packageId];
}

export function getPackageModuleSlugs(packageId: PackageId): Set<string> {
  return new Set(PACKAGE_MODULE_DEFS[packageId].map(m => m.slug));
}

export function getPackageCategoryDefs(packageId: PackageId): CategoryDef[] {
  return PACKAGE_CATEGORY_DEFS[packageId];
}

export function getPackageApiData(packageId: PackageId): ProjectReflection | null {
  return packageDataMap.get(packageId)?.apiData || null;
}

export function hasPackageData(packageId: PackageId): boolean {
  const data = packageDataMap.get(packageId);
  return !!(data?.apiData && data.allChildren.length > 0);
}

// ---------------------------------------------------------------------------
// Public API (package-aware versions)
// ---------------------------------------------------------------------------

export function getModuleDef(slug: string, packageId: PackageId = 'numwasm'): ModuleDef | undefined {
  return PACKAGE_MODULE_DEFS[packageId].find(m => m.slug === slug);
}

export function getModuleChildren(slug: string, packageId: PackageId = 'numwasm'): DeclarationReflection[] {
  const data = packageDataMap.get(packageId);
  return data?.moduleChildrenMap.get(slug) || [];
}

export function getModuleChild(
  moduleSlug: string,
  itemName: string,
  packageId: PackageId = 'numwasm'
): DeclarationReflection | undefined {
  const children = getModuleChildren(moduleSlug, packageId);
  return children.find(c => c.name === itemName);
}

export function getTopLevelItem(name: string, packageId: PackageId = 'numwasm'): DeclarationReflection | undefined {
  const data = packageDataMap.get(packageId);
  if (!data) return undefined;
  return data.allChildren.find(
    c => c.name === name && c.kind !== ReflectionKind.Reference
  );
}

export function resolveItemByPath(path: string, packageId: PackageId = 'numwasm'): DeclarationReflection | undefined {
  const parts = path.split('/');
  if (parts.length === 2) {
    return getModuleChild(parts[0], parts[1], packageId);
  }
  return getTopLevelItem(parts[0], packageId);
}

export function getCategorisedTopLevelItems(packageId: PackageId = 'numwasm'): Array<{
  title: string;
  items: DeclarationReflection[];
}> {
  const data = packageDataMap.get(packageId);
  if (!data) return [];

  const categoryDefs = PACKAGE_CATEGORY_DEFS[packageId];
  const { allChildren, moduleClaimedNames, categoryClaimedNames, referenceNames } = data;
  const result: Array<{ title: string; items: DeclarationReflection[] }> = [];

  for (const cat of categoryDefs) {
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
export function getItemDescription(item: DeclarationReflection, packageId: PackageId = 'numwasm'): string {
  const summary =
    item.comment?.summary
    || item.signatures?.[0]?.comment?.summary
    || (item.type as any)?.declaration?.signatures?.[0]?.comment?.summary;
  if (!summary) return `${packageId} ${item.name}`;
  return summary
    .map((p: any) => p.text)
    .join('')
    .replace(/\n/g, ' ')
    .slice(0, 160);
}

/** Generate all doc route paths (for sitemap + SSG) */
export function getAllRoutes(): string[] {
  const routes: string[] = ['/docs'];

  for (const pkg of PACKAGES) {
    const packageId = pkg.id;
    const moduleDefs = PACKAGE_MODULE_DEFS[packageId];
    const data = packageDataMap.get(packageId);

    // Package landing page
    routes.push(`/docs/${packageId}`);

    // Module overview pages
    for (const mod of moduleDefs) {
      routes.push(`/docs/${packageId}/${mod.slug}`);
      // Module children
      const children = getModuleChildren(mod.slug, packageId);
      for (const child of children) {
        routes.push(`/docs/${packageId}/${mod.slug}/${child.name}`);
      }
    }

    // Top-level items (non-module, non-reference) — only for packages with data
    if (data && data.allChildren.length > 0) {
      for (const child of data.allChildren) {
        if (
          child.kind !== ReflectionKind.Reference &&
          child.kind !== ReflectionKind.Namespace &&
          !data.moduleClaimedNames.has(child.name)
        ) {
          routes.push(`/docs/${packageId}/${child.name}`);
        }
      }
    }
  }

  return routes;
}
