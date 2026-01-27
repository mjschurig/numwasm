/**
 * Generate docs-index.json for the symwasm MCP server.
 * Reads api-symwasm.json and produces a structured, searchable JSON index.
 *
 * Usage: node generate-docs-index.js
 * Run after: typedoc (so api-symwasm.json exists)
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const outputPath = path.resolve(__dirname, 'docs-index.json');
const apiJsonPath = path.resolve(__dirname, '../../docs-site/public/api-symwasm.json');

// ---------------------------------------------------------------------------
// ReflectionKind constants
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
// Module definitions for symwasm
// ---------------------------------------------------------------------------

const MODULE_DEFS = [
  { slug: 'core', displayName: 'Core Expressions', varName: 'core' },
  { slug: 'simplify', displayName: 'Simplification', varName: 'simplify' },
  { slug: 'solvers', displayName: 'Equation Solvers', varName: 'solvers' },
  { slug: 'calculus', displayName: 'Calculus', varName: 'calculus' },
  { slug: 'matrices', displayName: 'Matrices', varName: 'matrices' },
  { slug: 'printing', displayName: 'Printing', varName: 'printing' },
];

// ---------------------------------------------------------------------------
// Helpers
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
    case 'intrinsic': return typeInfo.name || 'unknown';
    case 'reference':
      if (typeInfo.typeArguments?.length)
        return `${typeInfo.name}<${typeInfo.typeArguments.map(formatType).join(', ')}>`;
      return typeInfo.name || 'unknown';
    case 'array': return `${formatType(typeInfo.elementType)}[]`;
    case 'union': return typeInfo.types?.map(formatType).join(' | ') || 'unknown';
    case 'intersection': return typeInfo.types?.map(formatType).join(' & ') || 'unknown';
    case 'literal':
      if (typeInfo.value === null) return 'null';
      if (typeof typeInfo.value === 'string') return `"${typeInfo.value}"`;
      return String(typeInfo.value);
    case 'tuple': return `[${(typeInfo.types || typeInfo.elements || []).map(formatType).join(', ')}]`;
    case 'reflection': return 'object';
    default: return typeInfo.name || 'unknown';
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

function buildSection(item, moduleSlug, category) {
  const qualifiedName = moduleSlug ? `${moduleSlug}.${item.name}` : item.name;
  const sig = getSignature(item);
  const desc = getDescription(item);
  const kind = getKindLabel(item.kind);

  let signature = '';
  if (kind === 'class' || kind === 'interface') {
    signature = `${kind} ${item.name}`;
  } else if (sig) {
    signature = formatSignatureString(qualifiedName, sig);
  } else {
    signature = qualifiedName;
  }

  const lines = [];
  lines.push(`### ${item.name}\n`);
  lines.push(`\`${signature}\`\n`);
  if (desc) lines.push(`${desc}\n`);

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

  if (sig?.type) {
    lines.push(`**Returns:** \`${formatType(sig.type)}\`\n`);
  }

  return {
    name: item.name,
    module: moduleSlug || null,
    category: category || null,
    signature,
    description: desc,
    content: lines.join('\n'),
  };
}

function buildOverview() {
  const lines = [];
  lines.push('# symwasm');
  lines.push('');
  lines.push('> SymPy-inspired symbolic mathematics in TypeScript.');
  lines.push('');
  lines.push('symwasm provides symbolic computation capabilities including symbolic expressions,');
  lines.push('algebraic simplification, equation solving, calculus operations,');
  lines.push('matrix algebra, and expression formatting/printing.');
  lines.push('');
  lines.push('## Modules');
  lines.push('');
  for (const mod of MODULE_DEFS) {
    lines.push(`- ${mod.displayName} (${mod.slug})`);
  }
  lines.push('');
  return lines.join('\n');
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------

function main() {
  if (!fs.existsSync(apiJsonPath)) {
    console.warn('Warning: api-symwasm.json not found at', apiJsonPath);
    console.warn('Writing empty docs-index.json. Run typedoc first to generate full index.');
    fs.writeFileSync(outputPath, JSON.stringify({ overview: buildOverview(), sections: [] }, null, 2));
    return;
  }

  const apiData = JSON.parse(fs.readFileSync(apiJsonPath, 'utf-8'));
  const allChildren = apiData.children || [];

  const sections = [];

  for (const mod of MODULE_DEFS) {
    const varItem = allChildren.find(c => c.kind === ReflectionKind.Variable && c.name === mod.varName);
    if (varItem) {
      const decl = varItem.type?.declaration;
      if (decl?.children) {
        for (const child of decl.children) {
          sections.push(buildSection(child, mod.slug, mod.displayName));
        }
      }
    }

    const ns = allChildren.find(c => c.kind === ReflectionKind.Namespace && c.name === mod.slug);
    if (ns?.children) {
      for (const child of ns.children) {
        sections.push(buildSection(child, mod.slug, mod.displayName));
      }
    }
  }

  for (const child of allChildren) {
    if (child.kind === ReflectionKind.Reference) continue;
    const isModule = MODULE_DEFS.some(m => m.varName === child.name || m.slug === child.name);
    if (!isModule) {
      sections.push(buildSection(child, null, null));
    }
  }

  const index = {
    overview: buildOverview(),
    sections,
  };

  fs.writeFileSync(outputPath, JSON.stringify(index, null, 2));

  const sizeKB = (Buffer.byteLength(JSON.stringify(index)) / 1024).toFixed(1);
  console.log(`Generated docs-index.json: ${sections.length} sections (${sizeKB} KB)`);
}

main();
