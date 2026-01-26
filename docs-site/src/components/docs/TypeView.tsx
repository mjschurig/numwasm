import type { TypeInfo } from '../../types/typedoc';

interface TypeViewProps {
  type: TypeInfo | undefined;
}

export default function TypeView({ type }: TypeViewProps) {
  if (!type) {
    return <span className="param-type">any</span>;
  }

  const renderType = (t: TypeInfo): string => {
    switch (t.type) {
      case 'intrinsic':
        return t.name || 'unknown';

      case 'literal':
        if (typeof t.value === 'string') {
          return `"${t.value}"`;
        }
        return String(t.value);

      case 'reference':
        let result = t.name || 'unknown';
        if (t.typeArguments && t.typeArguments.length > 0) {
          result += `<${t.typeArguments.map(renderType).join(', ')}>`;
        }
        return result;

      case 'array':
        return `${renderType(t.elementType!)}[]`;

      case 'union':
        return t.types?.map(renderType).join(' | ') || 'unknown';

      case 'intersection':
        return t.types?.map(renderType).join(' & ') || 'unknown';

      case 'tuple':
        return `[${t.types?.map(renderType).join(', ') || ''}]`;

      case 'reflection':
        if (t.declaration?.signatures) {
          const sig = t.declaration.signatures[0];
          const params = sig.parameters?.map(p => `${p.name}: ${renderType(p.type!)}`).join(', ') || '';
          const returnType = renderType(sig.type!);
          return `(${params}) => ${returnType}`;
        }
        if (t.declaration?.children) {
          const props = t.declaration.children.map(c => `${c.name}: ${renderType(c.type!)}`).join('; ');
          return `{ ${props} }`;
        }
        return 'object';

      case 'indexedAccess':
        return `${renderType(t.types![0])}[${renderType(t.types![1])}]`;

      case 'conditional':
        return 'conditional';

      case 'mapped':
        return 'mapped';

      case 'templateLiteral':
        return 'template literal';

      case 'query':
        return `typeof ${t.name}`;

      default:
        return t.name || t.type || 'unknown';
    }
  };

  return <span className="param-type">{renderType(type)}</span>;
}
