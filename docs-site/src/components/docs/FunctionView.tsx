import type { DeclarationReflection } from '../../types/typedoc';
import TypeView from './TypeView';
import CommentView from './CommentView';

interface FunctionViewProps {
  reflection: DeclarationReflection;
}

export default function FunctionView({ reflection }: FunctionViewProps) {
  const signature = reflection.signatures?.[0];

  if (!signature) {
    return (
      <div className="mb-8">
        <h3 className="text-lg font-semibold text-gray-100">{reflection.name}</h3>
        <p className="text-gray-400">No signature available</p>
      </div>
    );
  }

  const params = signature.parameters || [];
  const typeParams = signature.typeParameter || [];

  // Build signature string
  const typeParamStr = typeParams.length > 0
    ? `<${typeParams.map(tp => tp.name).join(', ')}>`
    : '';

  const paramStr = params.map(p => {
    const optional = p.flags?.isOptional ? '?' : '';
    const rest = p.flags?.isRest ? '...' : '';
    return `${rest}${p.name}${optional}`;
  }).join(', ');

  return (
    <div className="mb-8">
      <h3 className="text-lg font-semibold mb-2 text-gray-100">{reflection.name}</h3>
      <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono text-sm overflow-x-auto mb-4 text-gray-200">
        {reflection.name}{typeParamStr}({paramStr}): <TypeView type={signature.type} />
      </div>

      <CommentView comment={signature.comment} />

      {params.length > 0 && (
        <div className="mt-4">
          <h4 className="text-sm font-semibold text-gray-400 mb-2">Parameters</h4>
          <table className="w-full">
            <thead>
              <tr className="border-b border-gray-700">
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Name</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
              </tr>
            </thead>
            <tbody>
              {params.map(param => (
                <tr key={param.id} className="border-b border-gray-800">
                  <td className="py-3 px-4">
                    <span className="font-mono text-primary">
                      {param.flags?.isRest && '...'}
                      {param.name}
                      {param.flags?.isOptional && '?'}
                    </span>
                  </td>
                  <td className="py-3 px-4"><TypeView type={param.type} /></td>
                  <td className="py-3 px-4 text-gray-400">
                    {param.comment?.summary?.map(p => p.text).join('') || '-'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}
    </div>
  );
}
