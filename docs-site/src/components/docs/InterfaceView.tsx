import type { DeclarationReflection } from '../../types/typedoc';
import { ReflectionKind } from '../../types/typedoc';
import TypeView from './TypeView';
import CommentView from './CommentView';

interface InterfaceViewProps {
  reflection: DeclarationReflection;
}

export default function InterfaceView({ reflection }: InterfaceViewProps) {
  const children = reflection.children || [];

  const properties = children.filter(c => c.kind === ReflectionKind.Property);
  const methods = children.filter(c => c.kind === ReflectionKind.Method);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h2 className="text-3xl font-bold">{reflection.name}</h2>
        {reflection.extendedTypes && reflection.extendedTypes.length > 0 && (
          <p className="text-gray-400 mt-1">
            extends{' '}
            {reflection.extendedTypes.map((t, i) => (
              <span key={i}>
                <TypeView type={t} />
                {i < reflection.extendedTypes!.length - 1 ? ', ' : ''}
              </span>
            ))}
          </p>
        )}
      </div>

      <CommentView comment={reflection.comment} />

      {properties.length > 0 && (
        <div className="mb-8">
          <h3 className="text-xl font-semibold text-primary mb-4">Properties</h3>
          <table className="w-full">
            <thead>
              <tr className="border-b border-gray-700">
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Name</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
              </tr>
            </thead>
            <tbody>
              {properties.map(prop => (
                <tr key={prop.id} className="border-b border-gray-800">
                  <td className="py-3 px-4 font-mono text-primary">
                    {prop.flags?.isReadonly && <span className="text-gray-500">readonly </span>}
                    {prop.name}
                    {prop.flags?.isOptional && '?'}
                  </td>
                  <td className="py-3 px-4"><TypeView type={prop.type} /></td>
                  <td className="py-3 px-4 text-gray-400">
                    {prop.comment?.summary?.map(p => p.text).join('') || '-'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {methods.length > 0 && (
        <div className="mb-8">
          <h3 className="text-xl font-semibold text-primary mb-4">Methods</h3>
          {methods.map(method => {
            const sig = method.signatures?.[0];
            if (!sig) return null;

            const params = sig.parameters || [];
            return (
              <div key={method.id} className="mb-4">
                <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono text-sm overflow-x-auto text-gray-200">
                  {method.name}({params.map((p, i) => (
                    <span key={p.id}>
                      {p.name}{p.flags?.isOptional && '?'}: <TypeView type={p.type} />
                      {i < params.length - 1 && ', '}
                    </span>
                  ))}): <TypeView type={sig.type} />
                </div>
                <p className="mt-2 text-gray-400">
                  {sig.comment?.summary?.map(p => p.text).join('') || ''}
                </p>
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}
