import type { DeclarationReflection } from '../../types/typedoc';
import { ReflectionKind } from '../../types/typedoc';
import TypeView from './TypeView';
import CommentView from './CommentView';

interface ClassViewProps {
  reflection: DeclarationReflection;
}

export default function ClassView({ reflection }: ClassViewProps) {
  const children = reflection.children || [];

  const constructor = children.find(c => c.kind === ReflectionKind.Constructor);
  const properties = children.filter(c => c.kind === ReflectionKind.Property);
  const methods = children.filter(c => c.kind === ReflectionKind.Method);
  const accessors = children.filter(c => c.kind === ReflectionKind.Accessor);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h2 className="text-3xl font-bold">{reflection.name}</h2>
        {reflection.extendedTypes && reflection.extendedTypes.length > 0 && (
          <p className="text-gray-400 mt-1">
            extends <TypeView type={reflection.extendedTypes[0]} />
          </p>
        )}
      </div>

      <CommentView comment={reflection.comment} />

      {constructor && constructor.signatures?.[0] && (
        <div className="mb-8">
          <h3 className="text-xl font-semibold text-primary mb-4">Constructor</h3>
          <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono text-sm overflow-x-auto text-gray-200">
            new {reflection.name}(
            {constructor.signatures[0].parameters?.map((p, i, arr) => (
              <span key={p.id}>
                {p.name}{p.flags?.isOptional ? '?' : ''}: <TypeView type={p.type} />
                {i < arr.length - 1 ? ', ' : ''}
              </span>
            ))}
            )
          </div>
          {constructor.signatures[0].parameters && constructor.signatures[0].parameters.length > 0 && (
            <table className="w-full mt-4">
              <thead>
                <tr className="border-b border-gray-700">
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Parameter</th>
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
                </tr>
              </thead>
              <tbody>
                {constructor.signatures[0].parameters.map(param => (
                  <tr key={param.id} className="border-b border-gray-800">
                    <td className="py-3 px-4 font-mono text-primary">{param.name}</td>
                    <td className="py-3 px-4"><TypeView type={param.type} /></td>
                    <td className="py-3 px-4 text-gray-400">
                      {param.comment?.summary?.map(p => p.text).join('') || '-'}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </div>
      )}

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

      {accessors.length > 0 && (
        <div className="mb-8">
          <h3 className="text-xl font-semibold text-primary mb-4">Accessors</h3>
          <table className="w-full">
            <thead>
              <tr className="border-b border-gray-700">
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Name</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
              </tr>
            </thead>
            <tbody>
              {accessors.map(acc => (
                <tr key={acc.id} className="border-b border-gray-800">
                  <td className="py-3 px-4 font-mono text-primary">{acc.name}</td>
                  <td className="py-3 px-4"><TypeView type={acc.getSignature?.type || acc.setSignature?.type} /></td>
                  <td className="py-3 px-4 text-gray-400">
                    {acc.comment?.summary?.map(p => p.text).join('') || '-'}
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
              <div key={method.id} className="mb-6 pb-6 border-b border-gray-800 last:border-0">
                <h4 className="font-semibold mb-2">{method.name}</h4>
                <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono text-sm overflow-x-auto mb-2 text-gray-200">
                  {method.name}({params.map((p, i) => (
                    <span key={p.id}>
                      {p.flags?.isRest && '...'}
                      {p.name}
                      {p.flags?.isOptional && '?'}
                      : <TypeView type={p.type} />
                      {i < params.length - 1 && ', '}
                    </span>
                  ))}): <TypeView type={sig.type} />
                </div>
                <CommentView comment={sig.comment} />
              </div>
            );
          })}
        </div>
      )}
    </div>
  );
}
