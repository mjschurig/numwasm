import type { DeclarationReflection } from '../../types/typedoc';
import { ReflectionKind, getKindString } from '../../types/typedoc';
import FunctionView from './FunctionView';
import ClassView from './ClassView';
import InterfaceView from './InterfaceView';
import TypeView from './TypeView';
import CommentView from './CommentView';

interface ModuleViewProps {
  reflection: DeclarationReflection;
}

export default function ModuleView({ reflection }: ModuleViewProps) {
  const children = reflection.children || [];

  const functions = children.filter(c => c.kind === ReflectionKind.Function);
  const classes = children.filter(c => c.kind === ReflectionKind.Class);
  const interfaces = children.filter(c => c.kind === ReflectionKind.Interface);
  const typeAliases = children.filter(c => c.kind === ReflectionKind.TypeAlias);
  const variables = children.filter(c => c.kind === ReflectionKind.Variable);
  const enums = children.filter(c => c.kind === ReflectionKind.Enum);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <span className="text-gray-400 text-sm uppercase tracking-wide">
          {getKindString(reflection.kind)}
        </span>
        <h1 className="text-3xl font-bold">{reflection.name}</h1>
      </div>

      <CommentView comment={reflection.comment} />

      {classes.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Classes</h2>
          {classes.map(cls => (
            <div key={cls.id} className="mb-12">
              <ClassView reflection={cls} />
            </div>
          ))}
        </div>
      )}

      {interfaces.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Interfaces</h2>
          {interfaces.map(iface => (
            <div key={iface.id} className="mb-8">
              <InterfaceView reflection={iface} />
            </div>
          ))}
        </div>
      )}

      {typeAliases.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Type Aliases</h2>
          <table className="w-full">
            <thead>
              <tr className="border-b border-gray-700">
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Name</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
              </tr>
            </thead>
            <tbody>
              {typeAliases.map(alias => (
                <tr key={alias.id} className="border-b border-gray-800">
                  <td className="py-3 px-4 font-mono text-primary">{alias.name}</td>
                  <td className="py-3 px-4"><TypeView type={alias.type} /></td>
                  <td className="py-3 px-4 text-gray-400">
                    {alias.comment?.summary?.map(p => p.text).join('') || '-'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {enums.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Enums</h2>
          {enums.map(en => (
            <div key={en.id} className="mb-8">
              <h3 className="text-lg font-semibold mb-2">{en.name}</h3>
              <CommentView comment={en.comment} />
              {en.children && (
                <table className="w-full mt-4">
                  <thead>
                    <tr className="border-b border-gray-700">
                      <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Member</th>
                      <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Value</th>
                    </tr>
                  </thead>
                  <tbody>
                    {en.children.map(member => (
                      <tr key={member.id} className="border-b border-gray-800">
                        <td className="py-3 px-4 font-mono text-primary">{member.name}</td>
                        <td className="py-3 px-4 font-mono text-gray-300">{member.defaultValue || '-'}</td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              )}
            </div>
          ))}
        </div>
      )}

      {variables.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Variables</h2>
          <table className="w-full">
            <thead>
              <tr className="border-b border-gray-700">
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Name</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Type</th>
                <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Description</th>
              </tr>
            </thead>
            <tbody>
              {variables.map(v => (
                <tr key={v.id} className="border-b border-gray-800">
                  <td className="py-3 px-4 font-mono text-primary">
                    {v.flags?.isConst && <span className="text-gray-500">const </span>}
                    {v.name}
                  </td>
                  <td className="py-3 px-4"><TypeView type={v.type} /></td>
                  <td className="py-3 px-4 text-gray-400">
                    {v.comment?.summary?.map(p => p.text).join('') || '-'}
                  </td>
                </tr>
              ))}
            </tbody>
          </table>
        </div>
      )}

      {functions.length > 0 && (
        <div className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">Functions</h2>
          {functions.map(fn => (
            <FunctionView key={fn.id} reflection={fn} />
          ))}
        </div>
      )}
    </div>
  );
}
