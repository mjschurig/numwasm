import { useState } from 'react';
import { useParams, Link } from 'react-router-dom';
import { Button } from 'react-aria-components';
import { Menu, X } from 'lucide-react';
import { useApiDocs } from '../hooks/useApiDocs';
import DocsSidebar from '../components/docs/DocsSidebar';
import ModuleView from '../components/docs/ModuleView';
import FunctionView from '../components/docs/FunctionView';
import ClassView from '../components/docs/ClassView';
import InterfaceView from '../components/docs/InterfaceView';
import { ReflectionKind, getKindString } from '../types/typedoc';
import type { DeclarationReflection } from '../types/typedoc';
import TypeView from '../components/docs/TypeView';
import CommentView from '../components/docs/CommentView';

function DocsOverview({ children }: { children: DeclarationReflection[] }) {
  const functions = children.filter(c => c.kind === ReflectionKind.Function);
  const classes = children.filter(c => c.kind === ReflectionKind.Class);
  const interfaces = children.filter(c => c.kind === ReflectionKind.Interface);
  const types = children.filter(c => c.kind === ReflectionKind.TypeAlias);
  const variables = children.filter(c => c.kind === ReflectionKind.Variable);
  const enums = children.filter(c => c.kind === ReflectionKind.Enum);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold mb-2">API Documentation</h1>
        <p className="text-gray-400">
          NumPy-inspired n-dimensional array operations in TypeScript/WebAssembly
        </p>
      </div>

      <div className="space-y-8">
        {classes.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Classes ({classes.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {classes.map(c => (
                <li key={c.id}>
                  <Link to={`/docs/${c.name}`} className="text-primary hover:text-primary-dark">{c.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}

        {interfaces.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Interfaces ({interfaces.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {interfaces.map(i => (
                <li key={i.id}>
                  <Link to={`/docs/${i.name}`} className="text-primary hover:text-primary-dark">{i.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}

        {types.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Type Aliases ({types.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {types.map(t => (
                <li key={t.id}>
                  <Link to={`/docs/${t.name}`} className="text-primary hover:text-primary-dark">{t.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}

        {enums.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Enums ({enums.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {enums.map(e => (
                <li key={e.id}>
                  <Link to={`/docs/${e.name}`} className="text-primary hover:text-primary-dark">{e.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}

        {variables.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Variables ({variables.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {variables.map(v => (
                <li key={v.id}>
                  <Link to={`/docs/${v.name}`} className="text-primary hover:text-primary-dark">{v.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}

        {functions.length > 0 && (
          <div>
            <h2 className="text-xl font-semibold text-primary mb-4">Functions ({functions.length})</h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {functions.map(f => (
                <li key={f.id}>
                  <Link to={`/docs/${f.name}`} className="text-primary hover:text-primary-dark">{f.name}</Link>
                </li>
              ))}
            </ul>
          </div>
        )}
      </div>
    </div>
  );
}

function ItemView({ item }: { item: DeclarationReflection }) {
  switch (item.kind) {
    case ReflectionKind.Function:
      return (
        <div className="text-gray-100">
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">Function</span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <FunctionView reflection={item} />
        </div>
      );

    case ReflectionKind.Class:
      return <ClassView reflection={item} />;

    case ReflectionKind.Interface:
      return <InterfaceView reflection={item} />;

    case ReflectionKind.Module:
    case ReflectionKind.Namespace:
      return <ModuleView reflection={item} />;

    case ReflectionKind.TypeAlias:
      return (
        <div className="text-gray-100">
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">Type Alias</span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono mb-4 overflow-x-auto">
            type {item.name} = <TypeView type={item.type} />
          </div>
          <CommentView comment={item.comment} />
        </div>
      );

    case ReflectionKind.Enum:
      return (
        <div className="text-gray-100">
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">Enum</span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <CommentView comment={item.comment} />
          {item.children && (
            <table className="w-full mt-4">
              <thead>
                <tr className="border-b border-gray-700">
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Member</th>
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">Value</th>
                </tr>
              </thead>
              <tbody>
                {item.children.map(member => (
                  <tr key={member.id} className="border-b border-gray-800">
                    <td className="py-3 px-4 font-mono text-primary">{member.name}</td>
                    <td className="py-3 px-4 font-mono text-gray-300">{member.defaultValue || '-'}</td>
                  </tr>
                ))}
              </tbody>
            </table>
          )}
        </div>
      );

    case ReflectionKind.Variable:
      return (
        <div className="text-gray-100">
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              {item.flags?.isConst ? 'Constant' : 'Variable'}
            </span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono mb-4 overflow-x-auto">
            {item.flags?.isConst && 'const '}{item.name}: <TypeView type={item.type} />
            {item.defaultValue && ` = ${item.defaultValue}`}
          </div>
          <CommentView comment={item.comment} />
        </div>
      );

    default:
      return (
        <div className="text-gray-100">
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              {getKindString(item.kind)}
            </span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <CommentView comment={item.comment} />
        </div>
      );
  }
}

export default function Docs() {
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const { '*': itemName } = useParams();
  const { data, loading, error } = useApiDocs();

  if (loading) {
    return (
      <div className="flex min-h-[calc(100vh-4rem)]">
        <div className="flex-1 p-4 sm:p-8">
          <div className="flex justify-center items-center min-h-[200px] text-gray-400">
            Loading documentation...
          </div>
        </div>
      </div>
    );
  }

  if (error || !data) {
    return (
      <div className="flex min-h-[calc(100vh-4rem)]">
        <div className="flex-1 p-4 sm:p-8">
          <div className="text-red-400 p-4 sm:p-8">
            <h2 className="text-xl font-bold mb-2">Error Loading Documentation</h2>
            <p>{error || 'Failed to load API documentation'}</p>
            <p className="mt-4 text-gray-400">
              Make sure to run <code className="bg-gray-800 px-2 py-1 rounded">pnpm run docs</code> to generate the API documentation.
            </p>
          </div>
        </div>
      </div>
    );
  }

  const children = data.children || [];
  const selectedItem = itemName
    ? children.find(c => c.name === itemName)
    : null;

  return (
    <div className="flex min-h-[calc(100vh-4rem)]">
      {/* Mobile sidebar toggle */}
      <Button
        onPress={() => setSidebarOpen(!sidebarOpen)}
        className="lg:hidden fixed bottom-4 right-4 z-50 p-3 bg-primary text-white rounded-full shadow-lg cursor-pointer hover:bg-primary/90 transition-colors"
        aria-label="Toggle sidebar"
      >
        {sidebarOpen ? <X className="w-6 h-6" /> : <Menu className="w-6 h-6" />}
      </Button>

      {/* Sidebar overlay for mobile */}
      {sidebarOpen && (
        <div
          className="lg:hidden fixed inset-0 bg-black/50 z-30"
          onClick={() => setSidebarOpen(false)}
        />
      )}

      {/* Sidebar */}
      <div
        className={`
          fixed lg:static inset-y-0 left-0 z-40 transform transition-transform duration-300 ease-in-out
          ${sidebarOpen ? 'translate-x-0' : '-translate-x-full lg:translate-x-0'}
        `}
      >
        <DocsSidebar data={data} onItemClick={() => setSidebarOpen(false)} />
      </div>

      <main className="flex-1 p-4 sm:p-8 max-w-full overflow-x-hidden">
        {selectedItem ? (
          <ItemView item={selectedItem} />
        ) : (
          <DocsOverview children={children} />
        )}
      </main>
    </div>
  );
}
