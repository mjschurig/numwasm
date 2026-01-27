import { useState } from 'react';
import { useParams, Link } from 'react-router-dom';
import { Button } from 'react-aria-components';
import { Menu, X } from 'lucide-react';
import { useApiDocs } from '../hooks/useApiDocs';
import DocsSidebar from '../components/docs/DocsSidebar';
import ModuleView from '../components/docs/ModuleView';
import ModuleOverviewPage from '../components/docs/ModuleOverviewPage';
import FunctionView from '../components/docs/FunctionView';
import ClassView from '../components/docs/ClassView';
import InterfaceView from '../components/docs/InterfaceView';
import { ReflectionKind, getKindString } from '../types/typedoc';
import type { DeclarationReflection } from '../types/typedoc';
import TypeView from '../components/docs/TypeView';
import CommentView from '../components/docs/CommentView';
import {
  MODULE_SLUGS,
  getModuleDef,
  getModuleChild,
  getTopLevelItem,
  getCategorisedTopLevelItems,
  getItemDescription,
} from '../utils/apiData';
import SEO from '../components/SEO';

// ---------------------------------------------------------------------------
// Docs Overview (landing page at /docs)
// ---------------------------------------------------------------------------

function DocsOverview() {
  const categories = getCategorisedTopLevelItems();

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold mb-2">API Documentation</h1>
        <p className="text-gray-400">
          NumPy-inspired n-dimensional array operations in TypeScript/WebAssembly
        </p>
      </div>

      <div className="space-y-8">
        {categories.map(cat => (
          <div key={cat.title}>
            <h2 className="text-xl font-semibold text-primary mb-4">
              {cat.title} ({cat.items.length})
            </h2>
            <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
              {cat.items.map(c => (
                <li key={c.id}>
                  <Link to={`/docs/${c.name}`} className="text-primary hover:text-primary/80">
                    {c.name}
                  </Link>
                </li>
              ))}
            </ul>
          </div>
        ))}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Single item view (function, class, interface, etc.)
// ---------------------------------------------------------------------------

function Breadcrumbs({ moduleName }: { moduleName?: string }) {
  const mod = moduleName ? getModuleDef(moduleName) : null;
  return (
    <nav className="text-sm text-gray-400 mb-4" aria-label="Breadcrumb">
      <ol className="flex items-center gap-1">
        <li>
          <Link to="/docs" className="hover:text-primary">Docs</Link>
        </li>
        {mod && (
          <li className="before:content-['/'] before:mx-1">
            <Link to={`/docs/${mod.slug}`} className="hover:text-primary">
              {mod.displayName}
            </Link>
          </li>
        )}
      </ol>
    </nav>
  );
}

function ItemView({ item, moduleName }: { item: DeclarationReflection; moduleName?: string }) {
  switch (item.kind) {
    case ReflectionKind.Function:
      return (
        <div className="text-gray-100">
          <Breadcrumbs moduleName={moduleName} />
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">Function</span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <FunctionView reflection={item} />
        </div>
      );

    case ReflectionKind.Class:
      return (
        <div>
          <Breadcrumbs moduleName={moduleName} />
          <ClassView reflection={item} />
        </div>
      );

    case ReflectionKind.Interface:
      return (
        <div>
          <Breadcrumbs moduleName={moduleName} />
          <InterfaceView reflection={item} />
        </div>
      );

    case ReflectionKind.Module:
    case ReflectionKind.Namespace:
      return <ModuleView reflection={item} />;

    case ReflectionKind.TypeAlias:
      return (
        <div className="text-gray-100">
          <Breadcrumbs moduleName={moduleName} />
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
          <Breadcrumbs moduleName={moduleName} />
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
          <Breadcrumbs moduleName={moduleName} />
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
          <Breadcrumbs moduleName={moduleName} />
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

// ---------------------------------------------------------------------------
// Not found
// ---------------------------------------------------------------------------

function NotFound({ name }: { name: string }) {
  return (
    <div className="text-gray-100">
      <Breadcrumbs />
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold">Not Found</h1>
      </div>
      <p className="text-gray-400">
        No documentation found for <code className="bg-gray-800 px-2 py-1 rounded">{name}</code>.
      </p>
      <Link to="/docs" className="text-primary hover:text-primary/80 mt-4 inline-block">
        Back to API Documentation
      </Link>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main Docs page
// ---------------------------------------------------------------------------

export default function Docs() {
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const { '*': wildcard } = useParams();
  const { data } = useApiDocs();

  // Parse the wildcard path: "" | "itemName" | "moduleName/itemName"
  const pathParts = (wildcard || '').split('/').filter(Boolean);

  let content: React.ReactNode;
  let seoNode: React.ReactNode = (
    <SEO
      title="API Documentation"
      description="Complete API reference for numwasm — NumPy-inspired n-dimensional array operations in TypeScript/WebAssembly."
      path="/docs"
      breadcrumbs={[
        { name: 'Home', url: '/' },
        { name: 'Docs', url: '/docs' },
      ]}
    />
  );

  if (pathParts.length === 0) {
    // /docs — overview
    content = <DocsOverview />;
  } else if (pathParts.length === 1) {
    const name = pathParts[0];
    if (MODULE_SLUGS.has(name)) {
      const mod = getModuleDef(name)!;
      content = <ModuleOverviewPage module={mod} />;
      seoNode = (
        <SEO
          title={`${mod.displayName} — API`}
          description={`numwasm ${mod.displayName} module reference. Functions and classes for ${mod.displayName.toLowerCase()} operations.`}
          path={`/docs/${mod.slug}`}
          breadcrumbs={[
            { name: 'Home', url: '/' },
            { name: 'Docs', url: '/docs' },
            { name: mod.displayName, url: `/docs/${mod.slug}` },
          ]}
        />
      );
    } else {
      const item = getTopLevelItem(name);
      if (item) {
        content = <ItemView item={item} />;
        seoNode = (
          <SEO
            title={`${item.name} — API`}
            description={getItemDescription(item)}
            path={`/docs/${item.name}`}
            breadcrumbs={[
              { name: 'Home', url: '/' },
              { name: 'Docs', url: '/docs' },
              { name: item.name, url: `/docs/${item.name}` },
            ]}
          />
        );
      } else {
        content = <NotFound name={name} />;
      }
    }
  } else if (pathParts.length === 2) {
    const [moduleName, itemName] = pathParts;
    if (MODULE_SLUGS.has(moduleName)) {
      const mod = getModuleDef(moduleName)!;
      const item = getModuleChild(moduleName, itemName);
      if (item) {
        content = <ItemView item={item} moduleName={moduleName} />;
        seoNode = (
          <SEO
            title={`${item.name} — ${mod.displayName} — API`}
            description={getItemDescription(item)}
            path={`/docs/${moduleName}/${item.name}`}
            breadcrumbs={[
              { name: 'Home', url: '/' },
              { name: 'Docs', url: '/docs' },
              { name: mod.displayName, url: `/docs/${mod.slug}` },
              { name: item.name, url: `/docs/${moduleName}/${item.name}` },
            ]}
          />
        );
      } else {
        content = <NotFound name={`${moduleName}/${itemName}`} />;
      }
    } else {
      content = <NotFound name={pathParts.join('/')} />;
    }
  } else {
    content = <NotFound name={pathParts.join('/')} />;
  }

  return (
    <div className="flex min-h-[calc(100vh-4rem)]">
      {seoNode}
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
        {content}
      </main>
    </div>
  );
}
