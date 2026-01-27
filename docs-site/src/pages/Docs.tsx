import { useState } from "react";
import { useParams, Link } from "react-router-dom";
import { Button } from "react-aria-components";
import { Menu, X } from "lucide-react";
import { useApiDocs } from "../hooks/useApiDocs";
import DocsSidebar from "../components/docs/DocsSidebar";
import ModuleView from "../components/docs/ModuleView";
import ModuleOverviewPage from "../components/docs/ModuleOverviewPage";
import FunctionView from "../components/docs/FunctionView";
import ClassView from "../components/docs/ClassView";
import InterfaceView from "../components/docs/InterfaceView";
import { ReflectionKind, getKindString } from "../types/typedoc";
import type { DeclarationReflection } from "../types/typedoc";
import TypeView from "../components/docs/TypeView";
import CommentView from "../components/docs/CommentView";
import {
  PACKAGES,
  PACKAGE_IDS,
  getPackageModuleSlugs,
  getModuleDef,
  getModuleChild,
  getTopLevelItem,
  getCategorisedTopLevelItems,
  getItemDescription,
  getPackageModuleDefs,
  hasPackageData,
} from "../utils/apiData";
import type { PackageId } from "../utils/apiData";
import SEO from "../components/SEO";

// ---------------------------------------------------------------------------
// Docs Landing (package picker at /docs)
// ---------------------------------------------------------------------------

function DocsLanding() {
  return (
    <div className="text-gray-100 max-w-4xl mx-auto py-8">
      <SEO
        title="Documentation — *wasm"
        description="API documentation for numwasm, sciwasm, and symwasm — scientific computing for TypeScript."
        path="/docs"
      />
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold mb-2">Documentation</h1>
        <p className="text-gray-400">
          Choose a package to browse its API reference.
        </p>
      </div>

      <div className="grid grid-cols-1 sm:grid-cols-3 gap-6">
        {PACKAGES.map((pkg) => (
          <Link
            key={pkg.id}
            to={`/docs/${pkg.id}`}
            className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors no-underline"
          >
            <h3 className="font-semibold text-primary mb-3 text-lg">
              {pkg.displayName}
            </h3>
            <p className="text-gray-400 text-sm">{pkg.description}</p>
            <div className="mt-4 text-xs text-gray-500">
              {getPackageModuleDefs(pkg.id).length} modules
              {hasPackageData(pkg.id) ? "" : " (coming soon)"}
            </div>
          </Link>
        ))}
      </div>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Package Overview (landing for /docs/:packageId with no further path)
// ---------------------------------------------------------------------------

function PackageOverview({ packageId }: { packageId: PackageId }) {
  const pkg = PACKAGES.find((p) => p.id === packageId)!;
  const moduleDefs = getPackageModuleDefs(packageId);
  const categories = getCategorisedTopLevelItems(packageId);
  const hasData = hasPackageData(packageId);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold mb-2">{pkg.displayName} API</h1>
        <p className="text-gray-400">{pkg.description}</p>
      </div>

      {!hasData && (
        <div className="bg-yellow-900/20 border border-yellow-700/50 rounded-lg p-4 mb-8">
          <p className="text-yellow-300 text-sm">
            This package is under development. Module stubs are scaffolded but
            full API documentation is not yet available.
          </p>
        </div>
      )}

      <div className="mb-8">
        <h2 className="text-xl font-semibold text-primary mb-4">Modules</h2>
        <div className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-3">
          {moduleDefs.map((mod) => (
            <Link
              key={mod.slug}
              to={`/docs/${packageId}/${mod.slug}`}
              className="p-4 bg-gray-800/30 border border-gray-700/50 rounded-lg hover:border-primary/30 transition-colors no-underline"
            >
              <span className="text-primary font-medium">
                {mod.displayName}
              </span>
            </Link>
          ))}
        </div>
      </div>

      {categories.length > 0 && (
        <div className="space-y-8">
          {categories.map((cat) => (
            <div key={cat.title}>
              <h2 className="text-xl font-semibold text-primary mb-4">
                {cat.title} ({cat.items.length})
              </h2>
              <ul className="grid grid-cols-[repeat(auto-fill,minmax(200px,1fr))] gap-2">
                {cat.items.map((c) => (
                  <li key={c.id}>
                    <Link
                      to={`/docs/${packageId}/${c.name}`}
                      className="text-primary hover:text-primary/80"
                    >
                      {c.name}
                    </Link>
                  </li>
                ))}
              </ul>
            </div>
          ))}
        </div>
      )}
    </div>
  );
}

// ---------------------------------------------------------------------------
// Single item view (function, class, interface, etc.)
// ---------------------------------------------------------------------------

function Breadcrumbs({
  packageId,
  moduleName,
}: {
  packageId: PackageId;
  moduleName?: string;
}) {
  const mod = moduleName ? getModuleDef(moduleName, packageId) : null;
  const pkg = PACKAGES.find((p) => p.id === packageId)!;
  return (
    <nav className="text-sm text-gray-400 mb-4" aria-label="Breadcrumb">
      <ol className="flex items-center gap-1">
        <li>
          <Link to="/docs" className="hover:text-primary">
            Docs
          </Link>
        </li>
        <li className="before:content-['/'] before:mx-1">
          <Link to={`/docs/${packageId}`} className="hover:text-primary">
            {pkg.displayName}
          </Link>
        </li>
        {mod && (
          <li className="before:content-['/'] before:mx-1">
            <Link
              to={`/docs/${packageId}/${mod.slug}`}
              className="hover:text-primary"
            >
              {mod.displayName}
            </Link>
          </li>
        )}
      </ol>
    </nav>
  );
}

function ItemView({
  item,
  packageId,
  moduleName,
}: {
  item: DeclarationReflection;
  packageId: PackageId;
  moduleName?: string;
}) {
  switch (item.kind) {
    case ReflectionKind.Function:
    case ReflectionKind.Method:
    case ReflectionKind.Property: {
      const isCallable =
        item.kind === ReflectionKind.Function ||
        item.kind === ReflectionKind.Method ||
        !!(item.type as any)?.declaration?.signatures;
      if (isCallable) {
        return (
          <div className="text-gray-100">
            <Breadcrumbs packageId={packageId} moduleName={moduleName} />
            <div className="border-b border-gray-700 pb-4 mb-8">
              <span className="text-gray-400 text-sm uppercase tracking-wide">
                Function
              </span>
              <h1 className="text-3xl font-bold">{item.name}</h1>
            </div>
            <FunctionView reflection={item} />
          </div>
        );
      }
      return (
        <div className="text-gray-100">
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              Property
            </span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <CommentView comment={item.comment} />
        </div>
      );
    }

    case ReflectionKind.Class:
      return (
        <div>
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <ClassView reflection={item} />
        </div>
      );

    case ReflectionKind.Interface:
      return (
        <div>
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <InterfaceView reflection={item} />
        </div>
      );

    case ReflectionKind.Module:
    case ReflectionKind.Namespace:
      return <ModuleView reflection={item} />;

    case ReflectionKind.TypeAlias:
      return (
        <div className="text-gray-100">
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              Type Alias
            </span>
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
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              Enum
            </span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <CommentView comment={item.comment} />
          {item.children && (
            <table className="w-full mt-4">
              <thead>
                <tr className="border-b border-gray-700">
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">
                    Member
                  </th>
                  <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase">
                    Value
                  </th>
                </tr>
              </thead>
              <tbody>
                {item.children.map((member) => (
                  <tr key={member.id} className="border-b border-gray-800">
                    <td className="py-3 px-4 font-mono text-primary">
                      {member.name}
                    </td>
                    <td className="py-3 px-4 font-mono text-gray-300">
                      {member.defaultValue || "-"}
                    </td>
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
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
          <div className="border-b border-gray-700 pb-4 mb-8">
            <span className="text-gray-400 text-sm uppercase tracking-wide">
              {item.flags?.isConst ? "Constant" : "Variable"}
            </span>
            <h1 className="text-3xl font-bold">{item.name}</h1>
          </div>
          <div className="bg-gray-800/50 border border-gray-700/50 p-4 rounded-lg font-mono mb-4 overflow-x-auto">
            {item.flags?.isConst && "const "}
            {item.name}: <TypeView type={item.type} />
            {item.defaultValue && ` = ${item.defaultValue}`}
          </div>
          <CommentView comment={item.comment} />
        </div>
      );

    default:
      return (
        <div className="text-gray-100">
          <Breadcrumbs packageId={packageId} moduleName={moduleName} />
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

function NotFound({
  name,
  packageId,
}: {
  name: string;
  packageId?: PackageId;
}) {
  return (
    <div className="text-gray-100">
      {packageId && <Breadcrumbs packageId={packageId} />}
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold">Not Found</h1>
      </div>
      <p className="text-gray-400">
        No documentation found for{" "}
        <code className="bg-gray-800 px-2 py-1 rounded">{name}</code>.
      </p>
      <Link
        to="/docs"
        className="text-primary hover:text-primary/80 mt-4 inline-block"
      >
        Back to Documentation
      </Link>
    </div>
  );
}

// ---------------------------------------------------------------------------
// Main Docs page
// ---------------------------------------------------------------------------

export default function Docs() {
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const { "*": wildcard } = useParams();

  // Parse the wildcard: "" | "packageId" | "packageId/moduleOrItem" | "packageId/module/item"
  const pathParts = (wildcard || "").split("/").filter(Boolean);

  // Determine if first segment is a package ID
  const firstSegment = pathParts[0] || "";
  const isPackageRoute = PACKAGE_IDS.has(firstSegment);
  const packageId: PackageId | null = isPackageRoute
    ? (firstSegment as PackageId)
    : null;
  const subParts = isPackageRoute ? pathParts.slice(1) : pathParts;

  const { data } = useApiDocs(packageId || "numwasm");

  let content: React.ReactNode;
  let seoNode: React.ReactNode = (
    <SEO
      title="Documentation — *wasm"
      description="API documentation for numwasm, sciwasm, and symwasm — scientific computing for TypeScript."
      path="/docs"
      breadcrumbs={[
        { name: "Home", url: "/" },
        { name: "Docs", url: "/docs" },
      ]}
    />
  );

  if (!packageId && pathParts.length === 0) {
    // /docs — landing page with package picker
    content = <DocsLanding />;
  } else if (packageId && subParts.length === 0) {
    // /docs/:packageId — package overview
    const pkg = PACKAGES.find((p) => p.id === packageId)!;
    content = <PackageOverview packageId={packageId} />;
    seoNode = (
      <SEO
        title={`${pkg.displayName} — API Documentation`}
        description={`API reference for ${pkg.displayName} — ${pkg.description}.`}
        path={`/docs/${packageId}`}
        breadcrumbs={[
          { name: "Home", url: "/" },
          { name: "Docs", url: "/docs" },
          { name: pkg.displayName, url: `/docs/${packageId}` },
        ]}
      />
    );
  } else if (packageId && subParts.length === 1) {
    const name = subParts[0];
    const moduleSlugs = getPackageModuleSlugs(packageId);
    if (moduleSlugs.has(name)) {
      const mod = getModuleDef(name, packageId)!;
      content = <ModuleOverviewPage module={mod} packageId={packageId} />;
      seoNode = (
        <SEO
          title={`${mod.displayName} — ${packageId} API`}
          description={`${packageId} ${mod.displayName} module reference.`}
          path={`/docs/${packageId}/${mod.slug}`}
          breadcrumbs={[
            { name: "Home", url: "/" },
            { name: "Docs", url: "/docs" },
            { name: packageId, url: `/docs/${packageId}` },
            { name: mod.displayName, url: `/docs/${packageId}/${mod.slug}` },
          ]}
        />
      );
    } else {
      const item = getTopLevelItem(name, packageId);
      if (item) {
        content = <ItemView item={item} packageId={packageId} />;
        seoNode = (
          <SEO
            title={`${item.name} — ${packageId} API`}
            description={getItemDescription(item, packageId)}
            path={`/docs/${packageId}/${item.name}`}
            breadcrumbs={[
              { name: "Home", url: "/" },
              { name: "Docs", url: "/docs" },
              { name: packageId, url: `/docs/${packageId}` },
              { name: item.name, url: `/docs/${packageId}/${item.name}` },
            ]}
          />
        );
      } else {
        content = <NotFound name={name} packageId={packageId} />;
      }
    }
  } else if (packageId && subParts.length === 2) {
    const [moduleName, itemName] = subParts;
    const moduleSlugs = getPackageModuleSlugs(packageId);
    if (moduleSlugs.has(moduleName)) {
      const mod = getModuleDef(moduleName, packageId)!;
      const item = getModuleChild(moduleName, itemName, packageId);
      if (item) {
        content = (
          <ItemView item={item} packageId={packageId} moduleName={moduleName} />
        );
        seoNode = (
          <SEO
            title={`${item.name} — ${mod.displayName} — ${packageId} API`}
            description={getItemDescription(item, packageId)}
            path={`/docs/${packageId}/${moduleName}/${item.name}`}
            breadcrumbs={[
              { name: "Home", url: "/" },
              { name: "Docs", url: "/docs" },
              { name: packageId, url: `/docs/${packageId}` },
              { name: mod.displayName, url: `/docs/${packageId}/${mod.slug}` },
              {
                name: item.name,
                url: `/docs/${packageId}/${moduleName}/${item.name}`,
              },
            ]}
          />
        );
      } else {
        content = (
          <NotFound name={`${moduleName}/${itemName}`} packageId={packageId} />
        );
      }
    } else {
      content = <NotFound name={subParts.join("/")} packageId={packageId} />;
    }
  } else {
    content = (
      <NotFound name={pathParts.join("/")} packageId={packageId || undefined} />
    );
  }

  // Show sidebar only when inside a specific package
  const showSidebar = !!packageId;

  return (
    <div className="flex min-h-[calc(100vh-4rem)]">
      {seoNode}
      {showSidebar && (
        <>
          {/* Mobile sidebar toggle */}
          <Button
            onPress={() => setSidebarOpen(!sidebarOpen)}
            className="lg:hidden fixed bottom-4 right-4 z-50 p-3 bg-primary text-white rounded-full shadow-lg cursor-pointer hover:bg-primary/90 transition-colors"
            aria-label="Toggle sidebar"
          >
            {sidebarOpen ? (
              <X className="w-6 h-6" />
            ) : (
              <Menu className="w-6 h-6" />
            )}
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
              ${sidebarOpen ? "translate-x-0" : "-translate-x-full lg:translate-x-0"}
            `}
          >
            <DocsSidebar
              data={data}
              packageId={packageId}
              onItemClick={() => setSidebarOpen(false)}
            />
          </div>
        </>
      )}

      <main className="flex-1 p-4 sm:p-8 max-w-full overflow-x-hidden">
        {content}
      </main>
    </div>
  );
}
