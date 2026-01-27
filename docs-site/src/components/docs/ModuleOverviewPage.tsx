import { Link } from 'react-router-dom';
import type { DeclarationReflection } from '../../types/typedoc';
import { ReflectionKind } from '../../types/typedoc';
import type { ModuleDef, PackageId } from '../../utils/apiData';
import { getModuleChildren, getItemDescription } from '../../utils/apiData';

interface ModuleOverviewPageProps {
  module: ModuleDef;
  packageId?: PackageId;
}

export default function ModuleOverviewPage({ module, packageId = 'numwasm' }: ModuleOverviewPageProps) {
  const children = getModuleChildren(module.slug, packageId);

  const functions = children.filter(
    c => c.kind === ReflectionKind.Function || c.kind === ReflectionKind.Method
  );
  const classes = children.filter(c => c.kind === ReflectionKind.Class);
  const interfaces = children.filter(c => c.kind === ReflectionKind.Interface);
  const types = children.filter(c => c.kind === ReflectionKind.TypeAlias);
  const variables = children.filter(
    c => c.kind === ReflectionKind.Variable || c.kind === ReflectionKind.Enum
  );
  // Items with other kinds (e.g. properties inside module-like objects rendered as functions)
  const other = children.filter(
    c =>
      c.kind !== ReflectionKind.Function &&
      c.kind !== ReflectionKind.Method &&
      c.kind !== ReflectionKind.Class &&
      c.kind !== ReflectionKind.Interface &&
      c.kind !== ReflectionKind.TypeAlias &&
      c.kind !== ReflectionKind.Variable &&
      c.kind !== ReflectionKind.Enum
  );

  const sections: Array<{ title: string; items: DeclarationReflection[] }> = [];
  if (classes.length > 0) sections.push({ title: 'Classes', items: classes });
  if (interfaces.length > 0) sections.push({ title: 'Interfaces', items: interfaces });
  if (types.length > 0) sections.push({ title: 'Type Aliases', items: types });
  if (variables.length > 0) sections.push({ title: 'Constants & Variables', items: variables });
  if (functions.length > 0) sections.push({ title: 'Functions', items: functions });
  if (other.length > 0) sections.push({ title: 'Other', items: other });

  return (
    <div className="text-gray-100">
      {/* Breadcrumbs */}
      <nav className="text-sm text-gray-400 mb-4" aria-label="Breadcrumb">
        <ol className="flex items-center gap-1">
          <li><Link to="/docs" className="hover:text-primary">Docs</Link></li>
          <li className="before:content-['/'] before:mx-1">
            <Link to={`/docs/${packageId}`} className="hover:text-primary">{packageId}</Link>
          </li>
          <li className="before:content-['/'] before:mx-1">{module.displayName}</li>
        </ol>
      </nav>

      <div className="border-b border-gray-700 pb-4 mb-8">
        <span className="text-gray-400 text-sm uppercase tracking-wide">Module</span>
        <h1 className="text-3xl font-bold">{module.displayName}</h1>
        <p className="text-gray-400 mt-2">
          {children.length} item{children.length !== 1 ? 's' : ''}
        </p>
      </div>

      {sections.map(section => (
        <div key={section.title} className="mb-8">
          <h2 className="text-xl font-semibold text-primary mb-4">{section.title}</h2>
          <div className="border border-gray-700/50 rounded-lg overflow-hidden">
            <table className="w-full">
              <thead>
                <tr className="border-b border-gray-700 bg-gray-800/30">
                  <th className="py-2 px-4 text-left text-sm font-semibold text-gray-400 uppercase w-48">
                    Name
                  </th>
                  <th className="py-2 px-4 text-left text-sm font-semibold text-gray-400 uppercase">
                    Description
                  </th>
                </tr>
              </thead>
              <tbody>
                {section.items.map(item => (
                  <tr
                    key={item.id || item.name}
                    id={`fn-${item.name}`}
                    className="border-b border-gray-800 hover:bg-gray-800/30 transition-colors"
                  >
                    <td className="py-2 px-4 font-mono">
                      <Link
                        to={`/docs/${packageId}/${module.slug}/${item.name}`}
                        className="text-primary hover:text-primary/80"
                      >
                        {item.name}
                      </Link>
                    </td>
                    <td className="py-2 px-4 text-gray-400 text-sm">
                      {getItemDescription(item, packageId)}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
          </div>
        </div>
      ))}
    </div>
  );
}
