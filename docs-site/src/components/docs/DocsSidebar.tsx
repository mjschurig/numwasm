import { useState, useMemo } from 'react';
import { NavLink, useParams } from 'react-router-dom';
import { SearchField, Input, Button } from 'react-aria-components';
import type { ProjectReflection } from '../../types/typedoc';
import { ReflectionKind } from '../../types/typedoc';
import {
  MODULE_DEFS,
  getModuleChildren,
  getCategorisedTopLevelItems,
} from '../../utils/apiData';

interface DocsSidebarProps {
  data: ProjectReflection;
  onItemClick?: () => void;
}

interface SidebarGroup {
  title: string;
  /** If this is a module, its slug for linking to /docs/<slug> */
  moduleSlug?: string;
  items: Array<{ name: string; kind: number; linkTo: string }>;
}

export default function DocsSidebar({ data, onItemClick }: DocsSidebarProps) {
  const [search, setSearch] = useState('');
  const { '*': wildcard } = useParams();
  const pathParts = (wildcard || '').split('/').filter(Boolean);
  const activeModule = pathParts.length >= 1 ? pathParts[0] : null;

  // Auto-expand the active module plus some defaults
  const [expanded, setExpanded] = useState<Record<string, boolean>>(() => {
    const initial: Record<string, boolean> = {};
    if (activeModule) initial[activeModule] = true;
    return initial;
  });

  const groups = useMemo((): SidebarGroup[] => {
    const result: SidebarGroup[] = [];

    // 1. Modules section
    for (const mod of MODULE_DEFS) {
      const children = getModuleChildren(mod.slug);
      result.push({
        title: mod.displayName,
        moduleSlug: mod.slug,
        items: children.map(c => ({
          name: c.name,
          kind: c.kind,
          linkTo: `/docs/${mod.slug}/${c.name}`,
        })),
      });
    }

    // 2. Categorised top-level items
    const categories = getCategorisedTopLevelItems();
    for (const cat of categories) {
      result.push({
        title: cat.title,
        items: cat.items.map(c => ({
          name: c.name,
          kind: c.kind,
          linkTo: `/docs/${c.name}`,
        })),
      });
    }

    return result;
  }, [data]);

  const filteredGroups = useMemo(() => {
    if (!search.trim()) return groups;
    const searchLower = search.toLowerCase();
    return groups
      .map(group => ({
        ...group,
        items: group.items.filter(item =>
          item.name.toLowerCase().includes(searchLower)
        ),
      }))
      .filter(group => group.items.length > 0);
  }, [groups, search]);

  const toggleGroup = (key: string) => {
    setExpanded(prev => ({ ...prev, [key]: !prev[key] }));
  };

  const getKindIcon = (kind: number): string => {
    switch (kind) {
      case ReflectionKind.Class: return 'C';
      case ReflectionKind.Interface: return 'I';
      case ReflectionKind.Function: return 'F';
      case ReflectionKind.Variable: return 'V';
      case ReflectionKind.TypeAlias: return 'T';
      case ReflectionKind.Enum: return 'E';
      case ReflectionKind.Namespace: return 'N';
      case ReflectionKind.Module: return 'M';
      case ReflectionKind.Property: return 'P';
      case ReflectionKind.Method: return 'M';
      default: return '';
    }
  };

  const groupKey = (group: SidebarGroup) => group.moduleSlug || group.title;

  return (
    <aside className="w-72 lg:w-70 border-r border-gray-700 p-4 sm:p-6 overflow-y-auto h-screen lg:sticky lg:top-16 lg:h-[calc(100vh-4rem)] bg-gray-900/95 lg:bg-gray-900/50 backdrop-blur-sm">
      <SearchField value={search} onChange={setSearch} className="mb-4">
        <Input
          placeholder="Search..."
          className="w-full px-3 py-2 border border-gray-700 bg-gray-800/50 text-gray-200 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-primary focus:border-transparent placeholder:text-gray-500"
        />
      </SearchField>

      {filteredGroups.map(group => {
        const key = groupKey(group);
        const isExpanded = expanded[key] || !!search.trim();

        return (
          <div key={key} className="mb-4">
            <div className="flex items-center gap-1">
              <Button
                onPress={() => toggleGroup(key)}
                className="flex items-center gap-2 cursor-pointer select-none py-1.5 flex-1 text-left"
              >
                <span
                  className={`text-xs transition-transform text-gray-400 ${isExpanded ? 'rotate-90' : ''}`}
                >
                  ▶
                </span>
                <span className="text-xs uppercase tracking-wide text-gray-400 font-semibold">
                  {group.title}
                </span>
                <span className="text-xs text-gray-600 ml-auto mr-1">
                  {group.items.length}
                </span>
              </Button>
              {group.moduleSlug && (
                <NavLink
                  to={`/docs/${group.moduleSlug}`}
                  end
                  onClick={onItemClick}
                  className={({ isActive }) =>
                    `text-xs px-1.5 py-0.5 rounded transition-colors no-underline ${
                      isActive
                        ? 'bg-primary text-white'
                        : 'text-gray-500 hover:text-primary'
                    }`
                  }
                  title={`${group.title} overview`}
                >
                  overview
                </NavLink>
              )}
            </div>

            {isExpanded && (
              <ul className="pl-4 space-y-0.5 mt-1">
                {group.items.map(item => (
                  <li key={item.linkTo}>
                    <NavLink
                      to={item.linkTo}
                      onClick={onItemClick}
                      className={({ isActive }) =>
                        `block px-3 py-1 rounded text-sm no-underline transition-colors ${
                          isActive
                            ? 'bg-primary text-white'
                            : 'text-gray-300 hover:bg-gray-800'
                        }`
                      }
                    >
                      <span className="inline-block w-5 font-mono text-xs font-bold text-primary">
                        {getKindIcon(item.kind)}
                      </span>
                      {item.name}
                    </NavLink>
                  </li>
                ))}
              </ul>
            )}
          </div>
        );
      })}
    </aside>
  );
}
