import { useState, useMemo } from 'react';
import { NavLink } from 'react-router-dom';
import { SearchField, Input, Button } from 'react-aria-components';
import type { ProjectReflection, DeclarationReflection } from '../../types/typedoc';
import { ReflectionKind } from '../../types/typedoc';

interface DocsSidebarProps {
  data: ProjectReflection;
  onItemClick?: () => void;
}

interface SidebarGroup {
  title: string;
  items: DeclarationReflection[];
}

export default function DocsSidebar({ data, onItemClick }: DocsSidebarProps) {
  const [search, setSearch] = useState('');
  const [expanded, setExpanded] = useState<Record<string, boolean>>({
    Core: true,
    'Linear Algebra': true,
  });

  const groups = useMemo((): SidebarGroup[] => {
    const children = data.children || [];

    // Categorize items
    const core: DeclarationReflection[] = [];
    const linalg: DeclarationReflection[] = [];
    const random: DeclarationReflection[] = [];
    const fft: DeclarationReflection[] = [];
    const types: DeclarationReflection[] = [];
    const other: DeclarationReflection[] = [];

    for (const child of children) {
      const name = child.name.toLowerCase();

      // Group by known categories
      if (name === 'ndarray' || name === 'array' || ['zeros', 'ones', 'arange', 'linspace', 'empty', 'full', 'eye', 'identity'].includes(name)) {
        core.push(child);
      } else if (name.includes('linalg') || ['matmul', 'dot', 'inv', 'det', 'solve', 'eig', 'svd', 'qr', 'cholesky', 'norm'].includes(name)) {
        linalg.push(child);
      } else if (name.includes('random')) {
        random.push(child);
      } else if (name.includes('fft')) {
        fft.push(child);
      } else if (child.kind === ReflectionKind.Interface || child.kind === ReflectionKind.TypeAlias || child.kind === ReflectionKind.Enum) {
        types.push(child);
      } else {
        other.push(child);
      }
    }

    const result: SidebarGroup[] = [];
    if (core.length > 0) result.push({ title: 'Core', items: core });
    if (linalg.length > 0) result.push({ title: 'Linear Algebra', items: linalg });
    if (random.length > 0) result.push({ title: 'Random', items: random });
    if (fft.length > 0) result.push({ title: 'FFT', items: fft });
    if (types.length > 0) result.push({ title: 'Types', items: types });
    if (other.length > 0) result.push({ title: 'Other', items: other });

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

  const toggleGroup = (title: string) => {
    setExpanded(prev => ({ ...prev, [title]: !prev[title] }));
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
      default: return '';
    }
  };

  return (
    <aside className="w-72 lg:w-70 border-r border-gray-700 p-4 sm:p-6 overflow-y-auto h-screen lg:sticky lg:top-16 lg:h-[calc(100vh-4rem)] bg-gray-900/95 lg:bg-gray-900/50 backdrop-blur-sm">
      <SearchField value={search} onChange={setSearch} className="mb-4">
        <Input
          placeholder="Search..."
          className="w-full px-3 py-2 border border-gray-700 bg-gray-800/50 text-gray-200 rounded-lg text-sm focus:outline-none focus:ring-2 focus:ring-primary focus:border-transparent placeholder:text-gray-500"
        />
      </SearchField>

      {filteredGroups.map(group => (
        <div key={group.title} className="mb-6">
          <Button
            onPress={() => toggleGroup(group.title)}
            className="flex items-center gap-2 cursor-pointer select-none py-2 w-full text-left"
          >
            <span className={`text-xs transition-transform text-gray-400 ${expanded[group.title] ? 'rotate-90' : ''}`}>
              ▶
            </span>
            <span className="text-xs uppercase tracking-wide text-gray-400 font-semibold">
              {group.title}
            </span>
          </Button>

          {expanded[group.title] && (
            <ul className="pl-4 space-y-1">
              {group.items.map(item => (
                <li key={item.id}>
                  <NavLink
                    to={`/docs/${item.name}`}
                    onClick={onItemClick}
                    className={({ isActive }) =>
                      `block px-3 py-1.5 rounded text-sm no-underline transition-colors ${
                        isActive
                          ? 'bg-primary text-white'
                          : 'text-gray-300 hover:bg-gray-800'
                      }`
                    }
                  >
                    <span className="inline-block w-6 font-mono text-xs font-bold text-primary">
                      {getKindIcon(item.kind)}
                    </span>
                    {item.name}
                  </NavLink>
                </li>
              ))}
            </ul>
          )}
        </div>
      ))}
    </aside>
  );
}
