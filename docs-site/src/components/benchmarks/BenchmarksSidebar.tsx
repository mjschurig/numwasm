import { useState } from 'react';
import { NavLink } from 'react-router-dom';
import { Button, SearchField, Input } from 'react-aria-components';
import { Search } from 'lucide-react';
import type { BenchmarkCategory } from '../../types/benchmarks';

interface BenchmarksSidebarProps {
  categories: BenchmarkCategory[];
  selectedCategory?: string;
  selectedOperation?: string;
  onItemClick?: () => void;
}

export default function BenchmarksSidebar({
  categories,
  selectedCategory: _selectedCategory,
  selectedOperation: _selectedOperation,
  onItemClick,
}: BenchmarksSidebarProps) {
  // selectedCategory and selectedOperation are passed for potential future use
  void _selectedCategory;
  void _selectedOperation;
  const [searchQuery, setSearchQuery] = useState('');
  const [expanded, setExpanded] = useState<Record<string, boolean>>(() => {
    // Expand all categories by default
    const initial: Record<string, boolean> = {};
    categories.forEach((c) => {
      initial[c.id] = true;
    });
    return initial;
  });

  const toggleExpand = (categoryId: string) => {
    setExpanded((prev) => ({ ...prev, [categoryId]: !prev[categoryId] }));
  };

  // Filter operations based on search query
  const filteredCategories = categories
    .map((category) => ({
      ...category,
      operations: category.operations.filter((op) =>
        op.name.toLowerCase().includes(searchQuery.toLowerCase())
      ),
    }))
    .filter((category) => category.operations.length > 0);

  return (
    <aside className="w-72 lg:w-70 border-r border-gray-700 p-4 sm:p-6 overflow-y-auto h-screen lg:sticky lg:top-16 lg:h-[calc(100vh-4rem)] bg-gray-900/95 lg:bg-gray-900/50 backdrop-blur-sm">
      <h2 className="text-lg font-semibold text-gray-200 mb-4">Benchmarks</h2>

      {/* Search */}
      <SearchField
        aria-label="Search operations"
        className="mb-6"
        value={searchQuery}
        onChange={setSearchQuery}
      >
        <div className="relative">
          <Search className="absolute left-3 top-1/2 -translate-y-1/2 w-4 h-4 text-gray-400" />
          <Input
            placeholder="Search operations..."
            className="w-full pl-10 pr-4 py-2 bg-gray-800 border border-gray-700 rounded-lg text-gray-200 placeholder-gray-500 focus:outline-none focus:border-primary focus:ring-1 focus:ring-primary"
          />
        </div>
      </SearchField>

      {/* Overview link */}
      <NavLink
        to="/benchmarks"
        end
        onClick={onItemClick}
        className={({ isActive }) =>
          `block px-3 py-2 mb-4 rounded text-sm no-underline transition-colors ${
            isActive
              ? 'bg-primary text-white'
              : 'text-gray-300 hover:bg-gray-800'
          }`
        }
      >
        Overview
      </NavLink>

      {/* Category navigation */}
      {filteredCategories.map((category) => (
        <div key={category.id} className="mb-4">
          <Button
            onPress={() => toggleExpand(category.id)}
            className="flex items-center gap-2 cursor-pointer select-none py-2 w-full text-left"
          >
            <span
              className={`text-xs transition-transform text-gray-400 ${
                expanded[category.id] ? 'rotate-90' : ''
              }`}
            >
              ▶
            </span>
            <span className="text-xs uppercase tracking-wide text-gray-400 font-semibold">
              {category.name}
            </span>
            <span className="text-xs text-gray-500">
              ({category.operations.length})
            </span>
          </Button>

          {expanded[category.id] && (
            <ul className="pl-4 space-y-1 mt-1">
              {category.operations.map((op) => (
                <li key={op.id}>
                  <NavLink
                    to={`/benchmarks/${category.id}/${op.id}`}
                    onClick={onItemClick}
                    className={({ isActive }) =>
                      `block px-3 py-1.5 rounded text-sm no-underline transition-colors ${
                        isActive
                          ? 'bg-primary text-white'
                          : 'text-gray-300 hover:bg-gray-800'
                      }`
                    }
                  >
                    {op.name}
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
