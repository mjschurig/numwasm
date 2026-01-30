import { useState, useCallback } from 'react';
import { Canvas } from '@react-three/fiber';
import { EffectComposer, Bloom } from '@react-three/postprocessing';

import { useGraphData } from './useGraphData';
import { useForceSimulation, type SimNode } from './useForceSimulation';
import { ForceGraph } from './ForceGraph';
import type { GraphFilters, PackageId } from './types';
import { DEFAULT_FILTERS, PACKAGE_COLORS, ALL_PACKAGES } from './types';

export function OverviewGraph() {
  const [filters, setFilters] = useState<GraphFilters>(DEFAULT_FILTERS);
  const [hoveredNode, setHoveredNode] = useState<SimNode | null>(null);
  const { elements, isLoading, error } = useGraphData(filters);
  const { nodes, edges, step, isSettled } = useForceSimulation(elements, filters.searchQuery);

  // Handle hover - use callback to avoid re-renders
  const handleHover = useCallback((node: SimNode | null) => {
    setHoveredNode(node);
  }, []);

  // Toggle package filter
  const togglePackage = (packageId: PackageId) => {
    setFilters((prev) => {
      const newPackages = new Set(prev.packages);
      if (newPackages.has(packageId)) {
        newPackages.delete(packageId);
      } else {
        newPackages.add(packageId);
      }
      return { ...prev, packages: newPackages };
    });
  };

  const toggleWasmFunctions = () => {
    setFilters((prev) => ({
      ...prev,
      showWasmFunctions: !prev.showWasmFunctions,
    }));
  };

  const toggleTsExports = () => {
    setFilters((prev) => ({
      ...prev,
      showTsExports: !prev.showTsExports,
    }));
  };

  const handleSearch = (query: string) => {
    setFilters((prev) => ({ ...prev, searchQuery: query }));
  };

  if (isLoading) {
    return (
      <div className="flex items-center justify-center h-full">
        <div className="text-gray-400">Loading exports data...</div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="flex items-center justify-center h-full">
        <div className="text-red-400">Error: {error}</div>
      </div>
    );
  }

  return (
    <div className="h-full flex flex-col">
      {/* Filter controls */}
      <div className="flex flex-wrap items-center gap-4 p-4 bg-gray-900/50 border-b border-gray-700">
        {/* Package toggles */}
        <div className="flex items-center gap-2 flex-wrap">
          <span className="text-sm text-gray-400">Packages:</span>
          {ALL_PACKAGES.map((pkg) => (
            <button
              key={pkg}
              onClick={() => togglePackage(pkg)}
              className={`px-2 py-1 rounded-full text-xs font-medium transition-colors ${
                filters.packages.has(pkg)
                  ? 'text-white'
                  : 'text-gray-500 bg-gray-800'
              }`}
              style={{
                backgroundColor: filters.packages.has(pkg)
                  ? PACKAGE_COLORS[pkg]
                  : undefined,
              }}
            >
              {pkg.replace('wasm', '')}
            </button>
          ))}
        </div>

        {/* Type toggles */}
        <div className="flex items-center gap-2">
          <span className="text-sm text-gray-400">Show:</span>
          <button
            onClick={toggleWasmFunctions}
            className={`px-3 py-1 rounded-full text-sm font-medium transition-colors ${
              filters.showWasmFunctions
                ? 'bg-emerald-600 text-white'
                : 'bg-gray-800 text-gray-500'
            }`}
          >
            WASM
          </button>
          <button
            onClick={toggleTsExports}
            className={`px-3 py-1 rounded-full text-sm font-medium transition-colors ${
              filters.showTsExports
                ? 'bg-blue-600 text-white'
                : 'bg-gray-800 text-gray-500'
            }`}
          >
            TypeScript
          </button>
        </div>

        {/* Search */}
        <div className="flex items-center gap-2 flex-1 min-w-[200px] max-w-[300px]">
          <input
            type="text"
            placeholder="Search functions..."
            value={filters.searchQuery}
            onChange={(e) => handleSearch(e.target.value)}
            className="flex-1 px-3 py-1.5 bg-gray-800 border border-gray-700 rounded-lg text-sm text-white placeholder-gray-500 focus:outline-none focus:border-primary"
          />
        </div>

        {/* Node count */}
        <div className="text-sm text-gray-500">
          {nodes.length} nodes · {edges.length} edges
          {!isSettled && ' · simulating...'}
        </div>
      </div>

      {/* 3D Graph */}
      <div className="flex-1 min-h-0 relative bg-[#0a0a12]">
        <Canvas
          camera={{ position: [0, 0, 80], fov: 60 }}
          dpr={[1, 2]}
          gl={{
            antialias: true,
            alpha: false,
            powerPreference: 'high-performance',
          }}
        >
          <color attach="background" args={['#0a0a12']} />
          <ForceGraph
            nodes={nodes}
            edges={edges}
            onStep={step}
            isSettled={isSettled}
            hoveredNode={hoveredNode}
            onHover={handleHover}
          />
          <EffectComposer>
            <Bloom
              luminanceThreshold={0.2}
              luminanceSmoothing={0.9}
              intensity={1.5}
              mipmapBlur
            />
          </EffectComposer>
        </Canvas>
      </div>

      {/* Legend */}
      <div className="flex flex-wrap items-center gap-6 p-3 bg-gray-900/50 border-t border-gray-700 text-xs text-gray-400">
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-emerald-500"></span>
          <span>WASM Function</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-primary"></span>
          <span>TypeScript Export</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="w-3 h-3 rounded-full bg-amber-500"></span>
          <span>Shared (multi-package)</span>
        </div>
        <div className="flex items-center gap-2">
          <span className="w-6 border-t-2 border-dashed border-green-500"></span>
          <span>TS→WASM binding</span>
        </div>
        <div className="ml-auto text-gray-500">
          Drag to rotate · Scroll to zoom · Right-click to pan
        </div>
      </div>
    </div>
  );
}
