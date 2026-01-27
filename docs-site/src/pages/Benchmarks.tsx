import { useState } from 'react';
import { useParams, Link } from 'react-router-dom';
import { Button } from 'react-aria-components';
import { Menu, X } from 'lucide-react';
import { useBenchmarks } from '../hooks/useBenchmarks';
import BenchmarksSidebar from '../components/benchmarks/BenchmarksSidebar';
import SizeScalingChart from '../components/benchmarks/SizeScalingChart';
import SpeedupChart from '../components/benchmarks/SpeedupChart';
import BenchmarkTable from '../components/benchmarks/BenchmarkTable';
import SEO from '../components/SEO';
import type { BenchmarkCategory, BenchmarkOperation, BenchmarkResults } from '../types/benchmarks';

function BenchmarksOverview({ data }: { data: BenchmarkResults }) {
  // Calculate summary stats
  const totalOperations = data.categories.reduce(
    (sum, cat) => sum + cat.operations.length,
    0
  );
  const totalDataPoints = data.categories.reduce(
    (sum, cat) =>
      sum + cat.operations.reduce((s, op) => s + op.results.length, 0),
    0
  );

  // Calculate overall speedup stats
  let totalSpeedup = 0;
  let numjsFasterCount = 0;
  let numpyFasterCount = 0;

  data.categories.forEach((cat) => {
    cat.operations.forEach((op) => {
      op.results.forEach((r) => {
        totalSpeedup += r.speedup;
        if (r.speedup > 1) numjsFasterCount++;
        else numpyFasterCount++;
      });
    });
  });

  // avgSpeedup available for future use
  void (totalSpeedup / totalDataPoints);

  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <h1 className="text-3xl font-bold mb-2">Benchmark Results</h1>
        <p className="text-gray-400">
          Performance comparison between NumPy (native Python) and NumJS (TypeScript/WebAssembly)
        </p>
      </div>

      {/* Summary stats */}
      <div className="grid grid-cols-2 md:grid-cols-4 gap-4 mb-8">
        <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
          <div className="text-2xl font-bold text-primary">{data.categories.length}</div>
          <div className="text-sm text-gray-400">Categories</div>
        </div>
        <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
          <div className="text-2xl font-bold text-primary">{totalOperations}</div>
          <div className="text-sm text-gray-400">Operations</div>
        </div>
        <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
          <div className="text-2xl font-bold text-green-400">{numjsFasterCount}</div>
          <div className="text-sm text-gray-400">NumJS Faster</div>
        </div>
        <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
          <div className="text-2xl font-bold text-red-400">{numpyFasterCount}</div>
          <div className="text-sm text-gray-400">NumPy Faster</div>
        </div>
      </div>

      {/* Hardware info */}
      {data.hardware && (
        <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4 mb-8">
          <h2 className="text-lg font-semibold mb-3">Hardware</h2>
          <div className="grid grid-cols-1 md:grid-cols-2 lg:grid-cols-3 gap-4 text-sm">
            <div>
              <span className="text-gray-400">CPU:</span>{' '}
              <span className="text-gray-200">{data.hardware.cpu}</span>
            </div>
            <div>
              <span className="text-gray-400">Cores:</span>{' '}
              <span className="text-gray-200">{data.hardware.cores}</span>
            </div>
            <div>
              <span className="text-gray-400">Memory:</span>{' '}
              <span className="text-gray-200">{data.hardware.memory}</span>
            </div>
            <div>
              <span className="text-gray-400">OS:</span>{' '}
              <span className="text-gray-200">{data.hardware.os}</span>
            </div>
            <div>
              <span className="text-gray-400">Architecture:</span>{' '}
              <span className="text-gray-200">{data.hardware.arch}</span>
            </div>
          </div>
        </div>
      )}

      {/* Environment info */}
      <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4 mb-8">
        <h2 className="text-lg font-semibold mb-3">Software</h2>
        <div className="grid grid-cols-2 md:grid-cols-4 gap-4 text-sm">
          <div>
            <span className="text-gray-400">NumJS:</span>{' '}
            <span className="text-gray-200">v{data.environment.numjsVersion}</span>
          </div>
          <div>
            <span className="text-gray-400">NumPy:</span>{' '}
            <span className="text-gray-200">v{data.environment.numpyVersion}</span>
          </div>
          <div>
            <span className="text-gray-400">Node.js:</span>{' '}
            <span className="text-gray-200">{data.environment.nodeVersion}</span>
          </div>
          <div>
            <span className="text-gray-400">Python:</span>{' '}
            <span className="text-gray-200">{data.environment.pythonVersion}</span>
          </div>
        </div>
        <div className="text-xs text-gray-500 mt-2">
          Generated: {new Date(data.generatedAt).toLocaleString()}
        </div>
      </div>

      {/* Categories grid */}
      <h2 className="text-xl font-semibold mb-4">Categories</h2>
      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
        {data.categories.map((category) => {
          // Calculate category stats
          const catSpeedups = category.operations.flatMap((op) =>
            op.results.map((r) => r.speedup)
          );
          const avgCatSpeedup =
            catSpeedups.reduce((a, b) => a + b, 0) / catSpeedups.length;
          const numjsFaster = catSpeedups.filter((s) => s > 1).length;

          return (
            <Link
              key={category.id}
              to={`/benchmarks/${category.id}/${category.operations[0]?.id}`}
              className="block bg-gray-800/50 border border-gray-700/50 rounded-lg p-4 hover:bg-gray-800 transition-colors no-underline"
            >
              <h3 className="text-lg font-semibold text-primary mb-2">
                {category.name}
              </h3>
              <div className="text-sm text-gray-400 mb-2">
                {category.operations.length} operations
              </div>
              <div className="text-sm">
                <span className="text-gray-400">Avg speedup: </span>
                <span
                  className={avgCatSpeedup > 1 ? 'text-green-400' : 'text-red-400'}
                >
                  {avgCatSpeedup.toFixed(2)}x
                </span>
              </div>
              <div className="text-sm">
                <span className="text-green-400">{numjsFaster}</span>
                <span className="text-gray-500"> / </span>
                <span className="text-red-400">
                  {catSpeedups.length - numjsFaster}
                </span>
                <span className="text-gray-400"> (NumJS/NumPy faster)</span>
              </div>
            </Link>
          );
        })}
      </div>
    </div>
  );
}

function CategoryView({ category }: { category: BenchmarkCategory }) {
  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <span className="text-gray-400 text-sm uppercase tracking-wide">Category</span>
        <h1 className="text-3xl font-bold">{category.name}</h1>
      </div>

      <div className="grid gap-4 md:grid-cols-2 lg:grid-cols-3">
        {category.operations.map((op) => {
          const avgSpeedup =
            op.results.reduce((sum, r) => sum + r.speedup, 0) / op.results.length;

          return (
            <Link
              key={op.id}
              to={`/benchmarks/${category.id}/${op.id}`}
              className="block bg-gray-800/50 border border-gray-700/50 rounded-lg p-4 hover:bg-gray-800 transition-colors no-underline"
            >
              <h3 className="text-lg font-semibold text-primary mb-2">{op.name}</h3>
              <div className="text-sm text-gray-400 mb-2">
                {op.results.length} data points
              </div>
              <div className="text-sm">
                <span className="text-gray-400">Avg speedup: </span>
                <span
                  className={avgSpeedup > 1 ? 'text-green-400' : 'text-red-400'}
                >
                  {avgSpeedup.toFixed(2)}x
                </span>
              </div>
            </Link>
          );
        })}
      </div>
    </div>
  );
}

function OperationView({
  operation,
  categoryName,
}: {
  operation: BenchmarkOperation;
  categoryName: string;
}) {
  return (
    <div className="text-gray-100">
      <div className="border-b border-gray-700 pb-4 mb-8">
        <span className="text-gray-400 text-sm uppercase tracking-wide">
          {categoryName}
        </span>
        <h1 className="text-3xl font-bold">{operation.name}</h1>
      </div>

      {/* Charts */}
      <div className="grid gap-6 lg:grid-cols-2 mb-8">
        <SizeScalingChart operation={operation} />
        <SpeedupChart operation={operation} />
      </div>

      {/* Data table */}
      <h2 className="text-xl font-semibold mb-4">Detailed Results</h2>
      <BenchmarkTable operation={operation} />
    </div>
  );
}

export default function Benchmarks() {
  const [sidebarOpen, setSidebarOpen] = useState(false);
  const { '*': path } = useParams();
  const { data, loading, error } = useBenchmarks();

  if (loading) {
    return (
      <div className="flex min-h-[calc(100vh-4rem)]">
        <div className="flex-1 p-4 sm:p-8">
          <div className="flex justify-center items-center min-h-[200px] text-gray-400">
            Loading benchmarks...
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
            <h2 className="text-xl font-bold mb-2">Error Loading Benchmarks</h2>
            <p>{error || 'Failed to load benchmark data'}</p>
            <p className="mt-4 text-gray-400">
              Make sure to run <code className="bg-gray-800 px-2 py-1 rounded">npm run benchmark</code> to generate the benchmark data.
            </p>
          </div>
        </div>
      </div>
    );
  }

  // Parse path: /benchmarks/category/operation
  const pathParts = (path || '').split('/').filter(Boolean);
  const categoryId = pathParts[0];
  const operationId = pathParts[1];

  const selectedCategory = categoryId
    ? data.categories.find((c) => c.id === categoryId)
    : null;

  const selectedOperation =
    selectedCategory && operationId
      ? selectedCategory.operations.find((o) => o.id === operationId)
      : null;

  return (
    <div className="flex min-h-[calc(100vh-4rem)]">
      <SEO
        title="Benchmarks"
        description="Performance benchmarks comparing numwasm WebAssembly operations against native JavaScript. See speedup charts and detailed results."
        path="/benchmarks"
        breadcrumbs={[
          { name: 'Home', url: '/' },
          { name: 'Benchmarks', url: '/benchmarks' },
        ]}
      />
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
        <BenchmarksSidebar
          categories={data.categories}
          selectedCategory={categoryId}
          selectedOperation={operationId}
          onItemClick={() => setSidebarOpen(false)}
        />
      </div>

      <main className="flex-1 p-4 sm:p-8 max-w-full overflow-x-hidden">
        {selectedOperation ? (
          <OperationView
            operation={selectedOperation}
            categoryName={selectedCategory?.name || ''}
          />
        ) : selectedCategory ? (
          <CategoryView category={selectedCategory} />
        ) : (
          <BenchmarksOverview data={data} />
        )}
      </main>
    </div>
  );
}
