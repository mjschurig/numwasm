import type { BenchmarkOperation } from '../../types/benchmarks';

interface BenchmarkTableProps {
  operation: BenchmarkOperation;
}

function formatSize(size: number): string {
  if (size >= 1_000_000) {
    return `${size / 1_000_000}M`;
  } else if (size >= 1_000) {
    return `${size / 1_000}K`;
  }
  return size.toString();
}

function formatMs(value: number): string {
  if (value < 0.001) {
    return `${(value * 1000).toFixed(3)} μs`;
  } else if (value < 1) {
    return `${value.toFixed(4)} ms`;
  }
  return `${value.toFixed(2)} ms`;
}

export default function BenchmarkTable({ operation }: BenchmarkTableProps) {
  return (
    <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg overflow-hidden">
      <div className="overflow-x-auto">
        <table className="w-full">
          <thead>
            <tr className="border-b border-gray-700 bg-gray-900/50">
              <th className="py-3 px-4 text-left text-sm font-semibold text-gray-400 uppercase tracking-wide">
                Size
              </th>
              <th className="py-3 px-4 text-right text-sm font-semibold text-gray-400 uppercase tracking-wide">
                NumPy
              </th>
              <th className="py-3 px-4 text-right text-sm font-semibold text-gray-400 uppercase tracking-wide">
                NumJS
              </th>
              <th className="py-3 px-4 text-right text-sm font-semibold text-gray-400 uppercase tracking-wide">
                Speedup
              </th>
              <th className="py-3 px-4 text-right text-sm font-semibold text-gray-400 uppercase tracking-wide hidden sm:table-cell">
                NumPy Iters
              </th>
              <th className="py-3 px-4 text-right text-sm font-semibold text-gray-400 uppercase tracking-wide hidden sm:table-cell">
                NumJS Iters
              </th>
            </tr>
          </thead>
          <tbody>
            {operation.results.map((result, index) => {
              const speedup = result.speedup;
              const numjsFaster = speedup > 1;
              const speedupText = numjsFaster
                ? `${speedup.toFixed(2)}x faster`
                : `${(1 / speedup).toFixed(2)}x slower`;

              return (
                <tr
                  key={index}
                  className="border-b border-gray-800 hover:bg-gray-800/50 transition-colors"
                >
                  <td className="py-3 px-4 text-gray-200 font-mono">
                    {formatSize(result.params.size as number)}
                  </td>
                  <td className="py-3 px-4 text-right text-gray-300 font-mono">
                    {formatMs(result.numpy.meanMs)}
                  </td>
                  <td className="py-3 px-4 text-right text-gray-300 font-mono">
                    {formatMs(result.numjs.meanMs)}
                  </td>
                  <td
                    className={`py-3 px-4 text-right font-semibold ${
                      numjsFaster ? 'text-green-400' : 'text-red-400'
                    }`}
                  >
                    {speedupText}
                  </td>
                  <td className="py-3 px-4 text-right text-gray-400 font-mono hidden sm:table-cell">
                    {result.numpy.iterations.toLocaleString()}
                  </td>
                  <td className="py-3 px-4 text-right text-gray-400 font-mono hidden sm:table-cell">
                    {result.numjs.iterations.toLocaleString()}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>
    </div>
  );
}
