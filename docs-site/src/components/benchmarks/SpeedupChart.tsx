import {
  BarChart,
  Bar,
  XAxis,
  YAxis,
  Tooltip,
  ResponsiveContainer,
  Cell,
  ReferenceLine,
  CartesianGrid,
} from 'recharts';
import type { BenchmarkOperation } from '../../types/benchmarks';

interface SpeedupChartProps {
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

export default function SpeedupChart({ operation }: SpeedupChartProps) {
  const chartData = operation.results.map((r) => ({
    label: formatSize(r.params.size as number),
    speedup: r.speedup,
    numjsFaster: r.speedup > 1,
  }));

  return (
    <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
      <h3 className="text-lg font-semibold text-gray-200 mb-4">
        Speedup Ratio (NumJS vs NumPy)
      </h3>
      <p className="text-sm text-gray-400 mb-4">
        Values &gt; 1 mean NumJS is faster. Values &lt; 1 mean NumPy is faster.
      </p>
      <ResponsiveContainer width="100%" height={250}>
        <BarChart data={chartData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#374151" vertical={false} />
          <XAxis
            dataKey="label"
            stroke="#9ca3af"
            tick={{ fill: '#9ca3af', fontSize: 12 }}
          />
          <YAxis
            stroke="#9ca3af"
            tick={{ fill: '#9ca3af', fontSize: 12 }}
            tickFormatter={(v) => `${v}x`}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: '#1f2937',
              border: '1px solid #374151',
              borderRadius: '8px',
            }}
            labelStyle={{ color: '#f3f4f6' }}
            formatter={(value) => {
              const v = value as number;
              if (v > 1) {
                return [`NumJS is ${v.toFixed(2)}x faster`, 'Speedup'];
              } else {
                return [`NumPy is ${(1 / v).toFixed(2)}x faster`, 'Speedup'];
              }
            }}
          />
          <ReferenceLine y={1} stroke="#6b7280" strokeDasharray="3 3" />
          <Bar dataKey="speedup" radius={[4, 4, 0, 0]}>
            {chartData.map((entry, index) => (
              <Cell
                key={`cell-${index}`}
                fill={entry.numjsFaster ? '#22c55e' : '#ef4444'}
              />
            ))}
          </Bar>
        </BarChart>
      </ResponsiveContainer>
    </div>
  );
}
