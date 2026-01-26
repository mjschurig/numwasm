import {
  LineChart,
  Line,
  XAxis,
  YAxis,
  Tooltip,
  Legend,
  ResponsiveContainer,
  CartesianGrid,
} from 'recharts';
import type { BenchmarkOperation } from '../../types/benchmarks';

interface SizeScalingChartProps {
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
    return `${(value * 1000).toFixed(2)}μs`;
  } else if (value < 1) {
    return `${value.toFixed(3)}ms`;
  }
  return `${value.toFixed(2)}ms`;
}

export default function SizeScalingChart({ operation }: SizeScalingChartProps) {
  const chartData = operation.results.map((r) => ({
    size: r.params.size as number,
    sizeLabel: formatSize(r.params.size as number),
    numjs: r.numjs.meanMs,
    numpy: r.numpy.meanMs,
  }));

  return (
    <div className="bg-gray-800/50 border border-gray-700/50 rounded-lg p-4">
      <h3 className="text-lg font-semibold text-gray-200 mb-4">
        Execution Time by Array Size
      </h3>
      <ResponsiveContainer width="100%" height={300}>
        <LineChart data={chartData} margin={{ top: 5, right: 30, left: 20, bottom: 5 }}>
          <CartesianGrid strokeDasharray="3 3" stroke="#374151" />
          <XAxis
            dataKey="sizeLabel"
            stroke="#9ca3af"
            tick={{ fill: '#9ca3af', fontSize: 12 }}
          />
          <YAxis
            scale="log"
            domain={['auto', 'auto']}
            stroke="#9ca3af"
            tick={{ fill: '#9ca3af', fontSize: 12 }}
            tickFormatter={formatMs}
          />
          <Tooltip
            contentStyle={{
              backgroundColor: '#1f2937',
              border: '1px solid #374151',
              borderRadius: '8px',
            }}
            labelStyle={{ color: '#f3f4f6' }}
            formatter={(value, name) => [
              formatMs(value as number),
              name === 'numpy' ? 'NumPy' : 'NumJS',
            ]}
          />
          <Legend
            formatter={(value) => (value === 'numpy' ? 'NumPy' : 'NumJS')}
            wrapperStyle={{ color: '#9ca3af' }}
          />
          <Line
            type="monotone"
            dataKey="numpy"
            stroke="#3b82f6"
            strokeWidth={2}
            dot={{ fill: '#3b82f6', strokeWidth: 2 }}
            activeDot={{ r: 6 }}
          />
          <Line
            type="monotone"
            dataKey="numjs"
            stroke="#22c55e"
            strokeWidth={2}
            dot={{ fill: '#22c55e', strokeWidth: 2 }}
            activeDot={{ r: 6 }}
          />
        </LineChart>
      </ResponsiveContainer>
    </div>
  );
}
