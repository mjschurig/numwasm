import { useState, useEffect } from 'react';
import type { BenchmarkResults } from '../types/benchmarks';

interface UseBenchmarksResult {
  data: BenchmarkResults | null;
  loading: boolean;
  error: string | null;
}

export function useBenchmarks(): UseBenchmarksResult {
  const [data, setData] = useState<BenchmarkResults | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadBenchmarks = async () => {
      try {
        const response = await fetch(import.meta.env.BASE_URL + 'benchmarks.json');
        if (!response.ok) {
          throw new Error(`Failed to load benchmarks: ${response.status}`);
        }
        const json = await response.json();
        setData(json);
      } catch (e) {
        setError(e instanceof Error ? e.message : 'Failed to load benchmarks');
      } finally {
        setLoading(false);
      }
    };

    loadBenchmarks();
  }, []);

  return { data, loading, error };
}
