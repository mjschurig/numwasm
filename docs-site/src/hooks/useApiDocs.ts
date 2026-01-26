import { useState, useEffect } from 'react';
import type { ProjectReflection } from '../types/typedoc';

interface UseApiDocsResult {
  data: ProjectReflection | null;
  loading: boolean;
  error: string | null;
}

export function useApiDocs(): UseApiDocsResult {
  const [data, setData] = useState<ProjectReflection | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const loadDocs = async () => {
      try {
        const response = await fetch(import.meta.env.BASE_URL + 'api.json');
        if (!response.ok) {
          throw new Error(`Failed to load documentation: ${response.status}`);
        }
        const json = await response.json();
        setData(json);
      } catch (e) {
        setError(e instanceof Error ? e.message : 'Failed to load documentation');
      } finally {
        setLoading(false);
      }
    };

    loadDocs();
  }, []);

  return { data, loading, error };
}
