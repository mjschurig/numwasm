import { apiData } from '../utils/apiData';
import type { ProjectReflection } from '../types/typedoc';

interface UseApiDocsResult {
  data: ProjectReflection;
  loading: false;
  error: null;
}

/**
 * Synchronous access to the API documentation data.
 * The JSON is statically imported at build time via apiData.ts.
 */
export function useApiDocs(): UseApiDocsResult {
  return { data: apiData, loading: false, error: null };
}
