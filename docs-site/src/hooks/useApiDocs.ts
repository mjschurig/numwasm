import { apiData, getPackageApiData } from '../utils/apiData';
import type { PackageId } from '../utils/apiData';
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
export function useApiDocs(packageId: PackageId = 'numwasm'): UseApiDocsResult {
  const data = getPackageApiData(packageId) || apiData;
  return { data, loading: false, error: null };
}
