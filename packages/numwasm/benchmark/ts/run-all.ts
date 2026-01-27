/**
 * Master benchmark runner for NumJS (TypeScript/WASM).
 * Runs all benchmark categories and outputs results to JSON.
 */

import { writeFileSync, mkdirSync } from 'fs';
import { dirname, join } from 'path';
import { fileURLToPath } from 'url';

import { runReduceBenchmarks } from './reduce/aggregations.js';
import { runCoreBenchmarks } from './core/creation.js';
import { runUfuncBenchmarks } from './ufunc/arithmetic.js';
import { runLinalgBenchmarks } from './linalg/operations.js';
import { runRandomBenchmarks } from './random_bench/distributions.js';
import { CATEGORY_NAMES } from './lib/config.js';
import type { CategoryResults } from './lib/types.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

// Parse command line arguments
const args = process.argv.slice(2);
const categoryArg = args.find(a => a.startsWith('--category='));
const selectedCategory = categoryArg?.split('=')[1];

async function main() {
  console.log('='.repeat(60));
  console.log('NumJS Benchmark Suite');
  console.log('='.repeat(60));

  const categories: CategoryResults[] = [];

  // Run selected or all categories
  const runCategory = async (
    id: string,
    runner: () => Promise<any[]>
  ) => {
    if (selectedCategory && selectedCategory !== id) {
      return;
    }
    console.log(`\n[${CATEGORY_NAMES[id] || id}]`);
    const operations = await runner();
    categories.push({
      id,
      name: CATEGORY_NAMES[id] || id,
      operations,
    });
  };

  await runCategory('reduce', runReduceBenchmarks);
  await runCategory('core', runCoreBenchmarks);
  await runCategory('ufunc', runUfuncBenchmarks);
  await runCategory('linalg', runLinalgBenchmarks);
  await runCategory('random', runRandomBenchmarks);

  // Output results
  const projectRoot = join(__dirname, '..', '..', '..');
  const outputDir = join(projectRoot, 'benchmark', 'results', 'numjs');
  mkdirSync(outputDir, { recursive: true });

  // Write individual category files
  for (const category of categories) {
    const outputFile = join(outputDir, `${category.id}.json`);
    writeFileSync(outputFile, JSON.stringify(category, null, 2));
    console.log(`\nResults saved: ${outputFile}`);
  }

  console.log('\n' + '='.repeat(60));
  console.log('NumJS benchmarks complete!');
  console.log('='.repeat(60));
}

main().catch((err) => {
  console.error('Benchmark failed:', err);
  process.exit(1);
});
