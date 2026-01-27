/**
 * Combines NumPy and NumJS benchmark results into a single JSON file
 * for the docs-site.
 */

import { readFileSync, writeFileSync, readdirSync, existsSync, copyFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';
import { execSync } from 'child_process';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

interface BenchmarkDataPoint {
  params: Record<string, string | number>;
  meanMs: number;
  minMs: number;
  maxMs: number;
  stdMs: number;
  iterations: number;
  result?: number;
}

interface OperationResults {
  id: string;
  name: string;
  results: BenchmarkDataPoint[];
}

interface CategoryResults {
  id: string;
  name: string;
  operations: OperationResults[];
}

interface CombinedDataPoint {
  params: Record<string, string | number>;
  numpy: {
    meanMs: number;
    minMs: number;
    maxMs: number;
    stdMs: number;
    iterations: number;
  };
  numjs: {
    meanMs: number;
    minMs: number;
    maxMs: number;
    stdMs: number;
    iterations: number;
  };
  speedup: number;
}

interface CombinedOperation {
  id: string;
  name: string;
  results: CombinedDataPoint[];
}

interface CombinedCategory {
  id: string;
  name: string;
  operations: CombinedOperation[];
}

interface HardwareInfo {
  cpu: string;
  cores: number;
  memory: string;
  os: string;
  arch: string;
}

interface BenchmarkResults {
  version: string;
  generatedAt: string;
  environment: {
    numjsVersion: string;
    numpyVersion: string;
    nodeVersion: string;
    pythonVersion: string;
    platform: string;
  };
  hardware: HardwareInfo;
  categories: CombinedCategory[];
}

import * as os from 'os';

function getVersion(cmd: string): string {
  try {
    return execSync(cmd, { encoding: 'utf-8', stdio: ['pipe', 'pipe', 'pipe'] }).trim();
  } catch {
    return 'unknown';
  }
}

function getHardwareInfo(): HardwareInfo {
  // Get CPU info
  let cpuModel = 'Unknown';
  const cpus = os.cpus();
  if (cpus.length > 0) {
    cpuModel = cpus[0].model.trim();
  }

  // Get total memory in GB
  const totalMemGB = (os.totalmem() / (1024 * 1024 * 1024)).toFixed(1);

  // Get OS info
  let osInfo = `${os.type()} ${os.release()}`;

  // Try to get more detailed OS info on Linux
  if (process.platform === 'linux') {
    try {
      const osRelease = readFileSync('/etc/os-release', 'utf-8');
      const prettyName = osRelease.match(/PRETTY_NAME="([^"]+)"/);
      if (prettyName) {
        osInfo = prettyName[1];
      }
    } catch {
      // Fall back to basic info
    }
  } else if (process.platform === 'darwin') {
    try {
      const swVers = execSync('sw_vers -productVersion', { encoding: 'utf-8' }).trim();
      osInfo = `macOS ${swVers}`;
    } catch {
      osInfo = 'macOS';
    }
  } else if (process.platform === 'win32') {
    osInfo = `Windows ${os.release()}`;
  }

  return {
    cpu: cpuModel,
    cores: cpus.length,
    memory: `${totalMemGB} GB`,
    os: osInfo,
    arch: os.arch(),
  };
}

function combineResults() {
  const projectRoot = join(__dirname, '..', '..', '..');
  const resultsDir = join(projectRoot, 'benchmark', 'results');
  const numpyDir = join(resultsDir, 'numpy');
  const numjsDir = join(resultsDir, 'numjs');

  if (!existsSync(numpyDir)) {
    console.error('NumPy results directory not found. Run benchmark:numpy first.');
    process.exit(1);
  }

  if (!existsSync(numjsDir)) {
    console.error('NumJS results directory not found. Run benchmark:numjs first.');
    process.exit(1);
  }

  // Get list of categories from numpy results
  const categories = readdirSync(numpyDir)
    .filter(f => f.endsWith('.json'))
    .map(f => f.replace('.json', ''));

  console.log(`Found ${categories.length} categories: ${categories.join(', ')}`);

  // Get environment info
  const numpyVersion = getVersion('python3 -c "import numpy; print(numpy.__version__)"');
  const pythonVersion = getVersion('python3 --version').replace('Python ', '');

  const hardware = getHardwareInfo();
  console.log(`Hardware: ${hardware.cpu} (${hardware.cores} cores), ${hardware.memory} RAM`);

  const combined: BenchmarkResults = {
    version: '1.0',
    generatedAt: new Date().toISOString(),
    environment: {
      numjsVersion: '0.1.0',
      numpyVersion,
      nodeVersion: process.version,
      pythonVersion,
      platform: process.platform,
    },
    hardware,
    categories: [],
  };

  // Process each category
  for (const categoryId of categories) {
    const numpyFile = join(numpyDir, `${categoryId}.json`);
    const numjsFile = join(numjsDir, `${categoryId}.json`);

    if (!existsSync(numjsFile)) {
      console.warn(`Warning: NumJS results for '${categoryId}' not found, skipping.`);
      continue;
    }

    const numpyData: CategoryResults = JSON.parse(readFileSync(numpyFile, 'utf-8'));
    const numjsData: CategoryResults = JSON.parse(readFileSync(numjsFile, 'utf-8'));

    const combinedCategory: CombinedCategory = {
      id: categoryId,
      name: numpyData.name,
      operations: [],
    };

    // Process each operation
    for (const npOp of numpyData.operations) {
      const njOp = numjsData.operations.find(o => o.id === npOp.id);

      if (!njOp) {
        console.warn(`Warning: NumJS operation '${npOp.id}' not found in category '${categoryId}'.`);
        continue;
      }

      const combinedOp: CombinedOperation = {
        id: npOp.id,
        name: npOp.name,
        results: [],
      };

      // Match results by params
      for (const npResult of npOp.results) {
        const njResult = njOp.results.find(r =>
          JSON.stringify(r.params) === JSON.stringify(npResult.params)
        );

        if (!njResult) {
          console.warn(`Warning: NumJS result for params ${JSON.stringify(npResult.params)} not found.`);
          continue;
        }

        const speedup = njResult.meanMs > 0 ? npResult.meanMs / njResult.meanMs : 0;

        combinedOp.results.push({
          params: npResult.params,
          numpy: {
            meanMs: npResult.meanMs,
            minMs: npResult.minMs,
            maxMs: npResult.maxMs,
            stdMs: npResult.stdMs,
            iterations: npResult.iterations,
          },
          numjs: {
            meanMs: njResult.meanMs,
            minMs: njResult.minMs,
            maxMs: njResult.maxMs,
            stdMs: njResult.stdMs,
            iterations: njResult.iterations,
          },
          speedup,
        });
      }

      combinedCategory.operations.push(combinedOp);
    }

    combined.categories.push(combinedCategory);
  }

  // Write combined results
  const outputFile = join(resultsDir, 'benchmarks.json');
  writeFileSync(outputFile, JSON.stringify(combined, null, 2));
  console.log(`\nCombined results written to: ${outputFile}`);

  // Copy to docs-site/public
  const docsPublicDir = join(projectRoot, 'docs-site', 'public');
  const docsBenchmarksFile = join(docsPublicDir, 'benchmarks.json');
  copyFileSync(outputFile, docsBenchmarksFile);
  console.log(`Copied to: ${docsBenchmarksFile}`);

  // Print summary
  console.log('\nSummary:');
  console.log(`  Categories: ${combined.categories.length}`);
  const totalOps = combined.categories.reduce((sum, c) => sum + c.operations.length, 0);
  console.log(`  Operations: ${totalOps}`);
  const totalDataPoints = combined.categories.reduce(
    (sum, c) => sum + c.operations.reduce((s, o) => s + o.results.length, 0),
    0
  );
  console.log(`  Data points: ${totalDataPoints}`);
}

combineResults();
