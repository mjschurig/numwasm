/**
 * Generate HTML benchmark report with Chart.js visualizations.
 *
 * Combines NumPy and NumJS benchmark results into an interactive
 * HTML report with charts and data tables.
 */

import { readFileSync, writeFileSync } from 'fs';
import { join, dirname } from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = dirname(__filename);

interface BenchmarkEntry {
  size: number;
  mean_ms: number;
  min_ms: number;
  max_ms: number;
  std_ms: number;
  iterations: number;
  result: number;
}

interface BenchmarkData {
  library: string;
  version: string;
  benchmarks: BenchmarkEntry[];
}

function generateHTML(numpy: BenchmarkData, numjs: BenchmarkData): string {
  const sizes = numpy.benchmarks.map((b) => b.size);
  const numpyTimes = numpy.benchmarks.map((b) => b.mean_ms);
  const numjsTimes = numjs.benchmarks.map((b) => b.mean_ms);
  const speedups = numpy.benchmarks.map(
    (b, i) => b.mean_ms / numjs.benchmarks[i].mean_ms
  );

  // Generate table rows
  const tableRows = numpy.benchmarks
    .map((np, i) => {
      const nj = numjs.benchmarks[i];
      const speedup = np.mean_ms / nj.mean_ms;
      const speedupClass = speedup > 1 ? 'speedup-positive' : 'speedup-negative';
      const speedupText =
        speedup > 1
          ? `${speedup.toFixed(2)}x faster`
          : `${(1 / speedup).toFixed(2)}x slower`;

      return `
      <tr>
        <td>${np.size.toLocaleString()}</td>
        <td>${np.mean_ms.toFixed(4)}</td>
        <td>${nj.mean_ms.toFixed(4)}</td>
        <td class="${speedupClass}">${speedupText}</td>
        <td>${np.iterations}</td>
        <td>${nj.iterations}</td>
      </tr>`;
    })
    .join('');

  return `<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NumPy vs NumJS sum() Benchmark</title>
  <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.1/dist/chart.umd.min.js"></script>
  <style>
    * { box-sizing: border-box; margin: 0; padding: 0; }
    body {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
      padding: 2rem;
      background: #f5f5f5;
      color: #333;
    }
    h1 { margin-bottom: 0.5rem; color: #1a1a1a; }
    .subtitle { color: #666; margin-bottom: 2rem; font-size: 1.1rem; }
    .charts {
      display: grid;
      grid-template-columns: 1fr 1fr;
      gap: 2rem;
      margin-bottom: 2rem;
    }
    .chart-container {
      background: white;
      padding: 1.5rem;
      border-radius: 12px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    }
    .chart-container h2 {
      font-size: 1.1rem;
      margin-bottom: 1rem;
      color: #444;
    }
    canvas { max-height: 400px; }
    table {
      width: 100%;
      border-collapse: collapse;
      background: white;
      border-radius: 12px;
      overflow: hidden;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    }
    th, td {
      padding: 1rem 1.25rem;
      text-align: right;
      border-bottom: 1px solid #eee;
    }
    th {
      background: #f8f9fa;
      font-weight: 600;
      color: #444;
      font-size: 0.9rem;
      text-transform: uppercase;
      letter-spacing: 0.5px;
    }
    td:first-child, th:first-child { text-align: left; }
    tr:last-child td { border-bottom: none; }
    tr:hover td { background: #f8f9fa; }
    .speedup-positive { color: #16a34a; font-weight: 600; }
    .speedup-negative { color: #dc2626; font-weight: 600; }
    .version-info {
      margin-top: 2rem;
      padding: 1rem;
      background: white;
      border-radius: 8px;
      color: #666;
      font-size: 0.875rem;
    }
    .legend {
      display: flex;
      gap: 2rem;
      margin-bottom: 1rem;
      font-size: 0.9rem;
    }
    .legend-item {
      display: flex;
      align-items: center;
      gap: 0.5rem;
    }
    .legend-color {
      width: 16px;
      height: 16px;
      border-radius: 4px;
    }
    .numpy-color { background: #3b82f6; }
    .numjs-color { background: #22c55e; }
    @media (max-width: 900px) {
      .charts { grid-template-columns: 1fr; }
      body { padding: 1rem; }
    }
  </style>
</head>
<body>
  <h1>NumPy vs NumJS sum() Performance Benchmark</h1>
  <p class="subtitle">Comparing array summation performance across different array sizes (10<sup>2</sup> to 10<sup>7</sup> elements)</p>

  <div class="legend">
    <div class="legend-item">
      <div class="legend-color numpy-color"></div>
      <span>NumPy (Python)</span>
    </div>
    <div class="legend-item">
      <div class="legend-color numjs-color"></div>
      <span>NumJS (TypeScript/WASM)</span>
    </div>
  </div>

  <div class="charts">
    <div class="chart-container">
      <h2>Execution Time (logarithmic scale)</h2>
      <canvas id="timeChart"></canvas>
    </div>
    <div class="chart-container">
      <h2>Speedup Ratio (NumPy time / NumJS time)</h2>
      <canvas id="speedupChart"></canvas>
    </div>
  </div>

  <table>
    <thead>
      <tr>
        <th>Array Size</th>
        <th>NumPy (ms)</th>
        <th>NumJS (ms)</th>
        <th>Speedup</th>
        <th>NumPy Iters</th>
        <th>NumJS Iters</th>
      </tr>
    </thead>
    <tbody>
      ${tableRows}
    </tbody>
  </table>

  <div class="version-info">
    <strong>Environment:</strong>
    NumPy v${numpy.version} |
    NumJS v${numjs.version} |
    Generated: ${new Date().toLocaleString()}
  </div>

  <script>
    const sizes = ${JSON.stringify(sizes.map((s) => s.toLocaleString()))};
    const numpyTimes = ${JSON.stringify(numpyTimes)};
    const numjsTimes = ${JSON.stringify(numjsTimes)};
    const speedups = ${JSON.stringify(speedups)};

    // Time comparison chart (log scale)
    new Chart(document.getElementById('timeChart'), {
      type: 'line',
      data: {
        labels: sizes,
        datasets: [
          {
            label: 'NumPy',
            data: numpyTimes,
            borderColor: '#3b82f6',
            backgroundColor: 'rgba(59, 130, 246, 0.1)',
            borderWidth: 2,
            tension: 0.1,
            fill: true,
            pointRadius: 4,
            pointHoverRadius: 6,
          },
          {
            label: 'NumJS',
            data: numjsTimes,
            borderColor: '#22c55e',
            backgroundColor: 'rgba(34, 197, 94, 0.1)',
            borderWidth: 2,
            tension: 0.1,
            fill: true,
            pointRadius: 4,
            pointHoverRadius: 6,
          }
        ]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        interaction: {
          intersect: false,
          mode: 'index',
        },
        scales: {
          y: {
            type: 'logarithmic',
            title: { display: true, text: 'Time (ms)', font: { weight: 'bold' } },
            grid: { color: 'rgba(0,0,0,0.05)' },
          },
          x: {
            title: { display: true, text: 'Array Size', font: { weight: 'bold' } },
            grid: { display: false },
          }
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            callbacks: {
              label: function(context) {
                return context.dataset.label + ': ' + context.parsed.y.toFixed(4) + ' ms';
              }
            }
          }
        }
      }
    });

    // Speedup chart
    new Chart(document.getElementById('speedupChart'), {
      type: 'bar',
      data: {
        labels: sizes,
        datasets: [{
          label: 'Speedup',
          data: speedups,
          backgroundColor: speedups.map(s => s > 1 ? 'rgba(34, 197, 94, 0.8)' : 'rgba(220, 38, 38, 0.8)'),
          borderColor: speedups.map(s => s > 1 ? '#16a34a' : '#dc2626'),
          borderWidth: 1,
          borderRadius: 4,
        }]
      },
      options: {
        responsive: true,
        maintainAspectRatio: true,
        scales: {
          y: {
            title: { display: true, text: 'Speedup Ratio', font: { weight: 'bold' } },
            grid: { color: 'rgba(0,0,0,0.05)' },
            suggestedMin: 0,
          },
          x: {
            title: { display: true, text: 'Array Size', font: { weight: 'bold' } },
            grid: { display: false },
          }
        },
        plugins: {
          legend: { display: false },
          tooltip: {
            callbacks: {
              label: function(context) {
                const val = context.parsed.y;
                if (val > 1) {
                  return 'NumJS is ' + val.toFixed(2) + 'x faster';
                } else {
                  return 'NumPy is ' + (1/val).toFixed(2) + 'x faster';
                }
              }
            }
          }
        }
      }
    });
  </script>
</body>
</html>`;
}

function main() {
  // Use project root benchmark folder (not dist)
  const projectRoot = join(__dirname, '..', '..', '..');
  const resultsDir = join(projectRoot, 'benchmark', 'results');

  console.log('Loading benchmark results...');

  let numpyData: BenchmarkData;
  let numjsData: BenchmarkData;

  try {
    numpyData = JSON.parse(
      readFileSync(join(resultsDir, 'numpy_results.json'), 'utf-8')
    );
    console.log(`  NumPy: ${numpyData.benchmarks.length} benchmarks`);
  } catch {
    console.error('Error: numpy_results.json not found. Run benchmark:numpy first.');
    process.exit(1);
  }

  try {
    numjsData = JSON.parse(
      readFileSync(join(resultsDir, 'numjs_results.json'), 'utf-8')
    );
    console.log(`  NumJS: ${numjsData.benchmarks.length} benchmarks`);
  } catch {
    console.error('Error: numjs_results.json not found. Run benchmark:numjs first.');
    process.exit(1);
  }

  console.log('Generating HTML report...');
  const html = generateHTML(numpyData, numjsData);

  const outputFile = join(projectRoot, 'benchmark', 'index.html');
  writeFileSync(outputFile, html);

  console.log(`Report generated: ${outputFile}`);
}

main();
