/**
 * Pre-render all doc routes to static HTML files.
 * Run after: vite build && vite build --ssr
 *
 * Usage: node scripts/prerender.js
 */

import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const __dirname = path.dirname(fileURLToPath(import.meta.url));
const distDir = path.resolve(__dirname, '../dist');
const ssrDir = path.resolve(__dirname, '../dist-ssr');

async function prerender() {
  // Read the built HTML template
  const template = fs.readFileSync(path.join(distDir, 'index.html'), 'utf-8');

  // Import the SSR entry
  const { render, getAllRoutes } = await import(path.join(ssrDir, 'entry-server.js'));

  // Static pages + all doc routes
  const staticRoutes = ['/', '/demo', '/docs', '/benchmarks'];
  const docRoutes = getAllRoutes();
  const allRoutes = [...staticRoutes, ...docRoutes];

  console.log(`Pre-rendering ${allRoutes.length} routes...`);

  let rendered = 0;
  for (const url of allRoutes) {
    try {
      const { html: appHtml, helmet } = render(url);

      // Build head tags from helmet
      const headTags = [
        helmet?.title?.toString() || '',
        helmet?.meta?.toString() || '',
        helmet?.link?.toString() || '',
        helmet?.script?.toString() || '',
      ]
        .filter(Boolean)
        .join('\n    ');

      // Inject app HTML and head tags into the template
      let pageHtml = template
        .replace('<div id="root"></div>', `<div id="root">${appHtml}</div>`);

      // Inject helmet tags before </head>
      if (headTags) {
        pageHtml = pageHtml.replace('</head>', `    ${headTags}\n  </head>`);
      }

      // Write to file
      const filePath = url === '/'
        ? path.join(distDir, 'index.html')
        : path.join(distDir, url, 'index.html');

      const dir = path.dirname(filePath);
      fs.mkdirSync(dir, { recursive: true });
      fs.writeFileSync(filePath, pageHtml);
      rendered++;
    } catch (err) {
      console.warn(`  Warning: Failed to render ${url}: ${err.message}`);
    }
  }

  console.log(`Done. Pre-rendered ${rendered}/${allRoutes.length} routes.`);

  // Generate sitemap.xml
  const baseUrl = 'https://numwasm.dev';
  const sitemap = `<?xml version="1.0" encoding="UTF-8"?>
<urlset xmlns="http://www.sitemaps.org/schemas/sitemap/0.9">
${allRoutes
  .map(
    url => `  <url>
    <loc>${baseUrl}${url}</loc>
    <changefreq>${url === '/' ? 'weekly' : 'monthly'}</changefreq>
    <priority>${url === '/' ? '1.0' : url === '/docs' ? '0.9' : '0.7'}</priority>
  </url>`
  )
  .join('\n')}
</urlset>`;

  fs.writeFileSync(path.join(distDir, 'sitemap.xml'), sitemap);
  console.log(`Generated sitemap.xml with ${allRoutes.length} URLs.`);
}

prerender().catch(err => {
  console.error('Pre-render failed:', err);
  process.exit(1);
});
