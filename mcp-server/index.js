#!/usr/bin/env node

/**
 * MCP server for searching numwasm API documentation.
 * Exposes two tools: search_numwasm_docs and list_numwasm_modules.
 *
 * The docs index is bundled as docs-index.json (generated from api.json).
 */

import { McpServer } from "@modelcontextprotocol/sdk/server/mcp.js";
import { StdioServerTransport } from "@modelcontextprotocol/sdk/server/stdio.js";
import { z } from "zod";
import { readFileSync } from "node:fs";
import { fileURLToPath } from "node:url";
import { dirname, join } from "node:path";

// ---------------------------------------------------------------------------
// Load bundled docs index
// ---------------------------------------------------------------------------

const __dirname = dirname(fileURLToPath(import.meta.url));
const indexPath = join(__dirname, "docs-index.json");
const docsIndex = JSON.parse(readFileSync(indexPath, "utf-8"));
const { overview, sections } = docsIndex;

// ---------------------------------------------------------------------------
// Search logic
// ---------------------------------------------------------------------------

function searchDocs(query, maxResults = 10) {
  const q = query.toLowerCase().trim();
  if (!q) return [];

  const results = [];
  const seen = new Set();

  function add(section) {
    const key = `${section.module || ""}/${section.name}`;
    if (seen.has(key)) return;
    seen.add(key);
    results.push(section);
  }

  // 1. Exact name match
  for (const s of sections) {
    if (s.name.toLowerCase() === q) add(s);
  }

  // 2. Name prefix/substring match
  if (results.length < maxResults) {
    for (const s of sections) {
      if (results.length >= maxResults) break;
      if (s.name.toLowerCase().includes(q)) add(s);
    }
  }

  // 3. Module name match — return all items in that module
  const moduleLimit = 20;
  if (results.length < maxResults) {
    let moduleCount = 0;
    for (const s of sections) {
      if (results.length >= maxResults || moduleCount >= moduleLimit) break;
      if (s.module && s.module.toLowerCase() === q) {
        add(s);
        moduleCount++;
      }
    }
  }

  // 4. Category match
  if (results.length < maxResults) {
    for (const s of sections) {
      if (results.length >= maxResults) break;
      if (s.category && s.category.toLowerCase().includes(q)) add(s);
    }
  }

  // 5. Full-text search on signature + description
  if (results.length < maxResults) {
    for (const s of sections) {
      if (results.length >= maxResults) break;
      const haystack = `${s.signature} ${s.description}`.toLowerCase();
      if (haystack.includes(q)) add(s);
    }
  }

  return results.slice(0, maxResults);
}

// ---------------------------------------------------------------------------
// MCP Server
// ---------------------------------------------------------------------------

const server = new McpServer({
  name: "numwasm-docs",
  version: "0.1.0",
});

server.tool(
  "search_numwasm_docs",
  "Search the numwasm API documentation by function name, module, category, or keyword. Returns matching function signatures, parameters, and descriptions.",
  {
    query: z.string().describe(
      "Search term: function name (e.g. 'matmul'), module name (e.g. 'linalg'), category (e.g. 'statistics'), or keyword (e.g. 'eigenvalue')"
    ),
  },
  async ({ query }) => {
    const results = searchDocs(query);

    if (results.length === 0) {
      return {
        content: [
          {
            type: "text",
            text: `No results found for "${query}". Try a different search term, or use the list_numwasm_modules tool to see available modules and categories.`,
          },
        ],
      };
    }

    const text = results
      .map((s) => s.content)
      .join("\n---\n\n");

    return {
      content: [
        {
          type: "text",
          text: `Found ${results.length} result(s) for "${query}":\n\n${text}`,
        },
      ],
    };
  }
);

server.tool(
  "list_numwasm_modules",
  "List all numwasm modules and API categories with their key functions. Use this to discover what's available before searching for specific functions.",
  {},
  async () => {
    return {
      content: [
        {
          type: "text",
          text: overview,
        },
      ],
    };
  }
);

// ---------------------------------------------------------------------------
// Start
// ---------------------------------------------------------------------------

async function main() {
  const transport = new StdioServerTransport();
  await server.connect(transport);
  console.error("numwasm-mcp server running on stdio");
}

main().catch((error) => {
  console.error("Fatal error:", error);
  process.exit(1);
});
