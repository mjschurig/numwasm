import SEO from '../components/SEO';

export default function MCP() {
  return (
    <main className="text-gray-100 max-w-3xl mx-auto px-4 sm:px-8 py-8 sm:py-12">
      <SEO
        title="MCP Servers"
        description="Set up MCP servers for numwasm, sciwasm, and symwasm to give AI coding assistants searchable access to the full API documentation."
        path="/mcp"
      />

      <h1 className="text-3xl sm:text-4xl font-bold hero-text mb-6">MCP Servers</h1>

      <p className="text-gray-300 text-lg mb-8 leading-relaxed">
        Each package has a dedicated{' '}
        <a
          href="https://modelcontextprotocol.io/"
          target="_blank"
          rel="noopener noreferrer"
          className="text-primary hover:underline"
        >
          Model Context Protocol
        </a>{' '}
        server that gives AI coding assistants searchable access to its API docs.
        They ship with bundled docs indexes — no network calls at runtime.
      </p>

      {/* Available Servers */}
      <section className="mb-10">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4">Available Servers</h2>
        <div className="space-y-4">
          <div className="p-4 bg-gray-800/30 border border-gray-700/50 rounded-lg backdrop-blur-sm">
            <h3 className="font-semibold text-primary mb-1">
              <a href="https://www.npmjs.com/package/numwasm-mcp" target="_blank" rel="noopener noreferrer" className="hover:underline">numwasm-mcp</a>
            </h3>
            <p className="text-gray-400 text-sm mb-2">
              Search numwasm API docs — n-dimensional arrays, linear algebra, FFT, random numbers, and 600+ functions.
            </p>
            <p className="text-gray-500 text-xs">
              Tools: <code>search_numwasm_docs</code>, <code>list_numwasm_modules</code>
            </p>
          </div>
          <div className="p-4 bg-gray-800/30 border border-gray-700/50 rounded-lg backdrop-blur-sm">
            <h3 className="font-semibold text-primary mb-1">
              <a href="https://www.npmjs.com/package/sciwasm-mcp" target="_blank" rel="noopener noreferrer" className="hover:underline">sciwasm-mcp</a>
            </h3>
            <p className="text-gray-400 text-sm mb-2">
              Search sciwasm API docs — optimization, integration, interpolation, statistics, signal processing, and more.
            </p>
            <p className="text-gray-500 text-xs">
              Tools: <code>search_sciwasm_docs</code>, <code>list_sciwasm_modules</code>
            </p>
          </div>
          <div className="p-4 bg-gray-800/30 border border-gray-700/50 rounded-lg backdrop-blur-sm">
            <h3 className="font-semibold text-primary mb-1">
              <a href="https://www.npmjs.com/package/symwasm-mcp" target="_blank" rel="noopener noreferrer" className="hover:underline">symwasm-mcp</a>
            </h3>
            <p className="text-gray-400 text-sm mb-2">
              Search symwasm API docs — symbolic expressions, simplification, equation solving, calculus, and printing.
            </p>
            <p className="text-gray-500 text-xs">
              Tools: <code>search_symwasm_docs</code>, <code>list_symwasm_modules</code>
            </p>
          </div>
        </div>
      </section>

      {/* Setup: Claude Desktop */}
      <section className="mb-10">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4">Claude Desktop</h2>
        <p className="text-gray-300 mb-3">
          Add the following to your{' '}
          <code className="text-sm bg-gray-800/50 px-1.5 py-0.5 rounded">
            claude_desktop_config.json
          </code>
          :
        </p>
        <pre className="bg-gray-800/50 border border-gray-700/50 p-4 sm:p-6 rounded-lg overflow-x-auto backdrop-blur-sm text-sm">
          <code className="text-gray-200">{`{
  "mcpServers": {
    "numwasm-docs": {
      "command": "npx",
      "args": ["-y", "numwasm-mcp"]
    },
    "sciwasm-docs": {
      "command": "npx",
      "args": ["-y", "sciwasm-mcp"]
    },
    "symwasm-docs": {
      "command": "npx",
      "args": ["-y", "symwasm-mcp"]
    }
  }
}`}</code>
        </pre>
      </section>

      {/* Setup: Claude Code */}
      <section className="mb-10">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4">Claude Code</h2>
        <p className="text-gray-300 mb-3">
          Add a{' '}
          <code className="text-sm bg-gray-800/50 px-1.5 py-0.5 rounded">.mcp.json</code>{' '}
          file to your project root:
        </p>
        <pre className="bg-gray-800/50 border border-gray-700/50 p-4 sm:p-6 rounded-lg overflow-x-auto backdrop-blur-sm text-sm">
          <code className="text-gray-200">{`{
  "mcpServers": {
    "numwasm-docs": {
      "command": "npx",
      "args": ["-y", "numwasm-mcp"]
    },
    "sciwasm-docs": {
      "command": "npx",
      "args": ["-y", "sciwasm-mcp"]
    },
    "symwasm-docs": {
      "command": "npx",
      "args": ["-y", "symwasm-mcp"]
    }
  }
}`}</code>
        </pre>
      </section>

      {/* Setup: Other MCP clients */}
      <section className="mb-10">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4">Other MCP Clients</h2>
        <p className="text-gray-300 mb-3">
          Any MCP-compatible client can use the servers. Run them as stdio transports:
        </p>
        <pre className="bg-gray-800/50 border border-gray-700/50 p-4 sm:p-6 rounded-lg overflow-x-auto backdrop-blur-sm text-sm">
          <code className="text-gray-200">{`npx -y numwasm-mcp
npx -y sciwasm-mcp
npx -y symwasm-mcp`}</code>
        </pre>
        <p className="text-gray-400 text-sm mt-3">
          Each server communicates over stdin/stdout using the{' '}
          <a
            href="https://spec.modelcontextprotocol.io/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-primary hover:underline"
          >
            MCP JSON-RPC protocol
          </a>
          .
        </p>
      </section>

      {/* llms.txt alternative */}
      <section className="mb-10">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4">Alternative: llms.txt</h2>
        <p className="text-gray-300 mb-3">
          For clients that support the{' '}
          <a
            href="https://llmstxt.org/"
            target="_blank"
            rel="noopener noreferrer"
            className="text-primary hover:underline"
          >
            llms.txt convention
          </a>
          , the documentation site also serves machine-readable files:
        </p>
        <ul className="list-disc list-inside text-gray-300 space-y-2">
          <li>
            <a
              href="https://numwasm.quebi.de/llms.txt"
              className="text-primary hover:underline"
            >
              llms.txt
            </a>{' '}
            — project overview with module links
          </li>
          <li>
            <a
              href="https://numwasm.quebi.de/llms-full.txt"
              className="text-primary hover:underline"
            >
              llms-full.txt
            </a>{' '}
            — complete API reference (all packages)
          </li>
        </ul>
      </section>
    </main>
  );
}
