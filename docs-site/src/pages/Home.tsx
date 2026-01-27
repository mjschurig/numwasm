import { Button } from "react-aria-components";
import { useNavigate } from "react-router-dom";
import SEO from "../components/SEO";

export default function Home() {
  const navigate = useNavigate();

  return (
    <main className="text-gray-100">
      <SEO
        title="*wasm - Scientific Computing for TypeScript"
        description="NumPy, SciPy, and SymPy for TypeScript. N-dimensional arrays, scientific computing, and symbolic math with WebAssembly acceleration."
        path="/"
      />
      <section className="text-center py-12 sm:py-20 px-4 sm:px-8 max-w-3xl mx-auto relative">
        <div className="hero-glow" />
        <h1 className="text-4xl sm:text-6xl font-bold hero-text mb-4 sm:mb-6">
          *wasm
        </h1>
        <p className="text-lg sm:text-xl text-gray-400 mb-8 sm:mb-10 leading-relaxed">
          Scientific computing for TypeScript. NumPy-style arrays, SciPy-style
          algorithms, and SymPy-style symbolic math — all with WebAssembly
          acceleration.
        </p>
        <div className="flex flex-col sm:flex-row gap-3 sm:gap-4 justify-center">
          <Button
            onPress={() => navigate("/demo")}
            className="cta-button px-6 sm:px-8 py-3 sm:py-4 rounded-lg font-semibold text-primary hover:text-white pressed:scale-95 transition-all cursor-pointer"
          >
            Try the Demo
          </Button>
          <Button
            onPress={() => navigate("/docs")}
            className="px-6 sm:px-8 py-3 sm:py-4 bg-gray-800/50 text-gray-200 border border-gray-700 rounded-lg font-semibold hover:bg-gray-700/50 hover:border-gray-600 pressed:scale-95 transition-all cursor-pointer"
          >
            View Documentation
          </Button>
        </div>
      </section>

      <section className="grid grid-cols-1 sm:grid-cols-3 gap-4 sm:gap-6 p-4 sm:p-8 max-w-6xl mx-auto">
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">numwasm</h3>
          <p className="text-gray-400 text-sm">
            NumPy-compatible n-dimensional arrays. Linear algebra, FFT, random
            numbers, broadcasting, and 600+ functions.
          </p>
          <Button
            onPress={() => navigate("/docs/numwasm")}
            className="mt-4 px-4 py-2 text-xs bg-gray-700/50 text-gray-200 border border-gray-600 rounded-md font-medium hover:bg-gray-600/50 cursor-pointer transition-colors"
          >
            Docs →
          </Button>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">sciwasm</h3>
          <p className="text-gray-400 text-sm">
            SciPy-compatible scientific computing. Optimization, integration,
            interpolation, signal processing, and statistics.
          </p>
          <Button
            onPress={() => navigate("/docs/sciwasm")}
            className="mt-4 px-4 py-2 text-xs bg-gray-700/50 text-gray-200 border border-gray-600 rounded-md font-medium hover:bg-gray-600/50 cursor-pointer transition-colors"
          >
            Docs →
          </Button>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">symwasm</h3>
          <p className="text-gray-400 text-sm">
            SymPy-compatible symbolic math. Symbolic expressions,
            simplification, equation solving, calculus, and printing.
          </p>
          <Button
            onPress={() => navigate("/docs/symwasm")}
            className="mt-4 px-4 py-2 text-xs bg-gray-700/50 text-gray-200 border border-gray-600 rounded-md font-medium hover:bg-gray-600/50 cursor-pointer transition-colors"
          >
            Docs →
          </Button>
        </div>
      </section>

      <section className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-4 sm:gap-6 px-4 sm:px-8 max-w-6xl mx-auto">
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">Python Compatible</h3>
          <p className="text-gray-400 text-sm">
            Familiar APIs following NumPy, SciPy, and SymPy conventions.
          </p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">
            WebAssembly Powered
          </h3>
          <p className="text-gray-400 text-sm">
            Performance-critical operations compiled to WASM for near-native
            speed.
          </p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">Type Safe</h3>
          <p className="text-gray-400 text-sm">
            Full TypeScript support with strict typing and excellent IDE
            integration.
          </p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">Universal</h3>
          <p className="text-gray-400 text-sm">
            Works in browsers and Node.js. ESM and CommonJS support.
          </p>
        </div>
      </section>

      <section className="max-w-3xl mx-auto py-8 sm:py-12 px-4 sm:px-8">
        <h2 className="text-xl sm:text-2xl font-bold text-gray-100 mb-4 sm:mb-6">
          Quick Start
        </h2>
        <pre className="bg-gray-800/50 border border-gray-700/50 p-4 sm:p-6 rounded-lg overflow-x-auto backdrop-blur-sm text-sm sm:text-base">
          <code className="text-gray-200">{`npm install numwasm sciwasm symwasm

import * as nw from 'numwasm';
import * as sci from 'sciwasm';
import * as sym from 'symwasm';

// numwasm: N-dimensional arrays
const a = nw.array([[1, 2], [3, 4]]);
const b = nw.linalg.matmul(a, a);

// sciwasm: Scientific computing
const result = sci.optimize.minimize(f, x0);

// symwasm: Symbolic math
const x = new sym.core.Symbol('x');
const deriv = sym.calculus.diff(expr, x);`}</code>
        </pre>
      </section>
    </main>
  );
}
