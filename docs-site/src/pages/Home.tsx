import { Button } from 'react-aria-components';
import { useNavigate } from 'react-router-dom';

export default function Home() {
  const navigate = useNavigate();

  return (
    <main className="text-gray-100">
      <section className="text-center py-20 px-8 max-w-3xl mx-auto relative">
        <div className="hero-glow" />
        <h1 className="text-6xl font-bold hero-text mb-6">numwasm</h1>
        <p className="text-xl text-gray-400 mb-10 leading-relaxed">
          NumPy-inspired n-dimensional array operations in TypeScript with WebAssembly acceleration.
          Full type safety, browser and Node.js support.
        </p>
        <div className="flex gap-4 justify-center">
          <Button
            onPress={() => navigate('/demo')}
            className="cta-button px-8 py-4 rounded-lg font-semibold text-primary hover:text-white pressed:scale-95 transition-all cursor-pointer"
          >
            Try the Demo
          </Button>
          <Button
            onPress={() => navigate('/docs')}
            className="px-8 py-4 bg-gray-800/50 text-gray-200 border border-gray-700 rounded-lg font-semibold hover:bg-gray-700/50 hover:border-gray-600 pressed:scale-95 transition-all cursor-pointer"
          >
            View Documentation
          </Button>
        </div>
      </section>

      <section className="grid grid-cols-1 sm:grid-cols-2 lg:grid-cols-4 gap-6 p-8 max-w-6xl mx-auto">
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">NumPy Compatible</h3>
          <p className="text-gray-400 text-sm">Familiar API following NumPy conventions. Easy migration for Python developers.</p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">WebAssembly Powered</h3>
          <p className="text-gray-400 text-sm">Performance-critical operations compiled to WASM for near-native speed.</p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">Type Safe</h3>
          <p className="text-gray-400 text-sm">Full TypeScript support with strict typing and excellent IDE integration.</p>
        </div>
        <div className="p-6 bg-gray-800/30 border border-gray-700/50 rounded-xl text-center backdrop-blur-sm hover:border-primary/30 transition-colors">
          <h3 className="font-semibold text-primary mb-3">Comprehensive</h3>
          <p className="text-gray-400 text-sm">Linear algebra, FFT, random numbers, broadcasting, and much more.</p>
        </div>
      </section>

      <section className="max-w-3xl mx-auto py-12 px-8">
        <h2 className="text-2xl font-bold text-gray-100 mb-6">Quick Start</h2>
        <pre className="bg-gray-800/50 border border-gray-700/50 p-6 rounded-lg overflow-x-auto backdrop-blur-sm">
          <code className="text-gray-200">{`npm install numwasm

import * as nw from 'numwasm';

// Create arrays
const a = nw.array([[1, 2], [3, 4]]);
const b = nw.array([[5, 6], [7, 8]]);

// Matrix multiplication
const c = nw.linalg.matmul(a, b);
console.log(c.tolist());
// [[19, 22], [43, 50]]`}</code>
        </pre>
      </section>
    </main>
  );
}
