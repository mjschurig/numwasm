import { useState, useCallback } from 'react';
import { Button } from 'react-aria-components';
import MatrixDisplay from './MatrixDisplay';

interface MatrixData {
  A: number[][];
  B: number[][];
  C: number[][];
}

function generateRandomMatrix(): number[][] {
  return Array.from({ length: 3 }, () =>
    Array.from({ length: 3 }, () => Math.random() * 10 - 5)
  );
}

function multiplyMatrices(a: number[][], b: number[][]): number[][] {
  const result: number[][] = [];
  for (let i = 0; i < 3; i++) {
    result[i] = [];
    for (let j = 0; j < 3; j++) {
      let sum = 0;
      for (let k = 0; k < 3; k++) {
        sum += a[i][k] * b[k][j];
      }
      result[i][j] = sum;
    }
  }
  return result;
}

export default function MatrixDemo() {
  const [matrices, setMatrices] = useState<MatrixData | null>(null);

  const generateMatrices = useCallback(() => {
    const A = generateRandomMatrix();
    const B = generateRandomMatrix();
    const C = multiplyMatrices(A, B);
    setMatrices({ A, B, C });
  }, []);

  return (
    <div className="max-w-5xl mx-auto p-4 sm:p-8 text-gray-100">
      <h1 className="text-2xl sm:text-3xl font-bold mb-3 sm:mb-4">Matrix Multiplication Demo</h1>
      <p className="text-gray-400 mb-6 sm:mb-8 text-sm sm:text-base">
        This demo generates two random 3x3 matrices and computes their product using
        the standard matrix multiplication algorithm. Click the button to generate new matrices.
      </p>

      <div className="text-center mb-6 sm:mb-8">
        <Button
          onPress={generateMatrices}
          className="cta-button px-5 sm:px-6 py-2.5 sm:py-3 rounded-lg font-semibold text-primary hover:text-white pressed:scale-95 transition-all cursor-pointer"
        >
          Generate Random Matrices
        </Button>
      </div>

      {matrices && (
        <div className="flex flex-col sm:flex-row items-center justify-center gap-3 sm:gap-4 flex-wrap my-6 sm:my-8">
          <MatrixDisplay data={matrices.A} label="Matrix A" />
          <span className="text-2xl sm:text-3xl font-bold text-gray-500">&times;</span>
          <MatrixDisplay data={matrices.B} label="Matrix B" />
          <span className="text-2xl sm:text-3xl font-bold text-gray-500">=</span>
          <MatrixDisplay data={matrices.C} label="Result (A &times; B)" />
        </div>
      )}

      {matrices && (
        <div className="mt-6 sm:mt-8 p-4 sm:p-6 bg-gray-800/50 border border-gray-700/50 rounded-lg backdrop-blur-sm">
          <h3 className="font-semibold mb-2 text-gray-100">How it works</h3>
          <p className="text-gray-400 mb-4 text-sm sm:text-base">
            Matrix multiplication computes each element C[i,j] as the dot product of row i from A
            and column j from B:
          </p>
          <pre className="text-xs sm:text-sm bg-gray-900/50 border border-gray-700/50 p-3 sm:p-4 rounded-lg overflow-x-auto">
            <code className="text-gray-200">{`C[i,j] = sum(A[i,k] * B[k,j]) for k = 0 to 2

// With numwasm:
import * as nw from 'numwasm';
const A = nw.random.random([3, 3]);
const B = nw.random.random([3, 3]);
const C = nw.linalg.matmul(A, B);`}</code>
          </pre>
        </div>
      )}
    </div>
  );
}
