/**
 * NumJS Einstein Summation Module
 *
 * Provides NumPy-compatible einsum and einsum_path operations.
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";
import { tensordot } from "./linalg.js";
import { diagonal as diagonalFn } from "./indexing.js";

/**
 * Error for einsum-specific validation errors.
 */
class EinsumError extends Error {
  constructor(message: string) {
    super(message);
    this.name = "EinsumError";
  }
}

/**
 * Parsed einsum subscripts.
 */
interface ParsedEinsum {
  inputSubs: string[];
  outputSubs: string | null;
}

/**
 * Parse einsum subscripts string.
 */
function parseEinsum(subscripts: string): ParsedEinsum {
  // Remove whitespace
  const clean = subscripts.replace(/\s/g, "");

  // Split on '->'
  const parts = clean.split("->");

  if (parts.length > 2) {
    throw new EinsumError(`Invalid subscripts: more than one '->' found`);
  }

  // Parse input subscripts
  const inputPart = parts[0];
  const inputSubs = inputPart.split(",");

  // Parse output subscripts
  const outputSubs = parts.length === 2 ? parts[1] : null;

  // Validate characters (only a-zA-Z allowed)
  const validChars = /^[a-zA-Z]*$/;
  for (const sub of inputSubs) {
    if (!validChars.test(sub)) {
      throw new EinsumError(
        `Invalid subscript character in '${sub}'. Only letters a-zA-Z are allowed.`,
      );
    }
  }

  if (outputSubs !== null && !validChars.test(outputSubs)) {
    throw new EinsumError(
      `Invalid output subscript character in '${outputSubs}'. Only letters a-zA-Z are allowed.`,
    );
  }

  return { inputSubs, outputSubs };
}

/**
 * Compute implicit output subscripts (alphabetical order of non-repeated labels).
 */
function implicitOutput(inputSubs: string[]): string {
  const counts = new Map<string, number>();

  for (const subs of inputSubs) {
    for (const label of subs) {
      counts.set(label, (counts.get(label) ?? 0) + 1);
    }
  }

  // Labels that appear exactly once, sorted alphabetically
  return Array.from(counts.entries())
    .filter(([_, count]) => count === 1)
    .map(([label, _]) => label)
    .sort()
    .join("");
}

/**
 * Build dimension map from operands and subscripts.
 */
function buildDimMap(
  inputSubs: string[],
  operands: NDArray[],
): Map<string, number> {
  const dimMap = new Map<string, number>();

  for (let i = 0; i < inputSubs.length; i++) {
    const subs = inputSubs[i];
    const arr = operands[i];

    if (subs.length !== arr.ndim) {
      throw new EinsumError(
        `Operand ${i} has ${arr.ndim} dimensions but subscripts specify ${subs.length}`,
      );
    }

    for (let j = 0; j < subs.length; j++) {
      const label = subs[j];
      const size = arr.shape[j];

      if (dimMap.has(label)) {
        if (dimMap.get(label) !== size) {
          throw new EinsumError(
            `Dimension mismatch for label '${label}': ${dimMap.get(label)} vs ${size}`,
          );
        }
      } else {
        dimMap.set(label, size);
      }
    }
  }

  return dimMap;
}

/**
 * Execute einsum for a single operand (trace, transpose, diagonal, sum).
 */
async function einsumSingle(
  arr: NDArray,
  inputSub: string,
  outputSub: string,
  _dimMap: Map<string, number>,
): Promise<NDArray> {
  // Find label positions
  const labelPositions = new Map<string, number[]>();
  for (let i = 0; i < inputSub.length; i++) {
    const label = inputSub[i];
    if (!labelPositions.has(label)) {
      labelPositions.set(label, []);
    }
    labelPositions.get(label)!.push(i);
  }

  let result = arr;
  let currentSub = inputSub;

  // Handle repeated indices (trace or diagonal)
  for (const [label, positions] of labelPositions.entries()) {
    if (positions.length > 1) {
      const inOutput = outputSub.includes(label);

      if (!inOutput) {
        // Trace: sum over diagonal
        // For now, handle simple 2D trace case
        if (result.ndim === 2 && positions[0] === 0 && positions[1] === 1) {
          const diag = await diagonalFn(result, 0, 0, 1);
          const diagData = diag.toArray() as number[];
          const sum = diagData.reduce((a, b) => a + b, 0);
          result = await NDArray.fromArray([sum]);
          result = result.reshape([]);
          currentSub = "";
        } else {
          throw new EinsumError(
            `Complex trace not yet supported for indices at positions ${positions}`,
          );
        }
      } else {
        // Diagonal: extract diagonal
        if (result.ndim === 2 && positions[0] === 0 && positions[1] === 1) {
          result = await diagonalFn(result, 0, 0, 1);
          currentSub = label;
        } else {
          throw new EinsumError(`Complex diagonal not yet supported`);
        }
      }
    }
  }

  // Handle sum over axes not in output
  const labelsToSum: string[] = [];
  for (const label of currentSub) {
    if (!outputSub.includes(label)) {
      labelsToSum.push(label);
    }
  }

  for (const label of labelsToSum) {
    const pos = currentSub.indexOf(label);
    if (pos >= 0) {
      // Sum over this axis
      const data = await result.toTypedArray();
      const shape = result.shape;
      const newShape = shape.filter((_, i) => i !== pos);

      if (newShape.length === 0) {
        // Summing to scalar
        const sum = Array.from(data).reduce((a, b) => a + b, 0);
        result = await NDArray.fromArray([sum]);
        result = result.reshape([]);
      } else {
        // Sum along axis
        const axisSize = shape[pos];
        const batchSize = result.size / axisSize;
        const resultData = new Float64Array(batchSize);

        // Move axis to last for easier computation
        const moved = result.moveaxis(pos, -1);
        const movedData = await moved.toTypedArray();

        for (let i = 0; i < batchSize; i++) {
          let sum = 0;
          for (let j = 0; j < axisSize; j++) {
            sum += movedData[i * axisSize + j];
          }
          resultData[i] = sum;
        }

        result = await NDArray.fromTypedArray(
          resultData,
          newShape,
          DType.Float64,
        );
      }
      currentSub = currentSub.replace(label, "");
    }
  }

  // Handle transpose if needed
  if (currentSub.length > 0 && currentSub !== outputSub) {
    const order = outputSub.split("").map((l) => currentSub.indexOf(l));
    if (order.some((i) => i === -1)) {
      throw new EinsumError(`Output label not found in remaining input labels`);
    }
    result = result.transpose(order);
  }

  return result;
}

/**
 * Execute einsum for a pair of operands.
 */
async function einsumPair(
  a: NDArray,
  aSub: string,
  b: NDArray,
  bSub: string,
  outputSub: string,
  _dimMap: Map<string, number>,
): Promise<NDArray> {
  const aLabels = new Set(aSub.split(""));
  const bLabels = new Set(bSub.split(""));

  // Find contracted indices (in both inputs, not in output)
  const contracted = Array.from(aLabels).filter(
    (l) => bLabels.has(l) && !outputSub.includes(l),
  );

  // Find contraction axes
  const axesA = contracted.map((l) => aSub.indexOf(l));
  const axesB = contracted.map((l) => bSub.indexOf(l));

  let result: NDArray;

  if (contracted.length === 0) {
    // Outer product
    result = await tensordot(a, b, 0);
  } else {
    // Use tensordot for contraction
    result = await tensordot(a, b, [axesA, axesB]);
  }

  // Build result subscripts after tensordot
  const remainingA = aSub.split("").filter((l) => !contracted.includes(l));
  const remainingB = bSub.split("").filter((l) => !contracted.includes(l));
  const resultSub = [...remainingA, ...remainingB].join("");

  // Transpose if needed
  if (resultSub !== outputSub && resultSub.length > 0) {
    const order = outputSub.split("").map((l) => resultSub.indexOf(l));
    if (order.some((i) => i === -1)) {
      throw new EinsumError(
        `Output label '${outputSub}' contains labels not in result '${resultSub}'`,
      );
    }
    result = result.transpose(order);
  }

  return result;
}

/**
 * Compute intermediate output subscripts during pairwise reduction.
 */
function intermediateOutput(
  sub1: string,
  sub2: string,
  finalOutput: string,
): string {
  const allLabels = new Set([...sub1.split(""), ...sub2.split("")]);
  const finalLabels = new Set(finalOutput.split(""));

  // Keep labels needed for final output or that appear in sub2 (for next contraction)
  return Array.from(allLabels)
    .filter((l) => finalLabels.has(l) || sub2.includes(l))
    .sort()
    .join("");
}

/**
 * Evaluates the Einstein summation convention on the operands.
 *
 * Using the Einstein summation convention, many common multi-dimensional,
 * linear algebraic array operations can be represented in a simple fashion.
 *
 * @param subscripts - Subscripts string in format "ij,jk->ik"
 * @param operands - Input arrays
 * @returns Result of the einstein summation
 *
 * @example
 * // Matrix multiplication
 * einsum('ij,jk->ik', A, B)
 *
 * // Trace
 * einsum('ii->', A)
 *
 * // Transpose
 * einsum('ij->ji', A)
 *
 * // Sum over axis
 * einsum('ij->i', A)  // sum over j
 *
 * // Outer product
 * einsum('i,j->ij', a, b)
 *
 * // Batch matrix multiply
 * einsum('bij,bjk->bik', A, B)
 */
export async function einsum(
  subscripts: string,
  ...operands: NDArray[]
): Promise<NDArray> {
  // Parse subscripts
  const { inputSubs, outputSubs: parsedOutput } = parseEinsum(subscripts);

  // Validate operand count
  if (operands.length !== inputSubs.length) {
    throw new EinsumError(
      `Number of operands (${operands.length}) doesn't match ` +
        `number of subscript groups (${inputSubs.length})`,
    );
  }

  // Build dimension map
  const dimMap = buildDimMap(inputSubs, operands);

  // Determine output subscripts
  const outputSub = parsedOutput ?? implicitOutput(inputSubs);

  // Validate output subscripts
  const allInputLabels = new Set(inputSubs.join("").split(""));
  for (const label of outputSub) {
    if (!allInputLabels.has(label)) {
      throw new EinsumError(
        `Output label '${label}' not found in input subscripts`,
      );
    }
  }

  // Execute based on number of operands
  if (operands.length === 1) {
    return einsumSingle(operands[0], inputSubs[0], outputSub, dimMap);
  }

  if (operands.length === 2) {
    return einsumPair(
      operands[0],
      inputSubs[0],
      operands[1],
      inputSubs[1],
      outputSub,
      dimMap,
    );
  }

  // Multiple operands: reduce pairwise
  let result = operands[0];
  let resultSubs = inputSubs[0];

  for (let i = 1; i < operands.length; i++) {
    // Compute intermediate output subscripts
    const intermediateSubs =
      i < operands.length - 1
        ? intermediateOutput(resultSubs, inputSubs[i], outputSub)
        : outputSub;

    result = await einsumPair(
      result,
      resultSubs,
      operands[i],
      inputSubs[i],
      intermediateSubs,
      dimMap,
    );
    resultSubs = intermediateSubs;
  }

  return result;
}

/**
 * Estimate contraction cost (product of all involved dimensions).
 */
function contractionCost(
  sub1: string,
  sub2: string,
  dimMap: Map<string, number>,
): number {
  const allLabels = new Set([...sub1.split(""), ...sub2.split("")]);
  let cost = 1;

  for (const label of allLabels) {
    cost *= dimMap.get(label) ?? 1;
  }

  return cost;
}

/**
 * Compute result subscripts after contracting two operands.
 */
function resultSubscripts(
  sub1: string,
  sub2: string,
  finalOutput: string,
): string {
  const labels1 = sub1.split("");
  const labels2 = sub2.split("");

  // Labels in both (contracted) and not in final output
  const contracted = labels1.filter(
    (l) => labels2.includes(l) && !finalOutput.includes(l),
  );

  // Remaining labels (either not contracted or kept in output)
  const remaining = [...labels1, ...labels2].filter(
    (l) => !contracted.includes(l) || finalOutput.includes(l),
  );

  // Remove duplicates, keep order
  const seen = new Set<string>();
  const unique: string[] = [];
  for (const l of remaining) {
    if (!seen.has(l)) {
      seen.add(l);
      unique.push(l);
    }
  }

  return unique.sort().join("");
}

/**
 * Evaluates the lowest cost contraction order for an einsum expression.
 *
 * @param subscripts - Subscripts string
 * @param operands - Input arrays
 * @returns [path, string_representation]
 *
 * @example
 * const [path, info] = einsum_path('ij,jk,kl->il', A, B, C);
 * // path: [[0, 1], [0, 1]]  // First contract A,B then result,C
 */
export function einsum_path(
  subscripts: string,
  ...operands: NDArray[]
): [number[][], string] {
  const { inputSubs, outputSubs: parsedOutput } = parseEinsum(subscripts);

  if (operands.length !== inputSubs.length) {
    throw new EinsumError(
      `Number of operands (${operands.length}) doesn't match ` +
        `number of subscript groups (${inputSubs.length})`,
    );
  }

  // Build dimension map
  const dimMap = buildDimMap(inputSubs, operands);

  // Determine output subscripts
  const outputSub = parsedOutput ?? implicitOutput(inputSubs);

  // Greedy path finding
  const path: number[][] = [];
  const subs = [...inputSubs];
  const indices = inputSubs.map((_, i) => i);

  while (subs.length > 1) {
    // Find best pair to contract (lowest cost)
    let bestPair = [0, 1];
    let bestCost = Infinity;

    for (let i = 0; i < subs.length; i++) {
      for (let j = i + 1; j < subs.length; j++) {
        const cost = contractionCost(subs[i], subs[j], dimMap);
        if (cost < bestCost) {
          bestCost = cost;
          bestPair = [i, j];
        }
      }
    }

    // Record the contraction using original indices
    path.push([indices[bestPair[0]], indices[bestPair[1]]]);

    // Update subscripts and indices
    const newSub = resultSubscripts(
      subs[bestPair[0]],
      subs[bestPair[1]],
      outputSub,
    );
    const newIdx = Math.max(...indices) + 1;

    // Remove contracted (higher index first to preserve lower index)
    subs.splice(bestPair[1], 1);
    subs.splice(bestPair[0], 1);
    indices.splice(bestPair[1], 1);
    indices.splice(bestPair[0], 1);

    // Add result
    subs.push(newSub);
    indices.push(newIdx);
  }

  // Generate info string
  const info = [
    `  Complete contraction:  ${subscripts}`,
    `         Naive scaling:  ${new Set(subscripts.replace(/[^a-zA-Z]/g, "")).size}`,
    `     Optimized scaling:  ${path.length}`,
    ``,
    `  Contraction path:`,
    ...path.map((p, i) => `    ${i}: ${JSON.stringify(p)}`),
  ].join("\n");

  return [path, info];
}
