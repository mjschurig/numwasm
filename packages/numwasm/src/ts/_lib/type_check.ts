/**
 * NumJS Type Checking and Complex Number Utilities
 *
 * Functions for extracting real/imaginary parts and computing angles
 * of complex arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { iscomplexobj } from "../logic.js";
import { arctan2, multiply } from "../ufunc.js";

/**
 * Return the angle of the complex argument.
 *
 * The angle is computed as the counterclockwise angle from the positive
 * real axis, using `arctan2(imag, real)`.
 *
 * @param z - A complex number or array of complex numbers
 * @param deg - Return angle in degrees if true, radians if false (default)
 * @returns The angle in range (-pi, pi] (or (-180, 180] if deg=true)
 *
 * @example
 * ```typescript
 * const c = await NDArray.fromArray([1+1j, 1-1j], { dtype: DType.Complex128 });
 * const rad = await angle(c);      // [0.785, -0.785] (pi/4 radians)
 * const deg = await angle(c, true); // [45, -45] degrees
 * ```
 */
export async function angle(
  z: NDArray,
  deg: boolean = false,
): Promise<NDArray> {
  let zimag: NDArray;
  let zreal: NDArray;

  if (iscomplexobj(z)) {
    zimag = z.imag;
    zreal = z.real;
  } else {
    // For real arrays, imaginary part is 0
    const shape = z.shape;
    const module = z._wasmModule;

    // Create zeros array for imaginary part (synchronously using the module)
    const shapePtr = module._malloc(shape.length * 4);
    for (let i = 0; i < shape.length; i++) {
      module.setValue(shapePtr + i * 4, shape[i], "i32");
    }
    const zerosPtr = module._ndarray_create(shape.length, shapePtr, z.dtype);
    module._free(shapePtr);

    if (zerosPtr === 0) {
      throw new Error("Failed to create zeros array for angle");
    }

    zimag = NDArray._fromPtr(zerosPtr, module);
    zreal = z.copy();
  }

  let a = arctan2(zimag, zreal);

  if (deg) {
    // Multiply by 180/pi using a scalar array that broadcasts
    const scale = await NDArray.full([1], 180 / Math.PI, { dtype: a.dtype });
    a = multiply(a, scale);
  }

  return a;
}

/**
 * Return the real part of the complex argument.
 *
 * For non-complex arrays, returns a copy of the input.
 *
 * @param val - Input array
 * @returns The real component of the complex argument
 *
 * @example
 * ```typescript
 * const c = await NDArray.fromArray([1+2j, 3+4j], { dtype: DType.Complex128 });
 * const r = real(c);  // [1, 3]
 * ```
 */
export function real(val: NDArray): NDArray {
  return val.real;
}

/**
 * Return the imaginary part of the complex argument.
 *
 * For non-complex arrays, returns an array of zeros with the same shape.
 *
 * @param val - Input array
 * @returns The imaginary component of the complex argument
 *
 * @example
 * ```typescript
 * const c = await NDArray.fromArray([1+2j, 3+4j], { dtype: DType.Complex128 });
 * const i = imag(c);  // [2, 4]
 * ```
 */
export function imag(val: NDArray): NDArray {
  return val.imag;
}
