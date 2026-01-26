# Phase 28: Complex Number Accessor Functions

## Overview

Implement accessor functions for complex-valued arrays. While Complex64 and Complex128 dtypes are already supported, the standard NumPy accessor functions are missing.

## Functions to Implement

### 1. `real(x)` - Extract real part
```typescript
function real(x: NDArray): NDArray
```
- Returns the real part of complex numbers
- For real-valued input, returns the input unchanged (or a view)
- Should work element-wise on arrays of any shape

**NumPy Reference:** `numpy/_core/fromnumeric.py` and `numpy/_core/src/multiarray/`

### 2. `imag(x)` - Extract imaginary part
```typescript
function imag(x: NDArray): NDArray
```
- Returns the imaginary part of complex numbers
- For real-valued input, returns zeros with same shape
- Should work element-wise on arrays of any shape

### 3. `conj(x)` / `conjugate(x)` - Complex conjugate
```typescript
function conj(x: NDArray): NDArray
function conjugate(x: NDArray): NDArray  // alias
```
- Returns the complex conjugate (a + bi → a - bi)
- For real-valued input, returns the input unchanged
- Should be a ufunc for broadcasting support

### 4. `angle(z, deg?)` - Phase angle
```typescript
function angle(z: NDArray, deg: boolean = false): NDArray
```
- Returns the angle (argument) of complex numbers in radians (or degrees if deg=true)
- Uses `atan2(imag(z), real(z))` internally
- Returns 0 for real positive numbers, π for real negative

## Implementation Strategy

### File Location
- Add to `src/ts/complex.ts` (new file)
- Export from `src/ts/index.ts`

### Dependencies
- Existing complex dtype support in `src/ts/types.ts`
- Complex data stored as interleaved [real, imag, real, imag, ...] in Float32Array/Float64Array

### WASM Acceleration (Optional)
- Could add WASM implementations for large arrays
- TypeScript implementation sufficient for initial version

## Implementation Details

### Accessing Complex Data

Complex arrays store data as interleaved pairs:
```typescript
// For Complex64: Float32Array with pairs [re0, im0, re1, im1, ...]
// For Complex128: Float64Array with pairs [re0, im0, re1, im1, ...]
```

### real() Implementation
```typescript
export function real(x: NDArray): NDArray {
  if (!isComplexDtype(x.dtype)) {
    return x;  // or x.copy() for consistency
  }

  const result = zeros(x.shape, getRealDtype(x.dtype));
  const data = x.data;
  const outData = result.data;

  for (let i = 0; i < x.size; i++) {
    outData[i] = data[i * 2];  // Real part at even indices
  }

  return result;
}
```

### imag() Implementation
```typescript
export function imag(x: NDArray): NDArray {
  if (!isComplexDtype(x.dtype)) {
    return zeros(x.shape, x.dtype);
  }

  const result = zeros(x.shape, getRealDtype(x.dtype));
  const data = x.data;
  const outData = result.data;

  for (let i = 0; i < x.size; i++) {
    outData[i] = data[i * 2 + 1];  // Imag part at odd indices
  }

  return result;
}
```

### conj() Implementation
```typescript
export function conj(x: NDArray): NDArray {
  if (!isComplexDtype(x.dtype)) {
    return x.copy();
  }

  const result = empty(x.shape, x.dtype);
  const data = x.data;
  const outData = result.data;

  for (let i = 0; i < x.size; i++) {
    outData[i * 2] = data[i * 2];      // Real part unchanged
    outData[i * 2 + 1] = -data[i * 2 + 1];  // Negate imaginary part
  }

  return result;
}

export const conjugate = conj;  // Alias
```

### angle() Implementation
```typescript
export function angle(z: NDArray, deg: boolean = false): NDArray {
  const re = real(z);
  const im = imag(z);

  // Use atan2 for proper quadrant handling
  const result = atan2(im, re);

  if (deg) {
    return multiply(result, 180 / Math.PI);
  }

  return result;
}
```

## Helper Functions Needed

```typescript
function isComplexDtype(dtype: DType): boolean {
  return dtype === 'complex64' || dtype === 'complex128';
}

function getRealDtype(complexDtype: DType): DType {
  return complexDtype === 'complex64' ? 'float32' : 'float64';
}
```

## Testing

Create `tests/ts/complex.test.ts`:

```typescript
describe('Complex accessors', () => {
  describe('real', () => {
    it('extracts real part from complex array', async () => {
      const z = createComplexArray([1, 2, 3, 4], [5, 6, 7, 8], 'complex128');
      const r = real(z);
      expect(await r.toArray()).toEqual([1, 2, 3, 4]);
    });

    it('returns input for real arrays', () => {
      const x = fromArray([1, 2, 3]);
      const r = real(x);
      expect(r.dtype).toBe('float64');
    });
  });

  describe('imag', () => {
    it('extracts imaginary part from complex array', async () => {
      const z = createComplexArray([1, 2], [3, 4], 'complex128');
      const i = imag(z);
      expect(await i.toArray()).toEqual([3, 4]);
    });

    it('returns zeros for real arrays', async () => {
      const x = fromArray([1, 2, 3]);
      const i = imag(x);
      expect(await i.toArray()).toEqual([0, 0, 0]);
    });
  });

  describe('conj', () => {
    it('conjugates complex array', async () => {
      const z = createComplexArray([1, 2], [3, -4], 'complex128');
      const c = conj(z);
      // Real parts: [1, 2], Imag parts: [-3, 4]
    });
  });

  describe('angle', () => {
    it('computes angle in radians', () => {
      // 1+1i has angle π/4
      // -1+0i has angle π
      // 1+0i has angle 0
    });

    it('computes angle in degrees when deg=true', () => {
      // 1+1i has angle 45°
    });
  });
});
```

## Exports to Add

In `src/ts/index.ts`:
```typescript
export { real, imag, conj, conjugate, angle } from './complex';
```

## Priority

Medium - Complex arrays are supported but these accessors are commonly needed for FFT results and signal processing workflows.

## Estimated Scope

- ~150 lines of TypeScript
- ~100 lines of tests
- No WASM changes required initially
