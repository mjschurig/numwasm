/**
 * Base Conversion Functions
 *
 * Functions for converting numbers to binary and arbitrary base representations.
 *
 * Reference: numpy/_core/numeric.py (binary_repr, base_repr)
 */

/**
 * Return the binary representation of a number as a string.
 *
 * For negative numbers, uses two's complement representation if width is provided,
 * otherwise returns the binary representation with a minus sign.
 *
 * @param num - Integer to convert
 * @param width - Fixed width for two's complement representation
 * @returns Binary string representation
 *
 * @example
 * binaryRepr(5)       // '101'
 * binaryRepr(5, 8)    // '00000101'
 * binaryRepr(-5)      // '-101'
 * binaryRepr(-5, 8)   // '11111011'
 *
 * // Two's complement for negative numbers
 * binaryRepr(-1, 8)   // '11111111'
 * binaryRepr(-128, 8) // '10000000'
 */
export function binaryRepr(num: number, width?: number): string {
  // Ensure integer
  const intNum = Math.trunc(num);

  if (intNum >= 0) {
    // Positive number: simple binary representation
    const binary = intNum.toString(2);

    if (width !== undefined) {
      if (binary.length > width) {
        throw new Error(
          `Width ${width} too small to represent ${intNum} in binary`
        );
      }
      return binary.padStart(width, '0');
    }

    return binary;
  }

  // Negative number
  if (width === undefined) {
    // Without width, return negative sign representation
    return '-' + Math.abs(intNum).toString(2);
  }

  // Two's complement representation
  // For width w, valid range is -2^(w-1) to 2^(w-1) - 1
  const minVal = -(1 << (width - 1));
  const maxVal = (1 << (width - 1)) - 1;

  if (intNum < minVal || intNum > maxVal) {
    throw new Error(
      `Width ${width} too small to represent ${intNum} in two's complement`
    );
  }

  // Two's complement: add 2^width to negative number
  const twosComplement = (1 << width) + intNum;
  return twosComplement.toString(2).padStart(width, '0');
}

/**
 * Return a string representation of a number in the given base.
 *
 * @param num - Integer to convert
 * @param base - Number base (2-36)
 * @param padding - Minimum number of digits
 * @returns String representation in base
 *
 * @example
 * baseRepr(10, 2)      // '1010'
 * baseRepr(10, 8)      // '12'
 * baseRepr(10, 16)     // 'A'
 * baseRepr(255, 16)    // 'FF'
 * baseRepr(10, 2, 8)   // '00001010'
 * baseRepr(-10, 16)    // '-A'
 * baseRepr(36, 36)     // '10'
 */
export function baseRepr(num: number, base: number = 10, padding: number = 0): string {
  // Validate base
  if (base < 2 || base > 36) {
    throw new Error(`Base must be between 2 and 36, got ${base}`);
  }

  // Ensure integer
  const intNum = Math.trunc(num);

  // Handle negative numbers
  const isNegative = intNum < 0;
  const absNum = Math.abs(intNum);

  // Convert to base
  let result = absNum.toString(base).toUpperCase();

  // Apply padding
  if (padding > 0) {
    result = result.padStart(padding, '0');
  }

  // Add sign for negative
  if (isNegative) {
    result = '-' + result;
  }

  return result;
}

/**
 * Convert a binary string to a number.
 *
 * @param binary - Binary string (optionally with '0b' prefix)
 * @returns Parsed integer
 *
 * @example
 * fromBinaryRepr('101')        // 5
 * fromBinaryRepr('0b101')      // 5
 * fromBinaryRepr('-101')       // -5
 * fromBinaryRepr('11111011', true, 8)  // -5 (two's complement)
 */
export function fromBinaryRepr(
  binary: string,
  twosComplement: boolean = false,
  width?: number
): number {
  // Remove '0b' prefix if present
  let str = binary.toLowerCase().replace(/^0b/, '');

  // Handle negative sign
  const isNegative = str.startsWith('-');
  if (isNegative) {
    str = str.slice(1);
  }

  // Parse as positive binary
  const value = parseInt(str, 2);

  if (isNaN(value)) {
    throw new Error(`Invalid binary string: ${binary}`);
  }

  if (isNegative) {
    return -value;
  }

  // Check for two's complement
  if (twosComplement && width !== undefined && str.length === width) {
    // Check if MSB is 1 (negative in two's complement)
    if (str[0] === '1') {
      // Convert from two's complement
      return value - (1 << width);
    }
  }

  return value;
}

/**
 * Convert a base-N string to a number.
 *
 * @param str - String representation
 * @param base - Number base (2-36)
 * @returns Parsed integer
 *
 * @example
 * fromBaseRepr('1010', 2)   // 10
 * fromBaseRepr('12', 8)     // 10
 * fromBaseRepr('A', 16)     // 10
 * fromBaseRepr('FF', 16)    // 255
 * fromBaseRepr('-A', 16)    // -10
 */
export function fromBaseRepr(str: string, base: number = 10): number {
  if (base < 2 || base > 36) {
    throw new Error(`Base must be between 2 and 36, got ${base}`);
  }

  const value = parseInt(str, base);

  if (isNaN(value)) {
    throw new Error(`Invalid base-${base} string: ${str}`);
  }

  return value;
}

/**
 * Return the hexadecimal representation of a number.
 *
 * @param num - Integer to convert
 * @param width - Minimum number of hex digits
 * @returns Hexadecimal string (uppercase)
 *
 * @example
 * hexRepr(255)      // 'FF'
 * hexRepr(255, 4)   // '00FF'
 * hexRepr(-1)       // '-1'
 * hexRepr(16)       // '10'
 */
export function hexRepr(num: number, width: number = 0): string {
  return baseRepr(num, 16, width);
}

/**
 * Return the octal representation of a number.
 *
 * @param num - Integer to convert
 * @param width - Minimum number of octal digits
 * @returns Octal string
 *
 * @example
 * octalRepr(8)      // '10'
 * octalRepr(64)     // '100'
 * octalRepr(7)      // '7'
 */
export function octalRepr(num: number, width: number = 0): string {
  return baseRepr(num, 8, width);
}

/**
 * Convert an array of numbers to binary representations.
 *
 * @param arr - Array of integers
 * @param width - Fixed width for all representations
 * @returns Array of binary strings
 *
 * @example
 * binaryReprArray([1, 2, 3, 4])     // ['1', '10', '11', '100']
 * binaryReprArray([1, 2, 3, 4], 4)  // ['0001', '0010', '0011', '0100']
 */
export function binaryReprArray(arr: number[], width?: number): string[] {
  return arr.map(num => binaryRepr(num, width));
}

/**
 * Convert an array of numbers to base-N representations.
 *
 * @param arr - Array of integers
 * @param base - Number base (2-36)
 * @param padding - Minimum digits per number
 * @returns Array of base-N strings
 *
 * @example
 * baseReprArray([10, 20, 30], 16)  // ['A', '14', '1E']
 */
export function baseReprArray(
  arr: number[],
  base: number = 10,
  padding: number = 0
): string[] {
  return arr.map(num => baseRepr(num, base, padding));
}
