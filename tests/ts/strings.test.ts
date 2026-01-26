/**
 * Tests for NumJS String Operations Module (Phase 16a)
 */

import { describe, it, expect, beforeAll } from 'vitest';
import {
  NDArray,
  DType,
  loadWasmModule,
  strings,
  ValueError,
} from '../../dist/numjs.mjs';

// Load WASM module before tests
beforeAll(async () => {
  await loadWasmModule();
});

describe('String Array Infrastructure', () => {
  describe('NDArray.fromStringArray', () => {
    it('creates a 1D string array', () => {
      const arr = NDArray.fromStringArray(['hello', 'world']);
      expect(arr.shape).toEqual([2]);
      expect(arr.dtype).toBe(DType.String);
      expect(arr.isStringArray).toBe(true);
      expect(arr.size).toBe(2);
    });

    it('creates a 2D string array', () => {
      const arr = NDArray.fromStringArray([
        ['a', 'b'],
        ['c', 'd'],
      ]);
      expect(arr.shape).toEqual([2, 2]);
      expect(arr.size).toBe(4);
    });

    it('creates string array with explicit shape', () => {
      const arr = NDArray.fromStringArray(['a', 'b', 'c', 'd'], [2, 2]);
      expect(arr.shape).toEqual([2, 2]);
    });

    it('throws on shape mismatch', () => {
      expect(() => {
        NDArray.fromStringArray(['a', 'b', 'c'], [2, 2]);
      }).toThrow();
    });
  });

  describe('NDArray.emptyString', () => {
    it('creates empty string array', () => {
      const arr = NDArray.emptyString([3]);
      expect(arr.shape).toEqual([3]);
      expect(arr.isStringArray).toBe(true);
      expect(arr.getStringFlat(0)).toBe('');
    });
  });

  describe('String element access', () => {
    it('getStringFlat returns correct values', () => {
      const arr = NDArray.fromStringArray(['hello', 'world']);
      expect(arr.getStringFlat(0)).toBe('hello');
      expect(arr.getStringFlat(1)).toBe('world');
    });

    it('setStringFlat modifies values', () => {
      const arr = NDArray.fromStringArray(['hello', 'world']);
      arr.setStringFlat(0, 'hi');
      expect(arr.getStringFlat(0)).toBe('hi');
    });

    it('toArray returns string array', () => {
      const arr = NDArray.fromStringArray(['hello', 'world']);
      expect(arr.toArray()).toEqual(['hello', 'world']);
    });

    it('copy creates independent copy', () => {
      const arr = NDArray.fromStringArray(['hello', 'world']);
      const copy = arr.copy();
      copy.setStringFlat(0, 'hi');
      expect(arr.getStringFlat(0)).toBe('hello');
      expect(copy.getStringFlat(0)).toBe('hi');
    });
  });
});

describe('Comparison Functions', () => {
  describe('equal', () => {
    it('compares string arrays element-wise', () => {
      const result = strings.equal(['a', 'b'], ['a', 'c']);
      expect(result.toArray()).toEqual([1, 0]);
    });

    it('works with NDArray inputs', () => {
      const a = NDArray.fromStringArray(['hello', 'world']);
      const b = NDArray.fromStringArray(['hello', 'test']);
      const result = strings.equal(a, b);
      expect(result.toArray()).toEqual([1, 0]);
    });
  });

  describe('not_equal', () => {
    it('returns true where strings differ', () => {
      const result = strings.not_equal(['a', 'b'], ['a', 'c']);
      expect(result.toArray()).toEqual([0, 1]);
    });
  });

  describe('less', () => {
    it('compares strings lexicographically', () => {
      const result = strings.less(['a', 'b'], ['b', 'a']);
      expect(result.toArray()).toEqual([1, 0]);
    });
  });

  describe('less_equal', () => {
    it('includes equality', () => {
      const result = strings.less_equal(['a', 'b'], ['a', 'a']);
      expect(result.toArray()).toEqual([1, 0]);
    });
  });

  describe('greater', () => {
    it('compares strings lexicographically', () => {
      const result = strings.greater(['b', 'a'], ['a', 'b']);
      expect(result.toArray()).toEqual([1, 0]);
    });
  });

  describe('greater_equal', () => {
    it('includes equality', () => {
      const result = strings.greater_equal(['a', 'b'], ['a', 'c']);
      expect(result.toArray()).toEqual([1, 0]);
    });
  });

  describe('compare_chararrays', () => {
    it('returns -1, 0, 1', () => {
      const result = strings.compare_chararrays(['a', 'b', 'c'], ['b', 'b', 'a']);
      expect(result.toArray()).toEqual([-1, 0, 1]);
    });
  });
});

describe('Property Testing Functions', () => {
  describe('isalpha', () => {
    it('returns true for alphabetic strings', () => {
      const result = strings.isalpha(['hello', '123', 'abc123', '']);
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });
  });

  describe('isdigit', () => {
    it('returns true for digit strings', () => {
      const result = strings.isdigit(['123', 'abc', '12.3', '']);
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });
  });

  describe('isalnum', () => {
    it('returns true for alphanumeric strings', () => {
      const result = strings.isalnum(['abc123', 'abc', '123', 'abc-123']);
      expect(result.toArray()).toEqual([1, 1, 1, 0]);
    });
  });

  describe('isspace', () => {
    it('returns true for whitespace strings', () => {
      const result = strings.isspace(['   ', '\t\n', 'a b', '']);
      expect(result.toArray()).toEqual([1, 1, 0, 0]);
    });
  });

  describe('islower', () => {
    it('returns true for lowercase strings', () => {
      const result = strings.islower(['hello', 'Hello', 'HELLO', '123']);
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });
  });

  describe('isupper', () => {
    it('returns true for uppercase strings', () => {
      const result = strings.isupper(['HELLO', 'Hello', 'hello', '123']);
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });
  });

  describe('istitle', () => {
    it('returns true for titlecased strings', () => {
      const result = strings.istitle(['Hello World', 'hello world', 'HELLO', '123']);
      expect(result.toArray()).toEqual([1, 0, 0, 0]);
    });
  });

  describe('str_len', () => {
    it('returns string lengths', () => {
      const result = strings.str_len(['hello', 'hi', '']);
      expect(result.toArray()).toEqual([5, 2, 0]);
    });
  });
});

describe('Search Functions', () => {
  describe('find', () => {
    it('returns lowest index of substring', () => {
      const result = strings.find(['hello world', 'test'], 'o');
      expect(result.toArray()).toEqual([4, -1]);
    });

    it('respects start parameter', () => {
      const result = strings.find(['hello world'], 'o', 5);
      expect(result.toArray()).toEqual([7]);
    });
  });

  describe('rfind', () => {
    it('returns highest index of substring', () => {
      const result = strings.rfind(['hello world'], 'o');
      expect(result.toArray()).toEqual([7]);
    });
  });

  describe('index', () => {
    it('returns index like find', () => {
      const result = strings.index(['hello'], 'l');
      expect(result.toArray()).toEqual([2]);
    });

    it('throws ValueError when not found', () => {
      expect(() => {
        strings.index(['hello'], 'x');
      }).toThrow(ValueError);
    });
  });

  describe('count', () => {
    it('counts non-overlapping occurrences', () => {
      const result = strings.count(['aabababa', 'foo'], 'aba');
      expect(result.toArray()).toEqual([2, 0]);
    });

    it('handles empty substring', () => {
      const result = strings.count(['hello'], '');
      expect(result.toArray()).toEqual([6]); // length + 1
    });
  });

  describe('startswith', () => {
    it('checks prefix', () => {
      const result = strings.startswith(['hello', 'world', 'help'], 'hel');
      expect(result.toArray()).toEqual([1, 0, 1]);
    });
  });

  describe('endswith', () => {
    it('checks suffix', () => {
      const result = strings.endswith(['hello', 'world', 'jello'], 'llo');
      expect(result.toArray()).toEqual([1, 0, 1]);
    });
  });
});

describe('Manipulation Functions', () => {
  describe('upper', () => {
    it('converts to uppercase', () => {
      const result = strings.upper(['hello', 'World']);
      expect(result.toArray()).toEqual(['HELLO', 'WORLD']);
    });
  });

  describe('lower', () => {
    it('converts to lowercase', () => {
      const result = strings.lower(['HELLO', 'World']);
      expect(result.toArray()).toEqual(['hello', 'world']);
    });
  });

  describe('swapcase', () => {
    it('swaps case', () => {
      const result = strings.swapcase(['Hello', 'WORLD']);
      expect(result.toArray()).toEqual(['hELLO', 'world']);
    });
  });

  describe('capitalize', () => {
    it('capitalizes first char', () => {
      const result = strings.capitalize(['hello', 'WORLD']);
      expect(result.toArray()).toEqual(['Hello', 'World']);
    });
  });

  describe('title', () => {
    it('titlecases words', () => {
      const result = strings.title(['hello world', 'HELLO WORLD']);
      expect(result.toArray()).toEqual(['Hello World', 'Hello World']);
    });
  });

  describe('add', () => {
    it('concatenates strings', () => {
      const result = strings.add(['hello', 'good'], [' world', 'bye']);
      expect(result.toArray()).toEqual(['hello world', 'goodbye']);
    });
  });

  describe('multiply', () => {
    it('repeats strings', () => {
      const result = strings.multiply(['ab', 'cd'], 3);
      expect(result.toArray()).toEqual(['ababab', 'cdcdcd']);
    });

    it('handles zero repetitions', () => {
      const result = strings.multiply(['ab'], 0);
      expect(result.toArray()).toEqual(['']);
    });
  });

  describe('strip', () => {
    it('removes whitespace by default', () => {
      const result = strings.strip(['  hi  ', '\tfoo\n']);
      expect(result.toArray()).toEqual(['hi', 'foo']);
    });

    it('removes specified chars', () => {
      const result = strings.strip(['xxhelloxx'], 'x');
      expect(result.toArray()).toEqual(['hello']);
    });
  });

  describe('lstrip', () => {
    it('removes leading whitespace', () => {
      const result = strings.lstrip(['  hi  ']);
      expect(result.toArray()).toEqual(['hi  ']);
    });
  });

  describe('rstrip', () => {
    it('removes trailing whitespace', () => {
      const result = strings.rstrip(['  hi  ']);
      expect(result.toArray()).toEqual(['  hi']);
    });
  });

  describe('replace', () => {
    it('replaces all occurrences by default', () => {
      const result = strings.replace(['aabbcc'], 'b', 'x');
      expect(result.toArray()).toEqual(['aaxxcc']);
    });

    it('respects count parameter', () => {
      const result = strings.replace(['aabbcc'], 'b', 'x', 1);
      expect(result.toArray()).toEqual(['aaxbcc']);
    });
  });

  describe('center', () => {
    it('centers string', () => {
      const result = strings.center(['ab', 'x'], 5);
      expect(result.toArray()).toEqual([' ab  ', '  x  ']);
    });

    it('uses custom fillchar', () => {
      const result = strings.center(['ab'], 5, '*');
      expect(result.toArray()).toEqual(['*ab**']);
    });
  });

  describe('ljust', () => {
    it('left justifies', () => {
      const result = strings.ljust(['ab'], 5);
      expect(result.toArray()).toEqual(['ab   ']);
    });
  });

  describe('rjust', () => {
    it('right justifies', () => {
      const result = strings.rjust(['ab'], 5);
      expect(result.toArray()).toEqual(['   ab']);
    });
  });

  describe('zfill', () => {
    it('pads with zeros', () => {
      const result = strings.zfill(['42', '7'], 5);
      expect(result.toArray()).toEqual(['00042', '00007']);
    });

    it('handles sign', () => {
      const result = strings.zfill(['-42', '+42'], 5);
      expect(result.toArray()).toEqual(['-0042', '+0042']);
    });
  });

  describe('partition', () => {
    it('splits on first occurrence', () => {
      const result = strings.partition(['a-b-c', 'foo'], '-');
      expect(result.toArray()).toEqual(['a', '-', 'b-c', 'foo', '', '']);
    });
  });

  describe('rpartition', () => {
    it('splits on last occurrence', () => {
      const result = strings.rpartition(['a-b-c'], '-');
      expect(result.toArray()).toEqual(['a-b', '-', 'c']);
    });
  });

  describe('expandtabs', () => {
    it('expands tabs to spaces', () => {
      const result = strings.expandtabs(['a\tb'], 4);
      expect(result.toArray()).toEqual(['a   b']);
    });
  });

  describe('encode/decode', () => {
    it('encodes strings to bytes', () => {
      const encoded = strings.encode(['hello']);
      expect(encoded[0]).toBeInstanceOf(Uint8Array);
      expect(encoded[0].length).toBe(5);
    });

    it('decodes bytes to strings', () => {
      const encoder = new TextEncoder();
      const bytes = [encoder.encode('hello'), encoder.encode('world')];
      const decoded = strings.decode(bytes);
      expect(decoded.toArray()).toEqual(['hello', 'world']);
    });
  });
});
