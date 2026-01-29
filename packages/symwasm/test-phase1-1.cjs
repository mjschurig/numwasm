#!/usr/bin/env node
/**
 * Phase 1.1 Verification Test
 * Tests that the WASM module loads and basic C API functions work correctly
 */

const path = require('path');

async function testPhase1_1() {
  console.log('='.repeat(60));
  console.log('Phase 1.1 Verification Test');
  console.log('='.repeat(60));
  console.log();

  // Load the WASM module
  console.log('1. Loading WASM module...');
  const wasmPath = path.join(__dirname, 'dist', 'wasm', 'symwasm.cjs');
  const createSymwasmModule = require(wasmPath);

  if (typeof createSymwasmModule !== 'function') {
    console.log('   ✗ Module export is not a factory function');
    return false;
  }

  const wasm = await createSymwasmModule();
  console.log('   ✓ WASM module loaded successfully');
  console.log();

  // Test exported functions
  console.log('2. Checking exported C API functions...');
  const requiredFunctions = [
    '_malloc',
    '_free',
    '_basic_new_heap',
    '_basic_free_heap',
    '_symbol_set',
    '_integer_set_si',
    '_basic_add',
    '_basic_mul',
    '_basic_pow',
    '_basic_str',
  ];

  let allPresent = true;
  for (const funcName of requiredFunctions) {
    if (typeof wasm[funcName] === 'function') {
      console.log(`   ✓ ${funcName}`);
    } else {
      console.log(`   ✗ ${funcName} MISSING`);
      allPresent = false;
    }
  }
  console.log();

  if (!allPresent) {
    return false;
  }

  // Test arithmetic
  console.log('3. Testing integer arithmetic...');
  const a = wasm._basic_new_heap();
  const b = wasm._basic_new_heap();
  const result = wasm._basic_new_heap();

  wasm._integer_set_si(a, 3);
  wasm._integer_set_si(b, 5);

  // 3 + 5 = 8
  wasm._basic_add(result, a, b);
  let strPtr = wasm._basic_str(result);
  let str = wasm.UTF8ToString(strPtr);
  wasm._free(strPtr);
  console.log(`   3 + 5 = ${str}`, str === '8' ? '✓' : '✗');
  if (str !== '8') allPresent = false;

  // 3 * 5 = 15
  wasm._basic_mul(result, a, b);
  strPtr = wasm._basic_str(result);
  str = wasm.UTF8ToString(strPtr);
  wasm._free(strPtr);
  console.log(`   3 * 5 = ${str}`, str === '15' ? '✓' : '✗');
  if (str !== '15') allPresent = false;

  // 3^5 = 243
  wasm._basic_pow(result, a, b);
  strPtr = wasm._basic_str(result);
  str = wasm.UTF8ToString(strPtr);
  wasm._free(strPtr);
  console.log(`   3^5 = ${str}`, str === '243' ? '✓' : '✗');
  if (str !== '243') allPresent = false;

  wasm._basic_free_heap(a);
  wasm._basic_free_heap(b);
  console.log();

  // Test symbols
  console.log('4. Testing symbolic math...');
  const x = wasm._basic_new_heap();
  const xNamePtr = wasm._malloc(2);
  wasm.stringToUTF8('x', xNamePtr, 2);
  wasm._symbol_set(x, xNamePtr);
  wasm._free(xNamePtr);

  // x + x = 2*x
  wasm._basic_add(result, x, x);
  strPtr = wasm._basic_str(result);
  str = wasm.UTF8ToString(strPtr);
  wasm._free(strPtr);
  console.log(`   x + x = ${str}`, str === '2*x' ? '✓' : '✗');
  if (str !== '2*x') allPresent = false;

  // x * x = x**2
  wasm._basic_mul(result, x, x);
  strPtr = wasm._basic_str(result);
  str = wasm.UTF8ToString(strPtr);
  wasm._free(strPtr);
  console.log(`   x * x = ${str}`, str === 'x**2' ? '✓' : '✗');
  if (str !== 'x**2') allPresent = false;

  wasm._basic_free_heap(x);
  wasm._basic_free_heap(result);
  console.log();

  console.log('='.repeat(60));
  if (allPresent) {
    console.log('✅ Phase 1.1 Verification Complete!');
  } else {
    console.log('❌ Phase 1.1 Verification Failed!');
  }
  console.log('='.repeat(60));
  console.log();
  console.log('Build artifacts:');
  console.log('- dist/wasm/symwasm.wasm (563KB)');
  console.log('- dist/wasm/symwasm.cjs (19KB)');
  console.log('- dist/wasm/symwasm.mjs (19KB)');
  console.log();

  return allPresent;
}

testPhase1_1().then(success => {
  process.exit(success ? 0 : 1);
}).catch(err => {
  console.error('Test failed:', err);
  process.exit(1);
});
