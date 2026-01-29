const { core } = require('./dist/symwasm.cjs');
const { loadWasmModule } = core;

(async () => {
  console.log('Loading SymWASM module with GMP support...');

  try {
    const wasm = await loadWasmModule();
    console.log('✓ WASM module loaded successfully!');

    // Test basic memory allocation
    const ptr = wasm._basic_new_stack();
    console.log(`✓ Allocated Basic object at pointer: ${ptr}`);

    // Test GMP-based integer creation
    const result = wasm._integer_set_si(ptr, 42);
    console.log(`✓ Created integer 42, result code: ${result}`);

    // Test string conversion
    const strPtr = wasm._basic_str(ptr);
    const str = wasm.UTF8ToString(strPtr);
    console.log(`✓ String representation: "${str}"`);

    wasm._free(strPtr);
    wasm._basic_free_stack(ptr);
    console.log('✓ Freed memory');

    console.log('\n=== GMP-powered SymWASM is working! ===');
  } catch (err) {
    console.error('Error:', err);
    process.exit(1);
  }
})();
