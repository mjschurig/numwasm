const createSymwasmModule = require('./dist/wasm/symwasm.cjs');

(async () => {
  console.log('Loading WASM module directly...');
  
  try {
    const wasm = await createSymwasmModule();
    console.log('✓ WASM module loaded!');
    
    // Test basic memory allocation
    const ptr = wasm._basic_new_stack();
    console.log(`✓ Allocated Basic at pointer: ${ptr}`);

    // Test integer creation
    wasm._integer_set_si(ptr, 42);
    const strPtr = wasm._basic_str(ptr);
    const str = wasm.UTF8ToString(strPtr);
    console.log(`✓ Integer value: "${str}"`);

    wasm._free(strPtr);
    wasm._basic_free_stack(ptr);
    console.log('✓ Freed memory\n');
    console.log('=== GMP-powered SymEngine WASM is working! ===');
  } catch (err) {
    console.error('Error:', err);
    process.exit(1);
  }
})();
