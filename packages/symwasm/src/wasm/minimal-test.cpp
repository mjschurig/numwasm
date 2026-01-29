/**
 * Minimal test implementation for Phase 1.1 WASM build verification
 * This proves the build system works before integrating full SymEngine
 */

#include <cstdint>
#include <cstring>
#include <cstdlib>

extern "C" {

// Basic memory management test
void* test_malloc(size_t size) {
    return malloc(size);
}

void test_free(void* ptr) {
    free(ptr);
}

// Basic arithmetic test
int64_t test_add(int64_t a, int64_t b) {
    return a + b;
}

int64_t test_multiply(int64_t a, int64_t b) {
    return a * b;
}

// String manipulation test
int test_string_length(const char* str) {
    return str ? strlen(str) : 0;
}

// Simple struct test
struct TestObject {
    int64_t value;
    int type;
};

TestObject* test_create_object(int64_t value, int type) {
    TestObject* obj = (TestObject*)malloc(sizeof(TestObject));
    if (obj) {
        obj->value = value;
        obj->type = type;
    }
    return obj;
}

int64_t test_get_value(TestObject* obj) {
    return obj ? obj->value : 0;
}

void test_free_object(TestObject* obj) {
    free(obj);
}

} // extern "C"
