# NumJS Makefile
#
# Build and test commands for the NumPy-to-TypeScript/WASM project

.PHONY: all clean build build-wasm build-ts test test-unit test-compare generate-fixtures install lint format help benchmark benchmark-numpy benchmark-numjs benchmark-report

# Default target
all: build

# Install dependencies
install:
	npm install

# Build everything
build: build-wasm build-lib

# Build WASM module
build-wasm:
	@echo "Building WASM module..."
	./scripts/build-wasm.sh

# Build library with Vite
build-lib:
	@echo "Building library with Vite..."
	npm run build:lib

# Generate test fixtures from Python/NumPy
generate-fixtures:
	@echo "Generating test fixtures with NumPy..."
	python3 tests/python/generate_test_cases.py

# Run all tests
test: generate-fixtures build
	npm test

# Run unit tests only
test-unit: build
	npm run test -- tests/ts/ndarray.test.ts

# Run comparison tests (requires fixtures)
test-compare: generate-fixtures build
	npm run test -- tests/ts/comparison.test.ts

# Run tests in watch mode
test-watch: build
	npm run test:watch

# Clean build artifacts
clean:
	rm -rf dist/
	rm -rf tests/fixtures/*.json
	rm -rf node_modules/.vitest

# Full clean including node_modules
clean-all: clean
	rm -rf node_modules/

# Lint TypeScript code
lint:
	npm run lint

# Format code with Prettier
format:
	npm run format

# Benchmark targets
benchmark: benchmark-numpy build benchmark-numjs benchmark-report
	@echo ""
	@echo "Benchmark complete! Open benchmark/index.html to view results."

benchmark-numpy:
	@echo "Running NumPy benchmark..."
	python3 benchmark/python/benchmark_sum.py

benchmark-numjs: build
	@echo "Building benchmark TypeScript..."
	npm run build:benchmark
	@echo "Running NumJS benchmark..."
	node dist/benchmark/ts/benchmark_sum.js

benchmark-report:
	@echo "Generating HTML report..."
	node dist/benchmark/report/generate_report.js

# Show help
help:
	@echo "NumJS Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all              Build everything (default)"
	@echo "  install          Install npm dependencies"
	@echo "  build            Build WASM and TypeScript"
	@echo "  build-wasm       Build only WASM module"
	@echo "  build-ts         Build only TypeScript"
	@echo "  generate-fixtures Generate test vectors with NumPy"
	@echo "  test             Run all tests"
	@echo "  test-unit        Run unit tests only"
	@echo "  test-compare     Run NumPy comparison tests"
	@echo "  test-watch       Run tests in watch mode"
	@echo "  benchmark        Run full benchmark suite"
	@echo "  benchmark-numpy  Run NumPy benchmark only"
	@echo "  benchmark-numjs  Run NumJS benchmark only"
	@echo "  benchmark-report Generate HTML report"
	@echo "  clean            Clean build artifacts"
	@echo "  clean-all        Clean everything including node_modules"
	@echo "  lint             Run ESLint"
	@echo "  format           Format code with Prettier"
	@echo "  help             Show this help"
	@echo ""
	@echo "Quick start:"
	@echo "  make install     # Install dependencies"
	@echo "  make test        # Build and run all tests"
