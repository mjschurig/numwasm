# Monorepo Build System
#
# Build and test commands for numwasm, sciwasm, and symwasm

.PHONY: all clean build build-numwasm build-sciwasm build-symwasm test test-numwasm test-sciwasm test-symwasm install lint format docs help

# Default target
all: build

# Install dependencies
install:
	pnpm install

# Build everything
build:
	pnpm --filter numwasm run build
	pnpm --filter sciwasm run build
	pnpm --filter symwasm run build

# Build individual packages
build-numwasm:
	pnpm --filter numwasm run build

build-sciwasm:
	pnpm --filter numwasm run build:lib
	pnpm --filter sciwasm run build

build-symwasm:
	pnpm --filter symwasm run build

# Run all tests
test: build
	pnpm -r --filter './packages/*' run test

# Test individual packages
test-numwasm: build-numwasm
	pnpm --filter numwasm run test

test-sciwasm: build-sciwasm
	pnpm --filter sciwasm run test

test-symwasm: build-symwasm
	pnpm --filter symwasm run test

# Generate API docs for all packages
docs:
	pnpm -r --filter './packages/*' run docs

# Clean build artifacts
clean:
	pnpm -r run clean

# Lint TypeScript code
lint:
	pnpm -r run lint

# Format code with Prettier
format:
	pnpm -r run format

# Show help
help:
	@echo "Monorepo Build System"
	@echo ""
	@echo "Packages: numwasm, sciwasm, symwasm"
	@echo ""
	@echo "Targets:"
	@echo "  all              Build everything (default)"
	@echo "  install          Install dependencies"
	@echo "  build            Build all packages"
	@echo "  build-numwasm    Build numwasm only"
	@echo "  build-sciwasm    Build sciwasm (and numwasm dep)"
	@echo "  build-symwasm    Build symwasm only"
	@echo "  test             Build and run all tests"
	@echo "  test-numwasm     Test numwasm only"
	@echo "  test-sciwasm     Test sciwasm only"
	@echo "  test-symwasm     Test symwasm only"
	@echo "  docs             Generate API docs for all packages"
	@echo "  clean            Clean build artifacts"
	@echo "  lint             Run ESLint"
	@echo "  format           Format code with Prettier"
	@echo "  help             Show this help"
	@echo ""
	@echo "Quick start:"
	@echo "  make install     # Install dependencies"
	@echo "  make test        # Build and run all tests"
