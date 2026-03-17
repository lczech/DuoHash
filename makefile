.PHONY: build clean

# All source files that should trigger a rebuild when changed
SOURCES := $(shell find src include example -type f \( -name "*.cpp" -o -name "*.cc" -o -name "*.c" -o -name "*.h" -o -name "*.hpp" \))

build: build/libDuoHash.a

# Re-run cmake + make if any source/header changes
build/libDuoHash.a: $(SOURCES) CMakeLists.txt
	mkdir -p build
	cd build && cmake .. && $(MAKE)

clean:
	rm -rf build
