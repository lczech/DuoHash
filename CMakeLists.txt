cmake_minimum_required(VERSION 3.10)
project(DuoHash)

# Imposta il compilatore g++
# Nope, hard coding a compiler is not a good idea.
# Let's leave that for the caller to decide.
# set(CMAKE_CXX_COMPILER g++-13)

# Specifica la versione di C++ da utilizzare
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Aggiungi le opzioni di ottimizzazione
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Aggiungi i file header e source
include_directories(${PROJECT_SOURCE_DIR}/include)

# Definisci le directory dei sorgenti
file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

# Crea la libreria statica
add_library(DuoHash STATIC ${SOURCES})

# Aggiungi la libreria agli include
target_include_directories(DuoHash PUBLIC ${PROJECT_SOURCE_DIR}/include)

# ==================================================================================================
#     Custom additions to simplify include path finding
# ==================================================================================================

# Build example executable
add_executable(benchmark ${PROJECT_SOURCE_DIR}/example/benchmark.cpp)

# Link against the library built above
target_link_libraries(benchmark PRIVATE DuoHash)

# Ensure it also sees public headers (already provided by DuoHash PUBLIC includes)
target_include_directories(benchmark PRIVATE ${PROJECT_SOURCE_DIR}/include)

# Place the executable into the source example/ directory
set_target_properties(benchmark PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/example
)
