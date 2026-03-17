// Compiled via cmake for simplicity of finding link and include directories

#include <DuoHash.h>
#include <chrono>
#include <iomanip>
#include <string>
#include <cctype>
#include <random>
#include <array>

// ---------------------------------------------------------------------------------------------
//     Benchmark Helpers
// ---------------------------------------------------------------------------------------------

// Generate a random sequence of ACGT of a given length
inline std::string random_sequence(std::size_t len)
{
    static thread_local std::mt19937_64 rng{std::random_device{}()};
    static constexpr std::array<char, 4> bases{'A','C','G','T'};
    std::uniform_int_distribution<size_t> dist(0, 3);

    std::string sequence;
    sequence.resize(len);
    for( std::size_t i = 0; i < len; ++i ) {
        sequence[i] = bases[dist(rng)];
    }
    return sequence;
}

// Generate a set of random sequences of ACGT
inline std::vector<std::string> random_sequences( std::size_t num, std::size_t len )
{
    std::vector<std::string> result;
    result.resize( num );
    for( std::size_t i = 0; i < num; ++i ) {
        result[i] = random_sequence( len );
    }
    return result;
}

// Get the total number of bases in a set of sequences
inline size_t total_size( std::vector<std::string> const& sequences )
{
    size_t s = 0;
    for( auto const& seq : sequences ) {
        s += seq.length();
    }
    return s;
}

// Run a timed benchmark of a DuoHash method
template<class Container, class DuoHashC>
void DuoHash_timed_call(
    std::string const& label,
    DuoHashC& obj,
    void (DuoHashC::*method)(std::vector<Container>&)
) {
    // Run timed benchmark
    std::vector<Container> encodings;
    auto t1 = std::chrono::steady_clock::now();
    (obj.*method)(encodings);
    auto t2 = std::chrono::steady_clock::now();

    // Get time and number of bases processed
    auto const t = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1);
    auto const s = static_cast<double>(
        obj.getSpacedQmerCount() * total_size( obj.getSequences() )
    );

    // Print result
    std::cout << label << "," << ( static_cast<double>(t.count()) / s ) << "\n";
}

// ---------------------------------------------------------------------------------------------
//     Benchmark
// ---------------------------------------------------------------------------------------------

void benchmark(
    std::vector<std::string> const& sequences,
    std::vector<SpacedQmer> const& multi_spaced,
    std::string const& bench_case
) {
    // ------------------------------------------------------------
    //     Test with single seed
    // ------------------------------------------------------------

    // std::cout << "\nTest with single seed\n";
    DuoHash test_single(sequences, multi_spaced[0]);
    auto const sname = "single," + bench_case;

    // Original functions
    DuoHash_timed_call<Encoding_V>(
        sname + ",naive", test_single, &DuoHash::GetEncoding_naive
    );
    DuoHash_timed_call<Encoding_V>(
        sname + ",FSH", test_single, &DuoHash::GetEncoding_FSH
    );
    DuoHash_timed_call<Encoding_V>(
        sname + ",ISSH", test_single, &DuoHash::GetEncoding_ISSH
    );

    // ------------------------------------------------------------
    //     Test with multiple seed
    // ------------------------------------------------------------

    // std::cout << "\nTest with multiple seed\n";
    DuoHash_multi test_multi(sequences, multi_spaced);
    auto const mname = "multi," + bench_case;

    // Original functions
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",naive", test_multi, &DuoHash_multi::GetEncoding_naive
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",FSH", test_multi, &DuoHash_multi::GetEncoding_FSH
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",ISSH", test_multi, &DuoHash_multi::GetEncoding_ISSH
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",FSH_multi", test_multi, &DuoHash_multi::GetEncoding_FSH_multi
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",MISSH_v1", test_multi, &DuoHash_multi::GetEncoding_MISSH_v1
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",MISSH_col", test_multi, &DuoHash_multi::GetEncoding_MISSH_col
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",MISSH_col_parallel", test_multi, &DuoHash_multi::GetEncoding_MISSH_col_parallel
    );
    DuoHash_timed_call<Encoding_V_V>(
        mname + ",MISSH_row", test_multi, &DuoHash_multi::GetEncoding_MISSH_row
    );
}

// ---------------------------------------------------------------------------------------------
//     Main
// ---------------------------------------------------------------------------------------------

int main()
{
    // Generate input data. DuoHash makes a copy of this.
    // Bit wasteful, but ok for a benchmark of this size.
    size_t num_seq = 500000;
    size_t len_seq = 150;
    auto const sequences = random_sequences( num_seq, len_seq );

    // Use all seed sets offered by DuoHash, except the largest one
    // whose inital size is beyond 64bit words. Future work.
    std::string const seed_dir = "./Seeds/";
    std::vector<std::string> const seed_files = {{
        "W10L15", "W14L31", "W18L31", "W22L31", "W26L31", "W32L45",
    }};

    // Run the benchmarks on all seeds, printing a table to std out.
    std::cout << "suite,case,benchmark,ns_per_op\n";
    for( auto const& seed_file : seed_files ) {
        // std::cout << "\n=====================================\n";
        // std::cout << "Running benchmark on " << seed_file << " seeds\n";

        std::vector<SpacedQmer> multi_spaced;
        std::vector<std::string> tmp;
        if (!loadFile(seed_dir + seed_file + ".fna", tmp)) {
            std::exit(1);
        }
        for( auto seed : tmp ) {
            multi_spaced.push_back(SpacedQmer(seed, 0));
        }

        benchmark( sequences, multi_spaced, seed_file );
    }
}
