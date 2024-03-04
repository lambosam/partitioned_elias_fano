#define BOOST_TEST_MODULE uniform_partitioned_sequence

#include "test_generic_sequence.hpp"

#include "uniform_partitioned_sequence.hpp"
#include "strict_sequence.hpp"
#include <vector>
#include <cstdlib>

BOOST_AUTO_TEST_CASE(uniform_partitioned_sequence)
{
    quasi_succinct::global_parameters params;
    using quasi_succinct::indexed_sequence;
    using quasi_succinct::strict_sequence;

    // test singleton sequences
    std::vector<uint64_t> short_seq;
    short_seq.push_back(0);
    test_sequence(quasi_succinct::uniform_partitioned_sequence<indexed_sequence>(),
                  params, 1, short_seq);
    test_sequence(quasi_succinct::uniform_partitioned_sequence<strict_sequence>(),
                  params, 1, short_seq);
    short_seq[0] = 1;
    test_sequence(quasi_succinct::uniform_partitioned_sequence<indexed_sequence>(),
                  params, 2, short_seq);
    test_sequence(quasi_succinct::uniform_partitioned_sequence<strict_sequence>(),
                  params, 2, short_seq);

    std::cout << "Uniform EF Index" << std::endl;

    std::vector<double> avg_gaps = { 1.1, 1.9, 2.5, 3, 4, 5, 10 };
    for (auto avg_gap: avg_gaps) {
        uint64_t n = 10000;
        uint64_t universe = uint64_t(n * avg_gap);
        auto seq = random_sequence(universe, n, true);

        std::cout << "Gap size:" << avg_gap << std::endl;

        uint64_t seq_size = seq.size() * sizeof(seq[0]);
        double seq_size_bpk = 8.0 / seq.size() * seq_size;
        std::cout << "Original size: " << seq_size << " BPLK: " << seq_size_bpk << std::endl;

        test_sequence(quasi_succinct::uniform_partitioned_sequence<indexed_sequence>(),
                      params, universe, seq);
        test_sequence(quasi_succinct::uniform_partitioned_sequence<strict_sequence>(),
                      params, universe, seq);
    }
}
