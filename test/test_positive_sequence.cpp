#define BOOST_TEST_MODULE positive_sequence

#include "test_generic_sequence.hpp"

#include "positive_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "uniform_partitioned_sequence.hpp"
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <numeric>

template <typename BaseSequence>
void test_positive_sequence()
{
    srand(42);
    quasi_succinct::global_parameters params;
    size_t n = 50000;
    std::vector<uint64_t> values(n);
    std::generate(values.begin(), values.end(), []() { return (rand() % 256) + 1; });
    uint64_t universe = std::accumulate(values.begin(), values.end(), 0) + 1;

    typedef quasi_succinct::positive_sequence<BaseSequence> sequence_type;
    succinct::bit_vector_builder bvb;
    sequence_type::write(bvb, values.begin(), universe, values.size(), params);
    succinct::bit_vector bv(&bvb);
    typename sequence_type::enumerator r(bv, 0, universe, values.size(), params);

    std::cout << "Builder size: " << bvb.size() << " BPK: " << (1.0 * bvb.size() / values.size()) << std::endl;
    std::cout << "Index size: " << bv.size() << " BPK: " << (1.0 * bv.size() / values.size()) << std::endl;

    for (size_t i = 0; i < n; ++i) {
        auto val = r.move(i);
        MY_REQUIRE_EQUAL(i, val.first,
                         "i = " << i);
        MY_REQUIRE_EQUAL(values[i], val.second,
                         "i = " << i);
    }
}

BOOST_AUTO_TEST_CASE(positive_sequence)
{
    test_positive_sequence<quasi_succinct::strict_sequence>();
    test_positive_sequence<quasi_succinct::partitioned_sequence<quasi_succinct::strict_sequence>>();
    test_positive_sequence<quasi_succinct::uniform_partitioned_sequence<quasi_succinct::strict_sequence>>();
}
