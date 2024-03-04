#pragma once

#include <vector>
#include <cstdlib>
#include <iostream>

#include "succinct/test_common.hpp"
#include "succinct/bit_vector.hpp"
#include "util.hpp"

#include "../index_types.hpp"
#include "uniform_partitioned_sequence.hpp"
#include "partitioned_sequence.hpp"
#include "compact_elias_fano.hpp"

#include "test_generic_sequence.hpp"
#include "strict_sequence.hpp"

#include "positive_sequence.hpp"


namespace quasi_succinct {

    class partitioned_sequence_test {
    public:
        template <typename Enumerator>
        static void test_construction(Enumerator& r, std::vector<uint64_t> const& seq)
        {
            if (r.m_partitions == 1) { // nothing to test here
                return;
            }

            for (size_t p = 0; p < r.m_partitions; ++p) {
                r.switch_partition(p);

                uint64_t cur_begin = r.m_cur_begin;
                uint64_t cur_end = r.m_cur_end;

                uint64_t cur_base = p ? seq[cur_begin - 1] + 1 : seq[0];
                uint64_t cur_upper_bound = seq[cur_end - 1];
                MY_REQUIRE_EQUAL(cur_base, r.m_cur_base,
                                 "p = " << p);
                MY_REQUIRE_EQUAL(cur_upper_bound, r.m_cur_upper_bound,
                                 "p = " << p);

                for (uint64_t i = cur_begin; i < cur_end; ++i) {
                    auto val = r.m_partition_enum.move(i - cur_begin);
                    MY_REQUIRE_EQUAL(seq[i], cur_base + val.second,
                                     "p = " << p << " i = " << i);
                }
            }
        }
    };
}

template <typename BaseSequence>
double cal_test_partitioned_sequence(uint64_t universe,
                               std::vector<uint64_t> const& seq)
{
    quasi_succinct::global_parameters params;
    typedef quasi_succinct::partitioned_sequence<BaseSequence> sequence_type;

    succinct::bit_vector_builder bvb;
    sequence_type::write(bvb, seq.begin(), universe, seq.size(), params);
    succinct::bit_vector bv(&bvb);

    typename sequence_type::enumerator r(bv, 0, universe, seq.size(), params);
    quasi_succinct::partitioned_sequence_test::test_construction(r, seq);
    test_sequence(r, seq);

    return (1.0 * bv.size() / seq.size());
}

template <typename BaseSequence>
double cal_test_positive_sequence(uint64_t universe,
                            std::vector<uint64_t> const& seq)
{
    universe = std::accumulate(seq.begin(), seq.end(), 0) + 1;

    quasi_succinct::global_parameters params;

    typedef quasi_succinct::positive_sequence<BaseSequence> sequence_type;
    succinct::bit_vector_builder bvb;
    sequence_type::write(bvb, seq.begin(), universe, seq.size(), params);
    succinct::bit_vector bv(&bvb);
    typename sequence_type::enumerator r(bv, 0, universe, seq.size(), params);

    return (1.0 * bv.size() / seq.size());
}
