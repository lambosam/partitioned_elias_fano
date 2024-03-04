#define BOOST_TEST_MODULE my_test

#include <fstream>
#include <iostream>
#include <algorithm>
#include <thread>
#include <numeric>
#include <vector>
#include <cstdlib>

#include <succinct/mapper.hpp>

#include "my_test.hpp"


BOOST_AUTO_TEST_CASE(my_test)
{
    quasi_succinct::global_parameters params;
    using quasi_succinct::indexed_sequence;

    std::vector<double> avg_gaps = {2, 5, 10, 100, 1000};
    for (auto avg_gap: avg_gaps) {
        vector<double> res_ori;
        vector<double> res_ef;
        vector<double> res_pef_uni;
        vector<double> res_pef_opt;
        uint64_t n = 10000;
        uint64_t universe = uint64_t(n * avg_gap);
        std::cout << "N = " << n << ", M/N = " << avg_gap << std::endl;
        for (size_t i=0; i < 10; i++){
            auto seq = random_sequence(universe, n, i, true);
            uint64_t seq_size = seq.size() * sizeof(seq[0]);
            double seq_size_bpk = 8.0 / seq.size() * seq_size;

            res_ori.push_back(seq_size_bpk);
            res_ef.push_back(cal_test_positive_sequence<indexed_sequence>(universe, seq));
            res_pef_uni.push_back(cal_test_sequence(quasi_succinct::uniform_partitioned_sequence<indexed_sequence>(), params, universe, seq));
            res_pef_opt.push_back(cal_test_partitioned_sequence<indexed_sequence>(universe, seq));
        }
        
        std::cout << "Original size (BPK): " << (std::accumulate(res_ori.begin(), res_ori.end(), 0.0) / res_ori.size()) << std::endl;
        std::cout << "EF index size (BPK): " << (std::accumulate(res_ef.begin(), res_ef.end(), 0.0) / res_ef.size()) << std::endl;
        std::cout << "Uniform PEF index size (BPK): " << (std::accumulate(res_pef_uni.begin(), res_pef_uni.end(), 0.0) / res_pef_uni.size()) << std::endl;
        std::cout << "Opt PEF index size (BPK): " << (std::accumulate(res_pef_opt.begin(), res_pef_opt.end(), 0.0) / res_pef_opt.size()) << "\n" << std::endl;            

    }
}
