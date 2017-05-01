#pragma once

#include <vector>
#include "global_typedefs.hpp"
#include "host.hpp"
#include "mosquito.hpp"

namespace testing
{
    void test_diversity_counting();
    void long_diversity_count(unsigned int& uniqueCount, unsigned int& totalCount, const std::vector<Host>& hosts, const std::vector<Mosquito>& mosquitoes);
    void long_diversity_count_helper(std::vector<unsigned int>& antigenFreqs, unsigned int& antigenCounter, const Strain& strain, unsigned int &totalCount);

    void new_tests();

    void run_tests();

    void test_immunity();
    void test_host_infection();
}
