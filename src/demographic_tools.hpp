#pragma once
#include "global_parameters.hpp"
#include "host.hpp"
#include "mosquito.hpp"

PTABLE generate_host_ptable();

PTABLE generate_mosquito_ptable();

PTABLE calculate_host_cdf(const PTABLE& pDeath); //Returns cumulative density function for host pDeath tables (taking yearly format into account).

PTABLE calculate_mosquito_cdf(const PTABLE& pDeath); //Returns cumulative density function for mosquito pDeath tables (daily).

unsigned int random_host_equilibrum_age(const PTABLE& cdf);

unsigned int random_mosquito_equilibrium_age(const PTABLE& cdf);
