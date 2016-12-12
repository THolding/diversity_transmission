#include "demographic_tools.hpp"
#include "utilities.hpp"
#include <cmath>

PTABLE generate_host_ptable()
{
    PTABLE pDeath;
    float lambda = -0.0000015;
    for (unsigned int age=0; age<pDeath.size(); ++age)
        pDeath[age] = std::exp(-((float)age*lambda)) - 1.0;
    utilities::arrayToFile(pDeath, FILE_PATH+RUN_NAME+"_host_pdeath.csv");
    return pDeath;
}

PTABLE generate_mosquito_ptable()
{
    PTABLE pDeath;
    for (unsigned int age=0; age<pDeath.size(); ++age)
        pDeath[age] = 0.25 / (1.0+std::exp(-(0.5*((float)age - (float)MEAN_MOSQUITO_LIFE_EXPECTANCY))));
    utilities::arrayToFile(pDeath, FILE_PATH+RUN_NAME+"_mosquito_pdeath.csv");
    return pDeath;
}

//Returns cumulative density function for host pDeath tables (taking yearly format into account).
PTABLE calculate_host_cdf(const PTABLE& pDeath)
{
    PTABLE cdf;
    cdf[0] = std::pow(1.0-pDeath[0], 365);
    for (unsigned int i=1; i<pDeath.size(); ++i)
        cdf[i] = cdf[i-1] *= std::pow(((1.0 - pDeath[i])), 365);
    utilities::arrayToFile(cdf, FILE_PATH+RUN_NAME+"_host_cdf.csv");
    return cdf;
}

//Returns cumulative density function for mosquito pDeath tables (daily).
PTABLE calculate_mosquito_cdf(const PTABLE& pDeath)
{
    PTABLE cdf;
    cdf[0] = 1.0-pDeath[0];
    for (unsigned int i=1; i<pDeath.size(); ++i)
        cdf[i] = (cdf[i-1] *= (1.0 - pDeath[i]));
    utilities::arrayToFile(cdf, FILE_PATH+RUN_NAME+"_mosquito_cdf.csv");
    return cdf;
}

unsigned int random_host_equilibrum_age(const PTABLE& cdf)
{
    float survivalP = utilities::random_float01();
    unsigned int ageYears = 0;
    while (cdf[ageYears]>=survivalP)
        ++ageYears;

    short ageDays = utilities::random(0, 365);
    return (ageYears*365) + ageDays;
}

unsigned int random_mosquito_equilibrium_age(const PTABLE& cdf)
{
    float survivalP = utilities::random_float01();
    unsigned int age = 0;
    while (cdf[age] >= survivalP)
        ++age;
    return age;
}
