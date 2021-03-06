#include "utilities.hpp"
#include "param_manager.hpp"
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <stdexcept>

void utilities::initialise_random()
{
    auto seed = time(NULL);
    srand(seed);

    std::ofstream file;
    file.open(ParamManager::file_path()+ParamManager::run_name()+"_seed.txt", std::ofstream::out | std::ofstream::trunc);
    file << seed;

    file.flush();
    file.close();
}

int utilities::wrap(int number, int low, int high)
{
    if (number >= high)
        number -= high * (std::floor(number/high));
    else if (number < low)
    {
        //int toAdd = (std::floor(std::abs(number) / high)+1) * high;
        //number += toAdd;
        number+=(std::floor(std::abs(number) / high)+1) * high;
    }

    return number;
}

int utilities::random(int start, int end)
{
    return start + rand() % (end-start);
}


unsigned int utilities::urandom(unsigned int start, unsigned int end)
{
    return start + rand() % (end-start);
}

//Returns a uniform random float on the interval [0, 1]
float utilities::random_float01()
{
    return ((float)rand() / (RAND_MAX));
}

//Returns a number on the interval [-1, 1]
float utilities::random_float_m1_1()
{
    return (((float)rand() / (RAND_MAX))-0.5)*2.0;
}

long utilities::factorial(unsigned int k)
{
    if (k==0)
        return 1;

    if (k>20)
        throw std::runtime_error("utilities::factorial(k): k shouldn't ever be this big...");

    unsigned long output = k;
    while (k > 2)
    {
        output *= --k;
    }
    return (double)output;
}
