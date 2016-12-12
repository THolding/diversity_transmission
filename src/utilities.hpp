#pragma once
#include <array>
#include <fstream>
#include <string>

#include <iostream>

namespace utilities
{
    void initialise_random();

    int random(int start, int end);
    unsigned int urandom(unsigned int start, unsigned int end);

    float random_float01(); //Returns a number on the interval [0, 1]
    float random_float_m1_1(); //Returns a number on the interval [-1, 1]

    long factorial(unsigned int k);

    template <typename T>
    void arrayToFile(T arr, std::string filename)
    {
        std::ofstream file;
        file.open(filename, std::ofstream::out | std::ofstream::trunc);

        auto itr = arr.begin();
        while (itr != arr.end())
        {
            file << *itr << "\n";
            ++itr;
        }

        file.flush();
        file.close();
    }

    template <typename T>
    void matrixToFile(T matrix, std::string filename, std::string delim)
    {
        std::ofstream file;
        file.open(filename, std::ofstream::out | std::ofstream::trunc);

        for (unsigned int row=0; row<matrix[0].size(); ++row)
        {
            for (unsigned int col=0; col<matrix.size(); ++col)
            {
                if (col != 0)
                    file << delim;
                if (matrix[col].size() <= row)
                    file << 0;
                else
                    file << matrix[col][row];
            }
            file << "\n";
        }

        file.flush();
        file.close();
    }
}
