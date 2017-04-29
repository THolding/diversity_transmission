#include "strain.hpp"
#include "param_manager.hpp"
#include "utilities.hpp"
#include <sstream>
#include <algorithm>


const Antigen GENOTYPE_MASK = init_genotype_mask();

Antigen init_genotype_mask()
{
    Antigen mask = 0;
    for (unsigned int i=0; i<ParamManager::num_genotype_only_bits; ++i)
    {
        mask = mask << 1;
        mask+=1;
    }
    return mask;
}

//Returns just the first NUM_GENOTYPE_ONLY_BITS bits, referring to the non-phenotype coding portion of the antigen.
Antigen get_genotype_id(const Antigen antigen)
{
    return antigen & GENOTYPE_MASK;
}

//Takes the whole antigen and returns phenotype ID only.
Antigen get_phenotype_id(const Antigen antigen)
{
    return (antigen >> ParamManager::num_genotype_only_bits) % ParamManager::num_phenotypes;
}

std::string strain_phenotype_str(const Strain& strain)
{
    std::ostringstream oss;
    for (const Antigen antigen : strain)
        oss << get_phenotype_id(antigen) << " ";
    return oss.str();
}

std::string strain_phenotype_str_ordered(const Strain& strain)
{
    //sort strain
    Strain orderedStrain = strain;
    std::sort(orderedStrain.begin(), orderedStrain.end());

    return strain_phenotype_str(orderedStrain);
}

//Returns a random antigen from the whole of genotypic / antigenic space.
Antigen random_antigen()
{
    return utilities::urandom(0, ParamManager::genotypic_space_size);
}

//Generates a strain from the given pool of antigens.
Strain strain_from_antigen_pool(const std::vector<Antigen>& pool)
{
    //Strain strain(REPERTOIRE_SIZE, 0);
    Strain strain;
    strain.reserve(ParamManager::repertoire_size);
    for (unsigned int a=0; a<ParamManager::repertoire_size; ++a)
        strain.push_back(pool[utilities::random(0, pool.size())]);
    return strain;
}

//intragenic recombination
Strain generate_recombinant_strain(const Strain& parent1)
{
    Strain recombinant = parent1;
    for (unsigned int i=0; i<ParamManager::repertoire_size; ++i)
    {
        float p = utilities::random_float01();
        if (p <= ParamManager::intragenic_recombination_p) //Intragenic (gene hybrid)
        {
            recombinant[i] = recombinant_antigen(parent1[i], parent1[utilities::random(0, parent1.size())]);
        }
    }
    return recombinant;
}

//Intergenic recombination
Strain generate_recombinant_strain(const Strain& parent1, const Strain& parent2)
{
    Strain recombinant = parent1;
    for (unsigned int i=0; i<ParamManager::repertoire_size; ++i)
    {
        float p = utilities::random_float01();
        if (p <= ParamManager::intergenic_recombination_p) //Intergenic (swap gene)
        {
            recombinant[i] = parent2[i];
        }
    }
    return recombinant;
}

Antigen recombinant_antigen(const Antigen a, const Antigen b)
{
    float antigenDifference = std::abs( (long)get_phenotype_id(a) - (long)get_phenotype_id(b) );
    long recombinant = a + (utilities::random_float_m1_1() * ParamManager::recombination_scale * antigenDifference);
    return (Antigen)(recombinant % ParamManager::genotypic_space_size);
}
