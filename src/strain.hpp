#pragma once
#include <vector>
#include "global_typedefs.hpp"

Antigen init_genotype_mask();

Antigen get_phenotype_id(const Antigen antigen);

Antigen get_genotype_id(const Antigen antigen); //Returns just the first NUM_GENOTYPE_ONLY_BITS bits, referring to the non-phenotype coding portion of the antigen.

std::string strain_phenotype_str(const Strain& strain);

std::string strain_phenotype_str_ordered(const Strain& strain);

//Returns a random antigen from the whole of genotypic / antigenic space.
Antigen random_antigen();

//Generates a strain from the given pool of antigens.
Strain strain_from_antigen_pool(const std::vector<Antigen>& pool);

//Intragenic recombination
Strain generate_recombinant_strain(const Strain& parent1);

//Intergenic recombination
Strain generate_recombinant_strain(const Strain& parent1, const Strain& parent2);

Antigen recombinant_antigen(const Antigen a, const Antigen b);
