#pragma once
#include <vector>
#include "global_typedefs.hpp"

//Keeps track of the number of antigens in circulation.
class DiversityMonitor
{
private:
    std::vector<unsigned int> antigenCounts;

    unsigned int totalAntigens;
    unsigned int uniqueAntigens;
    unsigned int numNewlyGenerated;
    unsigned int numExtinctions;

    DiversityMonitor(); //Singleton.
public:
    static DiversityMonitor& instance() //Singleton instance.
    {
        static DiversityMonitor diversityMonitor; //Implicit if flag ensures initialisation occurs only once!
        return diversityMonitor;
    }

    static void reset();

    static void register_antigen_gain(Antigen phenotypeID);

    static void register_antigen_loss(Antigen phenotypeID);

    static void register_new_strain(Strain geneList);

    static void register_lost_strain(Strain geneList);

    static void reset_loss_gen_count();

    static unsigned int get_antigen_count(const unsigned int phenotypeID);
    static const std::vector<unsigned int>& get_antigen_counts();
    static unsigned int get_total_antigens();
    static unsigned int get_num_unique_antigens();
    static unsigned int get_current_generation_count();
    static unsigned int get_current_loss_count();
};
