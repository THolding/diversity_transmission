#include "diversity_monitor.hpp"
#include "param_manager.hpp"
#include "strain.hpp"
#include <iostream>

DiversityMonitor::DiversityMonitor()
{
    //DiversityMonitor::reset(); Commented because constructor cannot rely on ParamManager. ParamManager calls DiversityManager::reset.
}

void DiversityMonitor::reset()
{
    instance().totalAntigens = 0;
    instance().uniqueAntigens = 0;
    instance().numExtinctions = 0;
    instance().numNewlyGenerated = 0;
    instance().antigenCounts = std::vector<unsigned int> (ParamManager::num_phenotypes , 0);
}

void DiversityMonitor::register_antigen_gain(Antigen phenotypeID, bool bypassGenerationRegister)
{
    if (instance().antigenCounts[phenotypeID] == 0) {
        #pragma omp atomic
        instance().uniqueAntigens++;
        if (bypassGenerationRegister == false) {
            #pragma omp atomic
            instance().numNewlyGenerated++;
        }
    }
    #pragma omp atomic
    instance().totalAntigens++;

    #pragma omp atomic
    instance().antigenCounts[phenotypeID] += 1;
}

void DiversityMonitor::register_antigen_loss(Antigen phenotypeID)
{
    #pragma omp atomic
    instance().totalAntigens--;

    #pragma omp atomic
    instance().antigenCounts[phenotypeID] -= 1;

    if (instance().antigenCounts[phenotypeID] == 0) {
        #pragma omp atomic
        instance().uniqueAntigens--;

        #pragma omp atomic
        instance().numExtinctions++;
    }
}

void DiversityMonitor::register_new_strain(Strain geneList, bool bypassGenerationRegister)
{
    for (unsigned int a=0; a<geneList.size(); ++a)
        DiversityMonitor::register_antigen_gain(get_phenotype_id(geneList[a]), bypassGenerationRegister);
}

void DiversityMonitor::register_lost_strain(Strain geneList)
{
    for (unsigned int a=0; a<geneList.size(); ++a)
        DiversityMonitor::register_antigen_loss(get_phenotype_id(geneList[a]));
}

void DiversityMonitor::reset_loss_gen_count()
{
    instance().numExtinctions = 0;
    instance().numNewlyGenerated = 0;
}

unsigned int DiversityMonitor::get_antigen_count(const unsigned int phenotypeID)
{
    return instance().antigenCounts[phenotypeID];
}

const std::vector<unsigned int>& DiversityMonitor::get_antigen_counts()
{
    return instance().antigenCounts;
}

unsigned int DiversityMonitor::get_total_antigens()
{
    return instance().totalAntigens;
}

unsigned int DiversityMonitor::get_num_unique_antigens()
{
    return instance().uniqueAntigens;
}


unsigned int DiversityMonitor::get_current_generation_count()
{
    return instance().numNewlyGenerated;
}

unsigned int DiversityMonitor::get_current_loss_count()
{
    return instance().numExtinctions;
}
