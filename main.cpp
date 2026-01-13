#include "main.h"


int main()
{
    ifstream myfile;

    // read file 
    myfile.open("GCA_000001405.15_GRCh38_no_alt_analysis_set.fna");
    if (myfile.is_open())
    {
        // Finding sequence
        printf("======================= Find Seq =======================\n");
        FindGenomeSequence(myfile, genomeSeq, genomeDStrandedSeq);
        myfile.close();


        // Computing markov
        printf("\n======================== Part 1 ========================\n");
        ComputeMarkov(genomeSeq, genomeDStrandedSeq);

        // Training HMM
        printf("\n======================== Part 2 ========================\n");
        thread emHMM3(TrainEMHMM, state3, ref(genomeSeq), ref(genomeDStrandedSeq));
        thread emHMM5(TrainEMHMM, state5, ref(genomeSeq), ref(genomeDStrandedSeq));
        
        emHMM3.join();
        emHMM5.join();

        PrintAllSeqLogLikelihood();
        printf("==========================================================\n");
    }
    else
    {
        printf("Failed to open file.\n");
    }

    return 0;
}