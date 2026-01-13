#ifndef HMM_H
#define HMM_H

#include <thread>
#include "genome.h"

#define seqSize      5000000
#define threshold    0.0001   // error ≤ threshold



static double init3SProb[3] = {0.35, 0.35, 0.3};
static double trans3SProb[3][3] = {
    {0.7, 0.15, 0.15}, 
    {0.15, 0.7, 0.15},
    {0.15, 0.15, 0.7}   
};

static double emission3bProb[3][4] = {
    {0.25, 0.25, 0.25, 0.25}, 
    {0.15, 0.25, 0.25, 0.35}, 
    {0.25, 0.25, 0.25, 0.25}
};

static double init5SProb[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
static double trans5SProb[5][5] = {
    {0.6, 0.1, 0.1, 0.1, 0.1}, 
    {0.1, 0.6, 0.1, 0.1, 0.1},
    {0.1, 0.1, 0.6, 0.1, 0.1}, 
    {0.1, 0.1, 0.1, 0.6, 0.1},
    {0.1, 0.1, 0.1, 0.1, 0.6}
};

static double emission5bProb[5][4] = {
    { 0.4,  0.3,  0.2,  0.1}, 
    { 0.3,  0.4,  0.2,  0.1}, 
    { 0.2,  0.3,  0.4,  0.1}, 
    {0.25, 0.25, 0.25, 0.25}, 
    { 0.3,  0.2,  0.4,  0.1}
};

static double alpha3State[seqSize][3];
static double alpha3StateRev[seqSize][3];
static double alpha5State[seqSize][5]; 
static double alpha5StateRev[seqSize][5];

static double beta3State[seqSize][3];
static double beta3StateRev[seqSize][3];
static double beta5State[seqSize][5]; 
static double beta5StateRev[seqSize][5];

static double model3Prob[3] ;
static double model5Prob[3] ;

static int decodeDNAHMM[130] = {0};

static double forwardProb;

static bool bCanPrint3State = true;

static void InitialDecodeDNA();
static double LogSumExp(double log_a, double log_b);
static void Forward(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining, double *alphaS, double *alphaSD);
static void Backward(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining,  double *betaS, double *betaSD);
static void ComputeGammaXiAndUpdateModel(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining, double *forwardP, double *forwardN, double *backwardP, double *backwardN, bool *bConvergen);
static void ComputeFinishModelLikelihood(int state, double *modelProb, double *init, double *trans, double *emission, GenomeSequence &genomeSeq);
static void PrintUpdateInfoMatrix(int state, double *init, double *trans, double *emission);
void PrintAllSeqLogLikelihood();
void TrainEMHMM(int state, GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq);

#endif