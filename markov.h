#ifndef MARKOV_H
#define MARKOV_H

#include <string>
#include <iomanip>
#include "genome.h"

#define dummy        10
#define seqArrSize   1048576 + dummy  // arr size = 4^10 + 10
#define seqSize      5000000


#define maxMer       10
#define totalKMers   maxMer + 1       // 1(dummy)
#define totalSeq     3                // training, test1, test2 

static int seqTraningCnt[seqArrSize];
static int seqTest1Cnt[seqArrSize];
static int seqTest2Cnt[seqArrSize];

static double seqTraningProb[seqArrSize];
static double probMarkov[totalKMers][totalSeq];

const char dnaBase[4] = {'A', 'T', 'C', 'G'};
static int decodeDNA[130] = {0};

static void InitializeArr(int kmers);
static void ComputeTrainingCntAndProb(int kmers, std::string sTraining, std::string sDTraining, double *prevProb, double *prevHeadAndTailProb);
static void ComputeTrainingProb(int kmers, std::string sTraining, double *prevProb, double *prevHeadAndTailProb);

static void ComputeTest1CntAndProb(int kmers, std::string sTest1, std::string sDTest1, double *prevProb, double *prevHeadAndTailProb);
static void ComputeTest1Prob(int kmers, double *prevProb, double *prevHeadAndTailProb);

static void ComputeTest2CntAndProb(int kmers, std::string sTest2, std::string sDTest2, double *prevProb, double *prevHeadAndTailProb);
static void ComputeTest2Prob(int kmers, double *prevProb, double *prevHeadAndTailProb);

static void PrintAllMarkovProb();

void ComputeMarkov(GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq);


#endif