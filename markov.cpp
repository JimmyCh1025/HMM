#include "markov.h"

static void InitialDecodeDNA()
{
    decodeDNA['A'] = 0;
    decodeDNA['T'] = 1;
    decodeDNA['C'] = 2;
    decodeDNA['G'] = 3;

    return;
}

static void InitializeArr(int kmers)
{
    memset(seqTraningCnt, 0, sizeof(seqTraningCnt));
    memset(seqTest1Cnt, 0, sizeof(seqTest1Cnt));
    memset(seqTest2Cnt, 0, sizeof(seqTest2Cnt));
    memset(seqTraningProb, 0, sizeof(seqTraningProb));

    return;
}

static void ComputeTrainingCntAndProb(int kmers, std::string sTraining, std::string sDTraining, double *prevProb, double *prevHeadAndTailProb)
{
    int charCnt = 0;
    int hashTableKey = 0;
    int maxKey = pow(4,kmers)-1;
    char ch = '\0';

    // DNA sequence 
    for (int i = 0 ; i < seqSize ; ++i)
    {
        ch = sTraining[i];
        if (ch == 'N') {
            charCnt = 0;
            hashTableKey = 0;
            continue;
        }

        hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
        if (++charCnt >= kmers)
        {
            ++seqTraningCnt[hashTableKey];
        }

    }


    charCnt = 0;
    hashTableKey = 0;
    
    // DNA double stranded 
    for (int i = 0 ; i < seqSize ; ++i)
    {
        ch = sDTraining[i];
        if (ch == 'N') {
            charCnt = 0;
            hashTableKey = 0;
            continue;
        }

        hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
        if (++charCnt >= kmers)
        {
            ++seqTraningCnt[hashTableKey];
        }
    }

    ComputeTrainingProb(kmers, sTraining, prevProb, prevHeadAndTailProb);

    return;
}

static void ComputeTrainingProb(int kmers, std::string sTraining, double *prevProb, double *prevHeadAndTailProb)
{
    int maxCombination = pow(4, kmers);
    int totalCount = 0;
    int count = 0;
    double trainingProbMarkov = 0.0;

    probMarkov[kmers-1][0] = -(*prevProb)+(*prevHeadAndTailProb);

    for (int i = 0 ; i < maxCombination ; ++i)
        totalCount += seqTraningCnt[i];

    for (int i = 0 ; i < maxCombination ; ++i)
    {
        seqTraningProb[i] = ((double)seqTraningCnt[i] + 0.5)/(0.5*maxCombination + totalCount);
        // trainingProbMarkov += seqTraningCnt[i]*log2(seqTraningProb[i]);
    }


    int charCnt = 0;
    int hashTableKey = 0;
    int kMerCount = 0;
    double tailProb ;
    int maxKey = pow(4,kmers)-1;
    char ch = '\0';

    // // DNA sequence 
    for (int i = 0 ; i < seqSize ; ++i)
    {
        ch = sTraining[i];
        if (ch == 'N') {
            charCnt = 0;
            hashTableKey = 0;
            continue;
        }

        hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
        if (++charCnt >= kmers)
        {
            trainingProbMarkov += log2(seqTraningProb[hashTableKey]);
            kMerCount++;
            
            if (kMerCount == 1)
            {
                *prevHeadAndTailProb = log2(seqTraningProb[hashTableKey]);
            }
            
            tailProb = log2(seqTraningProb[hashTableKey]);
        }

    }

    // printf("%d\n", kmers);
    // printf("%d\n", totalCount);
    // printf("%f\n",trainingProbMarkov);
    probMarkov[kmers-1][0] += trainingProbMarkov;
    *prevHeadAndTailProb += tailProb;
    *prevProb = trainingProbMarkov;
    return;
}

static void ComputeTest1CntAndProb(int kmers, std::string sTest1, std::string sDTest1, double *prevProb, double *prevHeadAndTailProb)
{
    int charCnt = 0;
    int hashTableKey = 0;
    int maxKey = pow(4,kmers)-1;
    int kMerCount = 0;
    double tailProb ;
    double test1ProbMarkov = 0.0;
    char ch = '\0';

    probMarkov[kmers-1][1] = -(*prevProb)+(*prevHeadAndTailProb);

    // DNA sequence 
    for (int i = 0 ; i < seqSize ; ++i)
    {
        ch = sTest1[i];
        if (ch == 'N') {
            charCnt = 0;
            hashTableKey = 0;
            continue;
        }

        hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
        if (++charCnt >= kmers)
        {
            test1ProbMarkov += log2(seqTraningProb[hashTableKey]);
            kMerCount++;
            
            if (kMerCount == 1)
            {
                *prevHeadAndTailProb = log2(seqTraningProb[hashTableKey]);
            }
            
            tailProb = log2(seqTraningProb[hashTableKey]);
            ++seqTest1Cnt[hashTableKey];
        }

    }
    
    probMarkov[kmers-1][1] += test1ProbMarkov;
    *prevHeadAndTailProb += tailProb;
    *prevProb = test1ProbMarkov;


    // charCnt = 0;
    // hashTableKey = 0;
    
    // DNA double stranded 
    // for (int i = 0 ; i < seqSize ; ++i)
    // {
    //     ch = sDTest1[i];
    //     if (ch == 'N') {
    //         charCnt = 0;
    //         hashTableKey = 0;
    //         continue;
    //     }

    //     hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
    //     if (++charCnt >= kmers)
    //     {
    //         ++seqTest1Cnt[hashTableKey];
    //     }

    // }

    // ComputeTest1Prob(kmers, prevProb, prevHeadAndTailProb);
    return;
}

static void ComputeTest1Prob(int kmers, double *prevProb, double *prevHeadAndTailProb)
{
    int maxCombination = pow(4, kmers);
    int count = 0;
    double test1ProbMarkov = 0.0;

    for (int i = 0 ; i < maxCombination ; ++i)
    {
        test1ProbMarkov += (log2(seqTraningProb[i])*seqTest1Cnt[i]);
    }

    probMarkov[kmers-1][1] = test1ProbMarkov;
    
    // probMarkov[kmers-1][1] -= (*prevProb + *prevHeadAndTailProb);
    return;
}

static void ComputeTest2CntAndProb(int kmers, std::string sTest2, std::string sDTest2, double *prevProb, double *prevHeadAndTailProb)
{
    int charCnt = 0;
    int hashTableKey = 0;
    int maxKey = pow(4,kmers)-1;
    int kMerCount = 0;
    double tailProb ;
    double test2ProbMarkov = 0.0;
    char ch = '\0';

    probMarkov[kmers-1][2] = -(*prevProb)+(*prevHeadAndTailProb);

    // DNA sequence 
    for (int i = 0 ; i < seqSize ; ++i)
    {
        ch = sTest2[i];
        if (ch == 'N') {
            charCnt = 0;
            hashTableKey = 0;
            continue;
        }

        hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
        if (++charCnt >= kmers)
        {
            test2ProbMarkov += log2(seqTraningProb[hashTableKey]);
            kMerCount++;
            
            if (kMerCount == 1)
            {
                *prevHeadAndTailProb = log2(seqTraningProb[hashTableKey]);
            }
            
            tailProb = log2(seqTraningProb[hashTableKey]);
            ++seqTest2Cnt[hashTableKey];
        }

    }
    
    probMarkov[kmers-1][2] += test2ProbMarkov;
    *prevHeadAndTailProb += tailProb;
    *prevProb = test2ProbMarkov;

    
    // charCnt = 0;
    // hashTableKey = 0;
    
    // DNA double stranded 
    // for (int i = 0 ; i < seqSize ; ++i)
    // {
    //     ch = sDTest2[i];
    //     if (ch == 'N') {
    //         charCnt = 0;
    //         hashTableKey = 0;
    //         continue;
    //     }

    //     hashTableKey = ((hashTableKey << 2) | decodeDNA[ch]) & maxKey;
    //     if (++charCnt >= kmers)
    //     {
    //         ++seqTest2Cnt[hashTableKey];
    //     }

    // }

    // ComputeTest2Prob(kmers, prevProb, prevHeadAndTailProb);

    return;
}

static void ComputeTest2Prob(int kmers, double *prevProb, double *prevHeadAndTailProb)
{
    int maxCombination = pow(4, kmers);
    int count = 0;
    double test2ProbMarkov = 0.0;

    for (int i = 0 ; i < maxCombination ; ++i)
    {
        test2ProbMarkov += (log2(seqTraningProb[i])*seqTest2Cnt[i]);
    }

    probMarkov[kmers-1][2] = test2ProbMarkov;

    return;
}

static void PrintAllMarkovProb()
{   
    printf("Order           Training           Test1           Test2\n");
    for (int i = 0 ; i < maxMer ; ++i)
    {
        std::cout << std::right << std::setw(5) << i << "   ";
        std::cout << std::right << std::setw(16) << probMarkov[i][0] ;
        std::cout << std::right << std::setw(16) << probMarkov[i][1] ;
        std::cout << std::right << std::setw(16) << probMarkov[i][2] << std::endl;
    }

    return;
}

void ComputeMarkov(GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq)
{
    double prevProb = 0, prevKmersHeadAndTailProb = 0;
    double prevProb1 = 0, prevKmersHeadAndTailProb1 = 0;
    double prevProb2 = 0, prevKmersHeadAndTailProb2 = 0;

    InitialDecodeDNA(); 
    
    for (int i = 1 ; i <= maxMer ; ++i)
    {
        InitializeArr(i);
        ComputeTrainingCntAndProb(i, genomeSeq.straining, genomeDSeq.straining, &prevProb, &prevKmersHeadAndTailProb);
        ComputeTest1CntAndProb(i, genomeSeq.stest1, genomeDSeq.stest1, &prevProb1, &prevKmersHeadAndTailProb1);
        ComputeTest2CntAndProb(i, genomeSeq.stest2, genomeDSeq.stest2, &prevProb2, &prevKmersHeadAndTailProb2);
    }

    PrintAllMarkovProb();

    return;
}