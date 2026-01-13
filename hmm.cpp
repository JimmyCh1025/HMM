#include "hmm.h"

static void InitialDecodeDNA()
{
    decodeDNAHMM['A'] = 0;
    decodeDNAHMM['T'] = 1;
    decodeDNAHMM['C'] = 2;
    decodeDNAHMM['G'] = 3;

    return;
}

static double LogSumExp(double log_a, double log_b) 
{
    if (log_a == -INFINITY) return log_b;
    if (log_b == -INFINITY) return log_a;
    double max_val = std::max(log_a, log_b);
    return max_val + log2(pow(2.0, log_a - max_val) + pow(2.0, log_b - max_val));
}

static void Forward(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining, double *alphaS, double *alphaSD)
{
    double alphaTemp = 0;
    double temp ;
    int genomePos ;
    double probTotal = 0.0, probPositive = 0.0, probNegative = 0.0;
        
    // Positive strand
    // step 1 : a1 = p1 * b1(O1)
    for (int i = 0 ; i < state ; ++i)
    {
        alphaS[i] = log2(init[i]) + log2(emission[4*i+decodeDNAHMM[sTraining[0]]]);
    }

    // step2 : ai = sum(ai-1(i)* aij * bj(Oi))
    for (int i = 1 ; i < seqSize ; ++i)
    {
        genomePos = decodeDNAHMM[sTraining[i]];
        for (int j = 0 ; j < state ; ++j)
        {
            alphaS[state*i+j] = -INFINITY;
        
            for (int k = 0 ; k < state ; ++k)
            {
                // Sk->Sj
                temp = alphaS[state*(i-1)+k] + log2(trans[state*k+j]) + log2(emission[4*j+genomePos]);
                alphaS[state*i+j] = LogSumExp(alphaS[state*i+j], temp);
            }
        }
    }


    // Negative strand
    // step 1 : a1 = p1 * b1(O1)
    for (int i = 0 ; i < state ; ++i)
    {
        alphaSD[i] = log2(init[i]) + log2(emission[4*i+decodeDNAHMM[sDTraining[0]]]);
    }

    // step2 : ai = sum(ai-1(i)* aij * bj(Oi))
    for (int i = 1 ; i < seqSize ; ++i)
    {
        genomePos = decodeDNAHMM[sDTraining[i]];
        for (int j = 0 ; j < state ; ++j)
        {
            alphaSD[state*i+j] = -INFINITY;
        
            for (int k = 0 ; k < state ; ++k)
            {
                // Sk->Sj
                temp = alphaSD[state*(i-1)+k] + log2(trans[state*k+j]) + log2(emission[4*j+genomePos]);
                alphaSD[state*i+j] = LogSumExp(alphaSD[state*i+j], temp);
            }
        }
    }

    return;
}

static void Backward(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining,  double *betaS, double *betaSD) 
{
    double temp ;
    int genomePos ;

    // Positive strand 
    // step1 : Bi(T) = 1  , for all i
    for (int i = 0 ; i < state ; ++i) 
    {
        betaS[state*(seqSize-1)+i] = 0;
    }


    // step2 : Bi(t) = sum from j = 0 to N-2( aij x bj(Ot+1) x Bj(t+1))
    for (int i = seqSize-2 ; i >= 0 ; --i)
    {
        genomePos = decodeDNAHMM[sTraining[i+1]];
        for (int j = 0 ; j < state ; ++j)
        {
            betaS[state*i+j] = -INFINITY;

            for (int k = 0 ; k < state ; ++k)
            {
                temp = log2(trans[state*j+k]) + log2(emission[4*k+genomePos]) + betaS[state*(i+1)+k];
                betaS[state*i+j] = LogSumExp(betaS[state*i+j], temp);
            }

        }

    }


    // Negative strand 
    // step1 : Bi(T) = 1  , for all i
    for (int i = 0 ; i < state ; ++i) 
    {
        betaSD[state*(seqSize-1)+i] = 0;
    }


    // step2 : Bi(t) = sum from j = 0 to N-2( aij x bj(Ot+1) x Bj(t+1))
    for (int i = seqSize-2 ; i >= 0 ; --i)
    {
        genomePos = decodeDNAHMM[sDTraining[i+1]];
        for (int j = 0 ; j < state ; ++j)
        {
            betaSD[state*i+j] = -INFINITY;

            for (int k = 0 ; k < state ; ++k)
            {
                temp = log2(trans[state*j+k]) + log2(emission[4*k+genomePos]) + betaSD[state*(i+1)+k];
                betaSD[state*i+j] = LogSumExp(betaSD[state*i+j], temp);
            }

        }

    }

    return ;
}

static void ComputeGammaXiAndUpdateModel(int state, double *init, double *trans, double *emission, std::string &sTraining, std::string &sDTraining, double *forwardP, double *forwardN, double *backwardP, double *backwardN, bool *bConvergen)
{
    double *initState = new double[state];
    
    // summation gamma only need N(state)
    double *sumGamma = new double[state];
    double *tempGamma = new double[state];

    // summation xi need N^2(state*state)
    double **sumXi = new double*[state];
    double **tempXi = new double*[state];

    // summation all state gamma
    double *sumStateGamma = new double[state];

    // To compute error(L2 norm)
    double **newTransProb = new double*[state];

    
    double **newEmissionProb = new double*[state];


    int genomePos ;
    int genomeObs ;

    double totalGammaProb = 0.0;
    double totalXiProb = 0.0;

    double rowSum = 0.0;
    double diffSquare;

    // allocate 2D dynamic array
    for (int i = 0 ; i < state ; ++i)
    {
        sumXi[i] = new double[state];
        tempXi[i] = new double[state];
        newTransProb[i] = new double[state];
        newEmissionProb[i] = new double[4];

        sumStateGamma[i] = -INFINITY;

        for (int j = 0; j < 4; ++j)
            newEmissionProb[i][j] = -INFINITY;
        
    }


    totalGammaProb = -INFINITY;
    totalXiProb = -INFINITY;

    // compute gamma and xi and update init state
    genomePos = decodeDNAHMM[sTraining[1]];
    genomeObs = decodeDNAHMM[sTraining[0]];

    for (int i = 0 ; i < state ; ++i)
    {
        // forward positive and backward positive (log2)
        // assign alpha*beta to gamma
        sumGamma[i] = (forwardP[i] + backwardP[i]);
        totalGammaProb = LogSumExp(totalGammaProb, sumGamma[i]);

        // sum all alpha*beta
        sumStateGamma[i] = LogSumExp(sumStateGamma[i], sumGamma[i]);
        newEmissionProb[i][genomeObs] = LogSumExp(newEmissionProb[i][genomeObs], sumGamma[i]);

        // assign alhpa*beta to xi
        for (int j = 0 ; j < state ; ++j)
        {
            sumXi[i][j] = forwardP[i] + trans[state*i+j] + emission[4*j+genomePos] + backwardP[state*1+j];
            totalXiProb = LogSumExp(totalXiProb, sumXi[i][j]) ;
        }
    }

    for (int i = 0 ; i < state ; ++i)
    {
        // gamma/sum => log2(gamma/sum) = log2(gamma) - log2(sum)
        // compute gamma1 prob and update initState 
        sumGamma[i] = sumGamma[i] - totalGammaProb;
        initState[i] = pow(2, sumGamma[i]);


        // compute xi1 prob
        for (int j = 0 ; j < state ; ++j)
        {
            sumXi[i][j] = sumXi[i][j] - totalXiProb;
        }
    }


    // compute gamma and xi, i from 1 to T-2, seq from 0 to T-1
    for (int i = 1 ; i < seqSize-1 ; ++i)
    {
        totalGammaProb = -INFINITY;
        totalXiProb = -INFINITY;

        genomePos = decodeDNAHMM[sTraining[i+1]];
        genomeObs = decodeDNAHMM[sTraining[i]];
        for (int j = 0 ; j < state ; ++j)
        {
            // forward positive and backward positive (log2)
            // assign alpha*beta to gamma
            tempGamma[j] = (forwardP[state*i+j] + backwardP[state*i+j]);
            totalGammaProb = LogSumExp(totalGammaProb, tempGamma[j]);

            // sum all alpha*beta
            sumStateGamma[j] = LogSumExp(sumStateGamma[j], tempGamma[j]);
            newEmissionProb[j][genomeObs] = LogSumExp(newEmissionProb[j][genomeObs], tempGamma[j]);
            
            // assign alhpa*beta to xi
            for (int k = 0 ; k < state ; ++k)
            {
                tempXi[j][k] = forwardP[state*i+j] + trans[state*j+k] + emission[4*k+genomePos] + backwardP[(state*(i+1))+k];
            }
        }

        // compute total xi prob
        for (int j = 0 ; j < state ; ++j)
        {
            for (int k = 0 ; k < state ; ++k)
            {
                totalXiProb = LogSumExp(totalXiProb, tempXi[j][k]);
            }
        }


        for (int j = 0 ; j < state ; ++j)
        {
            // gamma/sum => log2(gamma/sum) = log2(gamma) - log2(sum)
            // compute gammai prob
            sumGamma[j] = LogSumExp(sumGamma[j], (tempGamma[j] - totalGammaProb));

            // compute xii prob
            for (int k = 0 ; k < state ; ++k)
            {
                sumXi[j][k] = LogSumExp(sumXi[j][k], (tempXi[j][k] - totalXiProb));
            }
        }

    }

    // add final value
    genomeObs = decodeDNAHMM[sTraining[(seqSize-1)]];
    for (int i = 0 ; i < state ; ++i)
    {
        tempGamma[i] = (forwardP[state*(seqSize-1)+i] + backwardP[state*(seqSize-1)+i]);
        // sum all alpha*beta
        sumStateGamma[i] = LogSumExp(sumStateGamma[i], tempGamma[i]);
        newEmissionProb[i][genomeObs] = LogSumExp(newEmissionProb[i][genomeObs], tempGamma[i]);
    }

    // update transition probability matrix
    for (int i = 0 ; i < state ; ++i)
    {
        for (int j = 0 ; j < state ; ++j)
        {
            // sum xi(i,j)/ sum gamma(i), t from 0 to T-2, 
            newTransProb[i][j] = pow(2, (sumXi[i][j] - sumGamma[i]));
        }
    }


    // update emission probability matrix
    for (int i = 0; i < state; ++i) 
    {
        for (int j = 0; j < 4; ++j) 
        {
            emission[state*i + j] = pow(2, (newEmissionProb[i][j] - sumStateGamma[i]));
        }
    }

    // compute error 
    diffSquare = 0.0;

    for (int i = 0 ; i < state ; ++i)
    {
        for (int j = 0 ; j < state ; ++j)
        {
            // add all difference of square 
            diffSquare += pow((newTransProb[i][j]-trans[state*i+j]), 2);
        }
    }

    diffSquare = sqrt(diffSquare);
    // printf("%d => %f\n", state, diffSquare);
    // if error <= threshold
    if (diffSquare <= threshold)
        *bConvergen = true;

    
    // update initProb 
    for (int i = 0 ; i < state ; ++i)
    {
        init[i] = initState[i];
    }

    // update transProb
    for (int i = 0 ; i < state ; ++i)
    {
        for (int j = 0 ; j < state ; ++j)
        {
            trans[state*i+j] = newTransProb[i][j];
        }
    }


    // normalize probability matrix 
    rowSum = 0.0;

    // init
    for (int i = 0 ; i < state ; ++i)
        rowSum += init[i];
    
    for (int i = 0 ; i < state ; ++i)
        init[i] = init[i]/rowSum;
    
    // transition
    for (int i = 0 ; i < state ; ++i)
    {
        rowSum = 0.0;

        for (int j = 0 ; j < state ; ++j)
        {
            rowSum += trans[state*i+j];
        }

        for (int j = 0 ; j < state ; ++j)
        {
            trans[state*i+j] = trans[state*i+j]/rowSum;
        }
    }

    // emission
    for (int i = 0 ; i < state ; ++i)
    {
        rowSum = 0.0;

        for (int j = 0 ; j < 4 ; ++j)
        {
            rowSum += emission[state*i+j];
        }

        for (int j = 0 ; j < 4 ; ++j)
        {
            emission[state*i+j] = emission[state*i+j]/rowSum;
        }
    }


    delete []initState;
    initState = nullptr;

    delete []sumGamma;
    sumGamma = nullptr;

    delete []tempGamma;
    tempGamma = nullptr;

    delete []sumStateGamma;
    sumStateGamma = nullptr;

    // allocate 2D dynamic array
    for (int i = 0 ; i < state ; ++i)
    {
        delete [] sumXi[i];
        delete [] tempXi[i];
        delete [] newTransProb[i];
        delete [] newEmissionProb[i];
    }

    delete []sumXi;
    sumXi = nullptr;

    delete []tempXi;
    tempXi = nullptr;

    delete []newTransProb;
    newTransProb = nullptr;

    delete []newEmissionProb;
    newEmissionProb = nullptr;
    

    return;
}  


static void PrintUpdateInfoMatrix(int state, double *init, double *trans, double *emission)
{
    // print initProb
    std::cout << "Convergen " << state << " State Initial Probability = \n";

    for (int i = 0 ; i < state ; ++i)
    {
        if (i == 0)
            std::cout << "{" << std::right << std::fixed << std::setprecision(6) << std::setw(5) << init[i] << ", ";
        else if (i == (state-1))
            std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << init[i] << "}\n";
        else
            std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << init[i] << ", ";
    }


    // print transProb
    std::cout << "\nConvergen " << state << " State Transition Probability = \n";
    for (int i = 0 ; i < state ; ++i)
    {
        for (int j = 0 ; j < state ; ++j)
        {
            if (j == 0)
                std::cout << "{" << std::right << std::fixed << std::setprecision(6) << std::setw(5) << trans[state*i+j] << ", ";
            else if (j == (state-1))
                std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << trans[state*i+j] << "}\n";
            else
                std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << trans[state*i+j] << ", ";
        }
    }

    // print emissionProb
    std::cout << "\nConvergen " << state << " State Emission Probability = \n";
    for (int i = 0 ; i < state ; ++i)
    {
        for (int j = 0 ; j < 4 ; ++j)
        {
            if (j == 0)
                std::cout << "{" << std::right << std::fixed << std::setprecision(6) << std::setw(5) << emission[state*i+j] << ", ";
            else if (j == 3)
                std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << emission[state*i+j] << "}\n";
            else
                std::cout << std::right << std::fixed << std::setprecision(6) << std::setw(5) << emission[state*i+j] << ", ";
        }
    }

    printf("========================================================\n");

    return;

}

static void ComputeFinishModelLikelihood(int state, double *modelProb, double *init, double *trans, double *emission, GenomeSequence &genomeSeq)
{
    std::string sTraining = genomeSeq.straining;
    std::string sTest1 = genomeSeq.stest1;
    std::string sTest2 = genomeSeq.stest2;

    double *alphaTP = new double[state];
    double *talphaTP = new double[state];
    double *nextTP = new double[state];

    double *alpha1P = new double[state];
    double *talpha1P = new double[state];
    double *next1P = new double[state];

    double *alpha2P = new double[state];
    double *talpha2P = new double[state];
    double *next2P = new double[state];

    int location = 0;
    double tTrainProb = -INFINITY;
    double tTest1Prob = -INFINITY;
    double tTest2Prob = -INFINITY;
    
    int genomeTrainPos ;
    int genomeTest1Pos ;
    int genomeTest2Pos ;
        
    // Positive strand
    // step 1 : a1 = p1 * b1(O1)
    for (int i = 0 ; i < state ; ++i)
    {
        location = 4*i;
        alphaTP[i] = log2(init[i]) + log2(emission[location+decodeDNAHMM[sTraining[0]]]);
        alpha1P[i] = log2(init[i]) + log2(emission[location+decodeDNAHMM[sTest1[0]]]);
        alpha2P[i] = log2(init[i]) + log2(emission[location+decodeDNAHMM[sTest2[0]]]);
    }

    // step2 : ai = sum(ai-1(i)* aij * bj(Oi))
    for (int i = 1 ; i < seqSize ; ++i)
    {
        genomeTrainPos = decodeDNAHMM[sTraining[i]];
        genomeTest1Pos = decodeDNAHMM[sTest1[i]];
        genomeTest2Pos = decodeDNAHMM[sTest2[i]];

        for (int j = 0 ; j < state ; ++j)
        {
            talphaTP[j] = -INFINITY;
            talpha1P[j] = -INFINITY;
            talpha2P[j] = -INFINITY;

            for (int k = 0 ; k < state ; ++k)
            {
                // Sk->Sj
                tTrainProb = alphaTP[k] + log2(trans[state*k+j]) + log2(emission[4*j+genomeTrainPos]);
                tTest1Prob = alpha1P[k] + log2(trans[state*k+j]) + log2(emission[4*j+genomeTest1Pos]);
                tTest2Prob = alpha2P[k] + log2(trans[state*k+j]) + log2(emission[4*j+genomeTest2Pos]);
                talphaTP[j] = LogSumExp(talphaTP[j], tTrainProb);
                talpha1P[j] = LogSumExp(talpha1P[j], tTest1Prob);
                talpha2P[j] = LogSumExp(talpha2P[j], tTest2Prob);
            }

            nextTP[j] = talphaTP[j];
            next1P[j] = talpha1P[j];
            next2P[j] = talpha2P[j];
        }

        for (int j = 0 ; j < state ; ++j)
        {
            alphaTP[j] = nextTP[j];
            alpha1P[j] = next1P[j];
            alpha2P[j] = next2P[j];
        }
    }  

    tTrainProb = -INFINITY;
    tTest1Prob = -INFINITY;
    tTest2Prob = -INFINITY;

    for (int i = 0 ; i < state ; ++i)
    {
        tTrainProb = LogSumExp(alphaTP[i], tTrainProb);
        tTest1Prob = LogSumExp(alpha1P[i], tTest1Prob);
        tTest2Prob = LogSumExp(alpha2P[i], tTest2Prob);
    }

    modelProb[0] = tTrainProb;
    modelProb[1] = tTest1Prob;
    modelProb[2] = tTest2Prob;

    delete []alphaTP;
    alphaTP = nullptr;
    delete []alpha1P;
    alpha1P = nullptr;
    delete []alpha2P;
    alpha2P = nullptr;

    delete []talphaTP;
    talphaTP = nullptr;
    delete []talpha1P;
    talpha1P = nullptr;
    delete []talpha2P;
    talpha2P = nullptr;

    delete []nextTP;
    nextTP = nullptr;
    delete []next1P;
    next1P = nullptr;
    delete []next2P;
    next2P = nullptr;    
    
    return;

}

void PrintAllSeqLogLikelihood()
{
    printf("\n  Model           Training           Test1           Test2\n");
    std::cout << std::right << std::setw(7) << "state 3   ";
    for (int i = 0 ; i < 3 ; ++i)
        std::cout << std::right << std::scientific << std::setw(16) << model3Prob[i] ;
    
    std::cout << "\n" << std::right << std::setw(7) << "state 5   ";
    for (int i = 0 ; i < 3 ; ++i)
        std::cout << std::right << std::scientific << std::setw(16) << model5Prob[i] ;
    
    std::cout << std::endl ;

    return;
}

void TrainEMHMM(int state, GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq)
{

    InitialDecodeDNA();
    
    if (state == 3)
    {
        int count = 0;
        bool bConvergen = false;

        while (bConvergen == false)
        {
            ++count;

            std::thread forward(Forward, state, init3SProb, (double *)trans3SProb, (double *)emission3bProb, std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)alpha3State, (double *)alpha3StateRev);
            std::thread backward(Backward, state, init3SProb, (double *)trans3SProb, (double *)emission3bProb, std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)beta3State, (double *)beta3StateRev);

            forward.join();
            backward.join();

            ComputeGammaXiAndUpdateModel(state, init3SProb, (double *)trans3SProb, (double *)emission3bProb,  std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)alpha3State, (double *)alpha3StateRev, (double *)beta3State, (double *)beta3StateRev, &bConvergen);
        }

        while (bCanPrint3State == true)
        {
            printf("Count = %d\n", count);
            PrintUpdateInfoMatrix(state, init3SProb, (double *)trans3SProb, (double *)emission3bProb);
            bCanPrint3State = false;
        }

        ComputeFinishModelLikelihood(state, model3Prob, init3SProb, (double *)trans3SProb, (double *)emission3bProb, genomeSeq);
    }
    else if (state == 5)
    {
        int count = 0;
        bool bConvergen = false;

        while (bConvergen == false)
        {
            ++count;
            std::thread forward(Forward, state, init5SProb, (double *)trans5SProb, (double *)emission5bProb, std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)alpha5State, (double *)alpha5StateRev);
            std::thread backward(Backward, state, init5SProb, (double *)trans5SProb, (double *)emission5bProb, std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)beta5State, (double *)beta5StateRev);

            forward.join();
            backward.join();

            ComputeGammaXiAndUpdateModel(state, init5SProb, (double *)trans5SProb, (double *)emission5bProb,  std::ref(genomeSeq.straining), std::ref(genomeDSeq.straining), (double *)alpha5State, (double *)alpha5StateRev, (double *)beta5State, (double *)beta5StateRev, &bConvergen);
        }

        while (bCanPrint3State == true);

        printf("Count = %d\n", count);
        PrintUpdateInfoMatrix(state, init5SProb, (double *)trans5SProb, (double *)emission5bProb);
        bCanPrint3State = true;

        ComputeFinishModelLikelihood(state, model5Prob, init5SProb, (double *)trans5SProb, (double *)emission5bProb, genomeSeq);
        
    }
    
    return;
}




