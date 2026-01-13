#include "genome.h"

static void BuildStrandGenomeArr()
{
    for (int i = 0 ; i < 128 ; ++i)
    {
        if (i == 'A')
        {
            strandGenomeArr[i]='T';
        }
        else if (i == 'T')
        {
            strandGenomeArr[i]='A';
        }
        else if (i == 'C')
        {
            strandGenomeArr[i]='G';
        }
        else if (i == 'G')
        {
            strandGenomeArr[i]='C';
        }  
        else
            strandGenomeArr[i]=i;
    }

    return;
}


void FindGenomeSequence(std::ifstream &myfile, GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq)
{
    startPos = waiting;

    std::string str = "";

    char ch;

    long long charCnt = 0;


    BuildStrandGenomeArr();

    // if getline == eof, then break
    while (myfile.get(ch))
    {
        if (ch == 'c')
        {
            str = "";
            str += ch;
                
            while (myfile.get(ch) && ch != ' ')
            {
                str += ch;
            }
            
            while (myfile.get(ch) && ch != '\n');
            // printf("%s\n", str.c_str());

            if ((str == "chr1") && (startPos == waiting))
            {
                startPos = chr1Pos;
                str = "";
                charCnt = 0;
            }

            if ((str == "chr2") && (startPos == waiting))
            {
                startPos = chr2Pos;
                str = "";
                charCnt = 0;
            }
                
        }

        if (((startPos == chr1Pos) || (startPos == trainingPosition) || (startPos == chr2Pos)) 
             && ((ch == 'A') || (ch == 'T') || (ch == 'C') || (ch == 'G') || (ch == 'N')))
        {
            ++charCnt;
            if ((startPos == chr1Pos) && (charCnt == trainingStratPos))
            {
                startPos = trainingPosition;
                --charCnt;

                do{
                    if ((ch == 'A') || (ch == 'T') || (ch == 'C') || (ch == 'G') || (ch == 'N'))
                    {
                        ++charCnt;
                        genomeSeq.straining += ch;
                        genomeDSeq.straining += strandGenomeArr[ch];
                    }
                } while ((myfile.get(ch)) && (charCnt != trainingEndPos));

                reverse(genomeDSeq.straining.begin(), genomeDSeq.straining.end());

                printf("Traning = \n");
                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeSeq.straining[i]);
                }
                printf("\n");

                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeDSeq.straining[i]);
                }

                printf("\n");
            }
                
            if ((startPos == trainingPosition) && (charCnt == test1StartPos))
            {
                startPos = test1Position;
                --charCnt;

                do{
                    if ((ch == 'A') || (ch == 'T') || (ch == 'C') || (ch == 'G') || (ch == 'N'))
                    {
                        ++charCnt;
                        genomeSeq.stest1 += ch;
                        genomeDSeq.stest1 += strandGenomeArr[ch];
                    }
                } while ((myfile.get(ch)) && (charCnt != test1EndPos));

                reverse(genomeDSeq.stest1.begin(), genomeDSeq.stest1.end());

                printf("stest1 = \n");
                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeSeq.stest1[i]);
                }
                printf("\n");

                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeDSeq.stest1[i]);
                }
                printf("\n");

                startPos = waiting;
            }

            if ((startPos == chr2Pos) && (charCnt == test2StartPos))
            {
                startPos = test2Position;
                --charCnt;

                do{
                    if ((ch == 'A') || (ch == 'T') || (ch == 'C') || (ch == 'G') || (ch == 'N'))
                    {
                        ++charCnt;
                        genomeSeq.stest2 += ch;
                        genomeDSeq.stest2 += strandGenomeArr[ch];
                    }
                } while ((myfile.get(ch)) && (charCnt != test2EndPos));

                reverse(genomeDSeq.stest2.begin(), genomeDSeq.stest2.end());
                
                printf("stest2 = \n");
                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeSeq.stest2[i]);
                }
                printf("\n");

                for (int i = 0 ; i < 70 ; ++i)
                {
                    printf("%c", genomeDSeq.stest2[i]);
                }
                printf("\n");

                startPos = finish;
                break;
            }
       }
            

    }

    return ;

}
