#ifndef GENOME_H
#define GENOME_H

#include <iostream>
#include <fstream>
#include <string>

#define trainingStratPos 100000000
#define test1StartPos    145000000-1
#define test2StartPos    100000000-1

#define segSize          5000000

#define trainingEndPos   trainingStratPos+segSize-1
#define test1EndPos      test1StartPos+segSize-1
#define test2EndPos      test2StartPos+segSize-1


enum GenomeStartPosition{
    waiting,
    chr1Pos,
    trainingPosition,
    test1Position,
    chr2Pos,
    test2Position,
    finish
};

struct GenomeSequence{
    std::string straining;
    std::string stest1;
    std::string stest2;

    GenomeSequence() : straining(""), stest1(""), stest2("") {}
};

struct DoubleStranded{
    std::string straining;
    std::string stest1;
    std::string stest2;

    DoubleStranded() : straining(""), stest1(""), stest2("") {}
};

static char strandGenomeArr[130];

static GenomeStartPosition startPos;

static void BuildStrandGenomeArr();
void FindGenomeSequence(std::ifstream &myfile, GenomeSequence &genomeSeq, DoubleStranded &genomeDSeq);

#endif