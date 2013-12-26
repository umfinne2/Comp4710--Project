#ifndef ALIGN_LIB_H
#define ALIGN_LIB_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

//std::map<char, int> dnamap;

//NOTE: we will only be using the 4 nucleotide values, but this is the full matrix available from ncbi
//std::map<char,int> dnamap;
class AlignLib
{
    private:
        static int get_index(char x);

    public:
        static const int gapcost = -2;
        static int get_score(char a, char b);
        static float percent_match(TAlign &correct, TAlign &test);
        static void print_matrix(float **matrix, TSequence seq1, TSequence seq2, int rows, int cols, int offset, int len);
};

#endif
