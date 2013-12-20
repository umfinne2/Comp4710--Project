#ifndef LOCAL_ALIGNMENT_H
#define LOCAL_ALIGNMENT_H

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/align.h>

#define DIAGONAL 0
#define LEFT 1
#define TOP 2

using namespace seqan;

typedef String<char> TSequence;
typedef Align<TSequence,ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

class LocalAlignment {

  private:
    static int match_cost(TSequence seq1, TSequence seq2, int i, int j, Score<int, Simple> scores);
    static int max(int * array, int size);

  public:
    static int smith_waterman(TSequence seq1, TSequence seq2, Score<int, Simple> scores);
};

#endif
