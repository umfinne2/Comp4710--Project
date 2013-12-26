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
    static float max(float * array, int size);

  public:
    static int smith_waterman(TAlign &align, TSequence seq1, TSequence seq2);
};

#endif
