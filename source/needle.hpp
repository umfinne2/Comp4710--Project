#ifndef NEEDLE_H
#define NEEDLE_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

float max(float x, float y, float z);
int needle( TAlign &align, TSequence ref_seq, TSequence read_seq, seqan::Score<int, seqan::Simple> scheme);

#endif
