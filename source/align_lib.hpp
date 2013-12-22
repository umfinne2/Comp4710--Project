#ifndef ALIGN_LIB_H
#define ALIGN_LIB_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

float percent_match(TAlign &correct, TAlign &test);

#endif
