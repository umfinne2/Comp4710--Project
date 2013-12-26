#ifndef GLOBAL_ALIGNMENT_H
#define GLOBAL_ALIGNMENT_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

class GlobalAlignment
{
    private:
        static float max(float x, float y, float z);

    public:
        static int needle( TAlign &align, TSequence ref_seq, TSequence read_seq);
};
#endif
