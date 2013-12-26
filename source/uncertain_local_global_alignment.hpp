#ifndef UNCERTAIN_LOCAL_GLOBAL_ALIGNMENT_H
#define UNCERTAIN_LOCAL_GLOBAL_ALIGNMENT_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

typedef seqan::String<char> TSequence;
typedef seqan::Align<TSequence, seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

class UncertainLocalGlobalAlignment
{
    private:
        static float qscore_to_percent(char qscore);
        static float max(float x, float y, float z);

    public:
        static int ulga( TAlign &align, TSequence ref_seq, TSequence read_seq, seqan::CharString read_qual);
};
#endif
