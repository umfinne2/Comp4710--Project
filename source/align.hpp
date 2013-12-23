#ifndef ALIGN_H
#define ALIGN_H

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

class Align
{
    public:
        TAlign align;
        TSequence ref_seq;
        TSequence read_seq;

        Align(void);
        Align(TAlign align, TSequence ref_seq, TSequence read_seq);
        int run(void);
        virtual ~Align(void);

    protected:
        float max_index[2];     //stores the starting index of traceback
        float **matrix;
        int len1;
        int len2;

        void createMatrix(void);
        virtual void initMatrix(void) = 0;
        virtual float calcScore(int i, int j) = 0;
        virtual float calcElement(int i, int j) = 0;
        void traceback(void);

}
#endif
