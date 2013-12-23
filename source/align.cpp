#include "align.hpp"

using namespace seqan;
using namespace std;

Align::Align(void)
{
}

Align::Align(TAlign align, TSequence ref_seq, TSequence read_seq)
{
    this.align = align;
    this.ref_seq = ref_seq;
    this.read_seq = read_seq;
}

int Align::run(void)
{
    float score;

    //allocating memory for the matrix shouldn't change
    createMatrix(void);

    //this can change
    initMatrix(void);

    for (i = 1; i <= len1; i++)
    {
        for (j = 1; j <= len2; j++)
        {
            //since smith-waterman needs to keep track of max
            matrix[i][j] = this.calcElement(i, j);
        }
    }

    score = matrix[max_index[0]][max_index[1]];

    traceback(void);

    return score;
}

void traceback(void)
{

}
