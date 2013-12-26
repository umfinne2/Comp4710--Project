#include <fstream>
#include <iostream>

//#include <seqan/seq_io.h>
//#include <seqan/sequence.h>

#include "global_alignment.hpp"
#include "align_lib.hpp"

using namespace std;
using namespace seqan;

//typedef String<char> TSequence;
//typedef Align<TSequence, ArrayGaps> TAlign;
//typedef Row<TAlign>::Type TRow;

//ugly max function cause the std one was error out on me for some reason (Rory)
float GlobalAlignment::max(float x, float y, float z)
{
    if (x >= y && x >= z)
    {
        return x;
    }
    else if (y >= z)
    {
        return y;
    }
    else
    {
        return z;
    }
}

int GlobalAlignment::needle(    TAlign &align,
                                TSequence ref_seq,
                                TSequence read_seq)
{
    int score = 0;
    resize( rows(align), 2 );
    assignSource( row( align, 0 ), ref_seq);
    assignSource( row( align, 1 ), read_seq);
    int len1 = length(ref_seq);
    int len2 = length(read_seq);

    //create the DP Matrix/Table + could probably be its own function
    float **matrix = new float * [len1 + 1];
    for (int i = 0; i <= len1; i++)
    {
        matrix[i] = new float[len2 + 1];
    }

    //initialize (changes based on alg. could probably use function callback)
    for (int i = 0; i <= len1; i++)
    {
        matrix[i][0] = i * AlignLib::gapcost;
    }

    for (int j = 0; j <= len2; j++)
    {
        matrix[0][j] = j * AlignLib::gapcost;
    }

    float diagonal, vertical, horizontal;
    int i, j;
    for (i = 1; i <= len1; i++)
    {
        for (j = 1; j <= len2; j++)
        {
            //we may want to change this if we use a lookup table rather than simple values
            //separate function???
            diagonal = matrix[i-1][j-1] + AlignLib::get_score(ref_seq[i-1], read_seq[j-1]);
            vertical = matrix[i][j-1] + AlignLib::gapcost;
            horizontal = matrix[i-1][j] + AlignLib::gapcost;

            matrix[i][j] = max(diagonal, vertical, horizontal);
        }
    }

    //set i and j to corner index of matrix
    i--;
    j--;

    score = matrix[i][j];

    //print out the top corner of the matrix and sequences
    //AlignLib::print_matrix(matrix, ref_seq, read_seq, len2+1, len1+1, -10, 10);

    //traceback
    float pos, dmap, vgap, hgap;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    while ( i > 0 && j > 0 )
    {
        //if we write a function for computing values than this traceback routine should never change
        //between implementations
        pos = matrix[i][j];
        dmap = matrix[i-1][j-1] + AlignLib::get_score(ref_seq[i-1], read_seq[j-1]);
        //match = matrix[i-1][j-1] + scoreMatch(scheme);
        //mismatch = matrix[i-1][j-1] + scoreMismatch(scheme);
        vgap = matrix[i][j-1] + AlignLib::gapcost;
        hgap = matrix[i-1][j] + AlignLib::gapcost;

        if (pos == dmap)
        {
            i--;
            j--;
        }

        else if (pos == hgap)
        {
            insertGap(row2, j);
            i--;
        }

        else if (pos == vgap)
        {
            insertGap(row1, i);
            j--;
        }

        else
        {
            cout << "You distroyed the universe! You're on your own now :( \n";
            exit(1);
        }
    }

    while (i > 0)
    {
        insertGap(row2, j);
        i--;
    }

    while (j > 0)
    {
        insertGap(row1, i);
        j--;
    }

    //delete the dp matrix
    for (i = 0; i <= len1; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;

    return score;
}
