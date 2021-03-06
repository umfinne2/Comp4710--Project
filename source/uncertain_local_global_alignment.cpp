#include <fstream>
#include <iostream>

#include "uncertain_local_global_alignment.hpp"
#include "align_lib.hpp"

using namespace std;
using namespace seqan;

// function based on description @ http://en.wikipedia.org/wiki/Phred_quality_score
float UncertainLocalGlobalAlignment::qscore_to_percent(char qscore)
{
    //P = 10**(-Q/10)
    int q = (int)qscore;
    float p = (float)pow(10, -(q/10));
    float result = 1.0 - p;

    return result;
}

//ugly max function cause the std one was error out on me for some reason (Rory)
float UncertainLocalGlobalAlignment::max(float x, float y, float z)
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

int UncertainLocalGlobalAlignment::ulga(TAlign &align,
                                        TSequence ref_seq,
                                        TSequence read_seq,
                                        CharString read_qual)
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

    //initialize top row of matrix to 0s like local alignment
    for (int i = 0; i <= len1; i++)
    {
        matrix[i][0] = 0.0;
    }

    for (int j = 0; j <= len2; j++)
    {
        matrix[0][j] = j * AlignLib::gapcost;
    }

    int max_i, max_j = 0;
    float diagonal, vertical, horizontal;
    int i, j;
    for (i = 1; i <= len1; i++)
    {
        for (j = 1; j <= len2; j++)
        {
            diagonal = matrix[i-1][j-1] + (AlignLib::get_score(ref_seq[i-1], read_seq[j-1]) * qscore_to_percent(read_qual[j-1]));
            vertical = matrix[i][j-1] + AlignLib::gapcost;
            horizontal = matrix[i-1][j] + AlignLib::gapcost;

            matrix[i][j] = max(diagonal, vertical, horizontal);
        }

        //max_j should always be the length of the columns
        //keep track of max start column in bottom row.
        if (max_j != 0)
        {
            if (matrix[i][j-1] >= matrix[max_i][max_j])
            {
                max_i = i;
                max_j = j-1;
            }
        }
        else
        {
            max_i = i;
            max_j = j-1;
        }
    }

    //set i and j to corner index of matrix
    i--;
    j--;

    score = matrix[max_i][max_j];

    //traceback
    float pos, dmap, vgap, hgap;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    insertGaps(row2, j, (i - max_i));   //insert tailing gaps
    i = max_i;                          //start column iteration at index max_i
    while ( i > 0 && j > 0 )
    {
        //if we write a function for computing values than this traceback routine should never change
        //between implementations
        pos = matrix[i][j];
        dmap = matrix[i-1][j-1] + (AlignLib::get_score(ref_seq[i-1], read_seq[j-1]) * qscore_to_percent(read_qual[j-1]));
        //mismatch = matrix[i-1][j-1] + scoreMismatch(scheme);
        vgap = matrix[i][j-1] + AlignLib::gapcost;
        hgap = matrix[i-1][j] + AlignLib::gapcost;

        if (i > 0 && j > 0 && pos == dmap)
        {
            i--;
            j--;
        }

        else if (i > 0 && pos == hgap)
        {
            insertGap(row2, j);
            i--;
        }

        else if (j > 0 && pos == vgap)
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

    //maybe check that j == 0 and i >= 0
    insertGaps(row2, 0, (i));

    //delete the dp matrix
    for (i = 0; i <= len1; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;

    return score;
}
