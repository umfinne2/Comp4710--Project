#include <fstream>
#include <iostream>

//#include <seqan/seq_io.h>
//#include <seqan/sequence.h>

#include "global_alignment.hpp"

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
                                TSequence read_seq,
                                Score<int, Simple> scheme)
{
    int score = 0;
    resize( rows(align), 2 );

    assignSource( row( align, 0 ), ref_seq);
    assignSource( row( align, 1 ), read_seq);

    //use the built in one just to make sure it works
    //we don't actually want to use this cause they may have
    //optimization we don't care about.
    //score = globalAlignment(align, scheme, NeedlemanWunsch());

    //get the lengths of the sequences
    int len1 = length(ref_seq);
    int len2 = length(read_seq);

    cout << "Got sequence lengths\n";

    //create the DP Matrix/Table + could probably be its own function
    float **matrix = new float * [len1 + 1];
    for (int i = 0; i <= len1; i++)
    {
        matrix[i] = new float[len2 + 1];
    }

    cout << "Created matrix\n";

    //initialize (changes based on alg. could probably use function callback)
    for (int i = 0; i <= len1; i++)
    {
        matrix[i][0] = i * scoreGap(scheme);
    }

    for (int j = 0; j <= len2; j++)
    {
        matrix[0][j] = j * scoreGap(scheme);
    }

    cout << "Initialized the matrix\n";
    cout << "Gap=" << scoreGap(scheme) << " Match=" << scoreMatch(scheme) << " Mismatch=" << scoreMismatch(scheme) << "\n";
    //our values for storing each potential movement
    float diagonal, vertical, horizontal;
    int i, j;
    for (i = 1; i <= len1; i++)
    {
        for (j = 1; j <= len2; j++)
        {
            //we may want to change this if we use a lookup table rather than simple values
            //separate function???
            diagonal = matrix[i-1][j-1];
            vertical = matrix[i][j-1] + scoreGap(scheme);
            horizontal = matrix[i-1][j] + scoreGap(scheme);

            if (toupper(ref_seq[i-1]) == toupper(read_seq[j-1]))
            {
                diagonal += scoreMatch(scheme);
            }

            else
            {
                diagonal += scoreMismatch(scheme);
            }

            matrix[i][j] = max(diagonal, vertical, horizontal);
        }
    }

    //cout << align << "\n";

    //set i and j to corner index of matrix
    i--;
    j--;

    score = matrix[i][j];
    cout << "Score=" << score << "\n";

    //print out the top corner of the matrix and sequences
    /*
    cout << setw(5) << " ";
    for (int l = 1; l < 30; l++)
    {
        cout << setw(5) << ref_seq[l-1];
    }
    cout << endl;

    for (int k = 0; k < 30; k++)
    {
        if (k > 0)
        {
            cout << read_seq[k-1];
        }

        for (int l = 0; l < 30; l++)
        {
            cout << setw(5) << matrix[l][k];
        }
        cout << endl;
    }
    */
    //traceback
    float pos, match, mismatch, vgap, hgap;
    TRow &row1 = row(align,0);
    TRow &row2 = row(align,1);
    while ( i > 0 || j > 0 )
    {
        //if we write a function for computing values than this traceback routine should never change
        //between implementations
        pos = matrix[i][j];
        match = matrix[i-1][j-1] + scoreMatch(scheme);
        mismatch = matrix[i-1][j-1] + scoreMismatch(scheme);
        vgap = matrix[i][j-1] + scoreGap(scheme);
        hgap = matrix[i-1][j] + scoreGap(scheme);

        if (i > 0 && pos == hgap)
        {
            insertGap(row2, j);
            i--;
        }

        else if (j > 0 && pos == vgap)
        {
            insertGap(row1, i);
            j--;
        }

        else if (i > 0 && j > 0 && (pos == match || pos == mismatch))
        {
            i--;
            j--;
        }

        else
        {
            cout << "You distroyed the universe! You're on your own now :( \n";
            exit(1);
        }
    }

    //see if printing here works any better
    //cout << align << "\n";

    //delete the dp matrix
    for (i = 0; i <= len1; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;

    return score;
}
