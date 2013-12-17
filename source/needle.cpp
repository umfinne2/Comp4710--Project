#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

typedef String<char> TSequence;
typedef Align<TSequence, ArrayGaps> TAlign;
typedef Row<TAlign>::Type TRow;

//ugly max function cause the std one was error out on me for some reason (Rory)
float max(float x, float y, float z)
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

int needle( TAlign &align,
            TSequence ref_seq, 
            TSequence read_seq, 
            Score<int, Simple> scheme)
{
    int score = 0;
    //possibly inefficient copying object rather than using reference
    //needle(align, ref_seq, read_seq);
    resize( rows(align), 2 );
    cout << "Resized alignment\n";

    assignSource( row( align, 0 ), ref_seq);
    assignSource( row( align, 1 ), read_seq);
    cout << "Assigned source sequences\n";

    //use the built in one just to make sure it works
    //we don't actually want to use this cause they may have
    //optimization we don't care about. 
    //score = globalAlignment(align, scheme, NeedlemanWunsch());    
    
    //get the lengths of the sequences
    int len1 = length(ref_seq);
    int len2 = length(read_seq);

    cout << "Got sequence lengths\n";

    //create the DP Matrix/Table
    float **matrix = new float * [len1 + 1];
    for (int i = 0; i <= len1; i++)
    {
        matrix[i] = new float[len2 + 1];
    }

    cout << "Created matrix\n";
        
    //initialize (will change by implementation)
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
            //cout << "i=" << i << " j=" << j << "\n";
            diagonal = matrix[i-1][j-1];
            vertical = matrix[i][j-1] + scoreGap(scheme);
            horizontal = matrix[i-1][j] + scoreGap(scheme);

            //cout << "ref[i-1=" << i-1 << "] = " << ref_seq[i-1] << "\n";
            //not sure if elements in a Dna5String can be compared like this
            if (toupper(ref_seq[i-1]) == toupper(read_seq[j-1]))
            {
                //cout << "Match at i=" << i << " j=" << j << "\n";
                diagonal += scoreMatch(scheme);
            }

            else
            {
                diagonal += scoreMismatch(scheme);
            }

            matrix[i][j] = max(diagonal, vertical, horizontal);
        }

        //cout << "Diagonal=" << diagonal << " Vertical=" << vertical << " Horizontal=" << horizontal << "\n";
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
    while ( i > 0 || j > 0)
    {
        pos = matrix[i][j];
        match = matrix[i-1][j-1] + scoreMatch(scheme);
        mismatch = matrix[i-1][j-1] + scoreMismatch(scheme);
        vgap = matrix[i][j-1] + scoreGap(scheme);
        hgap = matrix[i-1][j] + scoreGap(scheme);

        if (j > 0 && pos == vgap)
        {
            //cout << "Inserting a gap in the reference sequence\n";
            insertGap(row1, i);
            j--;
        }

        else if (i > 0 && pos == hgap)
        {
            //cout << "Inserting a gap in the read\n";
            insertGap(row2, j);
            i--;
        }

        else if (i > 0 && (pos == match || pos == mismatch))
        {
            //cout << "Mapping bp from read to ref\n";
            i--;
            j--;
        }

    }

    //see if printing here works any better
    cout << align << "\n";

    //delete the dp matrix
    for (i = 0; i <= len1; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;
    
    return score;
}

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        cout << "Error: Invalid number of arguments\n";
        return 1;
    }

    fstream in_ref(argv[1], ios::binary | ios::in);
    fstream in_reads(argv[2], ios::binary | ios::in);
    RecordReader<fstream, SinglePass<> > reader_ref(in_ref);
    RecordReader<fstream, SinglePass<> > reader_reads(in_reads);

    CharString ref_id;
    TSequence ref_seq;
    CharString read_id;
    CharString read_qual;
    TSequence read_seq;

    TAlign align;

    if (!atEnd(reader_ref) && readRecord(ref_id, ref_seq, reader_ref, Fasta()) != 0)
    {
        cout << "Error: Couldn't read reference sequence\n";
        return 1;
    }

    cout << ref_id << "\t" << ref_seq << "\n";

    while (!atEnd(reader_reads))
    {
        if (readRecord(read_id, read_seq, read_qual, reader_reads, Fastq()) != 0)
        {
            cout << "Error: Couldn't read the sequence from the file\n";
            return 1;
        }

        Score<int, Simple> scoringScheme(0, -1, -1);
        int score = needle(align, ref_seq, read_seq, scoringScheme);
        
        //cout << align << "\n"; 

    }

    return 0;
}


