#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace std;
using namespace seqan;

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

int needle( Align<Dna5String> &align,
            Dna5String ref_seq, 
            Dna5String read_seq, 
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
    score = globalAlignment(align, scheme, NeedlemanWunsch());
    
    /*
    //get the lengths of the sequences
    int len1 = length(ref_seq);
    int len2 = length(read_seq);

    //create the DP Matrix/Table
    float **matrix = new float * [len1 + 1];
    for (int i = 0; i <= len1; i++)
    {
        matrix[i] = new float[len2];
    }
    
    //initialize (will change by implementation)
    for (int i = 0; i <= len1; i++)
    {
        matrix[i][0] = i * scoreGap(scheme);
    }

    for (int j = 0; j <= len2; j++)
    {
        matrix[0][j] = j * scoreGap(scheme);
    }

    //our values for storing each potential movement
    float diagonal, vertical, horizontal;
    int i, j;
    for (i = 1; i <= len1; i++)
    {
        for (j = 1; j <= len2; j++)
        {
            diagonal = matrix[i-1][j-1];
            vertical = matrix[i][j-1] + scoreGap(scheme);
            horizontal = matrix[i-1][j] + scoreGap(scheme);

            //not sure if elements in a Dna5String can be compared like this
            if (ref_seq[i] == read_seq[j])
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

    score = matrix[i][j];

    //traceback
    i--;
    j--;

    typedef Align<Dna5String> TAlign;
    typedef Row<TAlign>::Type TRow;

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

        if (i > 0 && pos == match)
        {
            i--;
            j--;
        }

        else if (i > 0 && pos == mismatch)
        {
            i--;
            j--;
        }

        else if (i > 0 && pos == vgap)
        {
            insertGap(row1, i);
            j--;
        }

        else if (i > 0 && pos == hgap)
        {
            insertGap(row2, j);
            i--;
        }
    }

    //delete the dp matrix
    for (i = 0; i <= len1; i++)
    {
        delete matrix[i];
    }
    delete[] matrix;
    */
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
    Dna5String ref_seq;
    CharString read_id;
    CharString read_qual;
    Dna5String read_seq;

    Align<Dna5String> align;

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
        
        cout << align << "\n"; 

    }

    return 0;
}



