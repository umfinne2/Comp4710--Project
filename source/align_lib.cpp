#include <seqan/align.h>
#include <map>
#include "align_lib.hpp"
#include "EDNAFULL.hpp"

using namespace seqan;
using namespace std;

typedef Iterator<TRow>::Type TRowIterator;

int AlignLib::get_index(char x)
{
    int result = -1;

    switch (x)
    {
        case 'A':
            result = 0;
            break;
        case 'T':
            result = 1;
            break;
        case 'G':
            result = 2;
            break;
        case 'C':
            result = 3;
            break;
        case 'S':
            result = 4;
            break;
        case 'W':
            result = 5;
            break;
        case 'R':
            result = 6;
            break;
        case 'Y':
            result = 7;
            break;
        case 'K':
            result = 8;
            break;
        case 'M':
            result = 9;
            break;
        case 'B':
            result = 10;
            break;
        case 'V':
            result = 11;
            break;
        case 'H':
            result = 12;
            break;
        case 'D':
            result = 13;
            break;
        case 'N':
            result = 14;
            break;
        default:
            result = -1;
    }
        /*
        dnamap['T'] = 1;
        dnamap['G'] = 2;
        dnamap['C'] = 3;
        dnamap['S'] = 4;
        dnamap['W'] = 5;
        dnamap['R'] = 6;
        dnamap['Y'] = 7;
        dnamap['K'] = 8;
        dnamap['M'] = 9;
        dnamap['B'] = 10;
        dnamap['V'] = 11;
        dnamap['H'] = 12;
        dnamap['D'] = 13;
        dnamap['N'] = 14;
        */
    return result;
}

int AlignLib::get_score(char a, char b)
{
    if (a == gapValue<char>() || b == gapValue<char>())
    {
        return gapcost;
    }

    else
    {
        int i = get_index(toupper(a));
        int j = get_index(toupper(b));

        return EDNAFULL_matrix[i][j];
    }
}

float AlignLib::percent_match(TAlign &correct, TAlign &test)
{
    float result = 0.0;
    TRow &correct_ref = row(correct, 0);
    TRow &correct_read = row(correct, 1);
    TRow &test_ref = row(test, 0);
    TRow &test_read = row(test, 1);

    //If I understand how the alignments are suppose to work the length
    //of the original read sequence should be the same for both.
    if (length(source(correct_read)) != length(source(test_read)))
    {
        //cout << source(correct_read) << endl;
        //cout << source(test_read) << endl;
        return -1.0;
    }

    //iterate over read from correct alignment
    //if we have a match
    float total_corr_matches = 0.0;
    float valid_test_matches = 0.0;
    for (unsigned i = 0; i < length(source(correct_read)); ++i)
    {
        //*********** MORE LENIENT APPROACH **************
        unsigned corr_read_view_pos = toViewPosition(correct_read, i);

        //if the correct alignment has a match, check the test alignment
        if ((length(correct_ref) > corr_read_view_pos) &&
             (toupper(correct_ref[corr_read_view_pos]) == toupper(correct_read[corr_read_view_pos])))
        {
            total_corr_matches += 1.0;

            //get the correct reference source position
            unsigned corr_ref_src_pos = toSourcePosition(correct_ref, corr_read_view_pos);

            //get the test reference view position which should in the ideal case should be the same as the corr_read_view_pos
            unsigned test_ref_view_pos = toViewPosition(test_ref, corr_ref_src_pos);

            //now if the reference view position exists and the test_read can be indexed and the value match then
            //it is probably mapped properly
            if ((length(test_read) > test_ref_view_pos) &&
                 (toupper(test_ref[test_ref_view_pos]) == toupper(test_read[test_ref_view_pos])))
            {
                valid_test_matches += 1.0;
                cout << test_ref[test_ref_view_pos] << ":" << test_read[test_ref_view_pos] << "  ";
            }
        }

        //   ************* TO STRICT !!! ****************
        //for the same source position, if the view positions match then our
        //bp has been mapped correctly.
        //if (toViewPosition(correct_read, i) == toViewPosition(test_read, i))
        //{
        //    matches += 1.0;
        //}
    }
    cout << endl;

    if (total_corr_matches != 0.0)
    {
        result = valid_test_matches / total_corr_matches;
    }

    return result;
    //return (float)(matches / (float)(length(source(correct_read))));
}

void AlignLib::print_matrix(float **matrix, TSequence seq1, TSequence seq2, int rows, int cols, int offset, int len)
{
    int i_start, i_end, j_start, j_end;

    //if offset is neg. subtract that from the end pos.
    if (offset < 0)
    {
        i_start = cols + offset;
        j_start = rows + offset;
    }

    else
    {
        i_start = offset + 1;
        j_start = offset + 1;
    }

    if (i_start + len >= cols)
    {
        i_end = cols - 1;
    }

    if (j_start + len >= rows)
    {
        j_end = rows - 1;
    }

    cout << setw(5) << " ";
    for (int i = i_start; i < i_end; i++)
    {
        cout << setw(5) << seq1[i-1];
    }

    cout << endl;

    for (int j = j_start; j < j_end; j++)
    {
        if (j > 0)
        {
            cout << seq2[j-1];
        }

        for (int i = i_start; i < i_end; i++)
        {
            cout << setw(5) << matrix[i][j];
        }
        cout << endl;
    }
}
