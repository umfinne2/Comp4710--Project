#include "align_lib.hpp"

using namespace seqan;
using namespace std;

typedef Iterator<TRow>::Type TRowIterator;

float percent_match(TAlign &correct, TAlign &test)
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
