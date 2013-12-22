#include "align_lib.hpp"

using namespace seqan;
using namespace std;

typedef Iterator<TRow>::Type TRowIterator;

float percent_match(TAlign &correct, TAlign &test)
{
    //TRow &cr1 = row(correct, 0);
    TRow &correct_read = row(correct, 1);
    //TRow &tr1 = row(test, 0);
    TRow &test_read = row(test, 1);

    //If I understand how the alignments are suppose to work the length
    //of the original read sequence should be the same for both.
    if (length(source(correct_read)) == length(source(test_read))) { return (float) -1;}

    //iterate over read from correct alignment
    //if we have a match
    float matches = 0.0;
    for (unsigned i = 0; i < length(source(correct_read)); ++i)
    {
        //for the same source position, if the view positions match then our
        //bp has been mapped correctly.
        if (toViewPosition(correct_read, i) == toViewPosition(test_read, i))
        {
            matches += 1.0;
        }
    }

    return (float)(matches / (float)(length(source(correct_read))));
}
