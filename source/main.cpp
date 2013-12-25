#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "global_alignment.hpp"
#include "local_alignment.hpp"
#include "local_global_alignment.hpp"
#include "align_lib.hpp"

using namespace std;
using namespace seqan;


int main(int argc, char **argv)
{
    if (argc != 4)
    {
        cout << "Error: Invalid number of arguments\n";
        return 1;
    }

    fstream in_ref(argv[1], ios::binary | ios::in);
    fstream in_reads(argv[2], ios::binary | ios::in);
    RecordReader<fstream, SinglePass<> > reader_ref(in_ref);
    RecordReader<fstream, SinglePass<> > reader_reads(in_reads);
    BamStream bamStreamIn(argv[3]);
    BamAlignmentRecord bam_record;

    CharString ref_id;
    TSequence ref_seq;
    CharString read_id;
    CharString read_qual;
    TSequence read_seq;
    TAlign corr_align;
    TAlign test_align;
    Score<int, Simple> scoringScheme(1, -1, -1);

    if (!atEnd(reader_ref) && readRecord(ref_id, ref_seq, reader_ref, Fasta()) != 0)
    {
        cout << "Error: Couldn't read reference sequence\n";
        return 1;
    }

    cout << ref_id << "\t" << ref_seq << "\n";

    while (!atEnd(reader_reads) && !atEnd(bamStreamIn))
    {
        if (readRecord(bam_record, bamStreamIn) != 0)
        {
            cout << "Error: Couldn't read next entry in sam file\n";
            return 1;
        }

        if (readRecord(read_id, read_seq, read_qual, reader_reads, Fastq()) != 0)
        {
            cout << "Error: Couldn't read the sequence from the file\n";
            return 1;
        }

        bamRecordToAlignment(corr_align, ref_seq, bam_record);
        //int score = GlobalAlignment::needle(test_align, ref_seq, read_seq, scoringScheme);
        //int score = LocalAlignment::smith_waterman(test_align, ref_seq, read_seq, scoringScheme);
        int score = LocalGlobalAlignment::lga(test_align, ref_seq, read_seq, scoringScheme);
        //int score = UncertainLocalGlobalAlignment::ulga(test_align, ref_seq, read_seq, read_qual, scoringScheme);
/*
        //START TEST BUILT IN CHUNK
        resize( rows(test_align), 2 );
        for (int i = 0; i < length(read_seq); i++)
        {
            read_seq[i] = toupper(read_seq[i]);
        }
        assignSource( row( test_align, 0 ), ref_seq);
        assignSource( row( test_align, 1 ), read_seq);
        //int score = globalAlignment(test_align, scoringScheme, AlignConfig<true, false, false, true>());
        int score = globalAlignment(test_align, scoringScheme, NeedlemanWunsch());
        //END
*/

        float pm = percent_match(corr_align, test_align);

        cout << "******************** Read ID: " << bam_record.qName << " ***************\n";
        cout << "Correct Alignment:\n" << corr_align << endl;
        cout << "Test Alignment: " << score << "\n" << test_align << endl;
        cout << "\tPercent Match = " << pm << endl;
        cout << endl;
    }

    return 0;
}
