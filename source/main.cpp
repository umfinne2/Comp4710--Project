#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "needle.hpp"
#include "local_alignment.hpp"
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
    Score<int, Simple> scoringScheme(0, -1, -1);

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
        //int score = needle(test_align, ref_seq, read_seq, scoringScheme);

        resize( rows(test_align), 2 );

        for (int i = 0; i < length(read_seq); i++)
        {
            read_seq[i] = toupper(read_seq[i]);
        }

        assignSource( row( test_align, 0 ), ref_seq);
        assignSource( row( test_align, 1 ), read_seq);
        int score = globalAlignment(test_align, scoringScheme, AlignConfig<true, false, false, true>());

        float pm = percent_match(corr_align, test_align);

        cout << "******************** Read ID: " << bam_record.qName << " ***************\n";
        cout << "Correct Alignment:\n" << corr_align << endl;
        cout << "Test Alignment: " << score << "\n" << test_align << endl;
        cout << "\tPercent Match = " << pm << endl;
        cout << endl;
    }

    return 0;
}

/*  We will need to integrate this main with the other for actually running experiments
    for now the other is simpler to use.
*/
/*
int main()
{
    // Open file and create RecordReader.
    std::fstream in("../data/test/hiv-part.fa", std::ios::binary | std::ios::in);
    seqan::RecordReader<std::fstream, seqan::SinglePass<> > reader(in);

    seqan::CharString id;
    TSequence seq;

    if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
            return 1;  // Could not record from file.

    std::cout << seq << "\n";

    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn("../data/test/SOLiD/hiv-part-sim.sam");
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        readRecord(record, bamStreamIn);
        Align<TSequence> align;

        bamRecordToAlignment(align, seq, record);

        std::cout << record.qName << "\n";
        std::cout << align << "\n";
    }

    return 0;
}
*/
