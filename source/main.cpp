#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "needle.hpp"
#include "local_alignment.hpp"

using namespace std;
using namespace seqan;

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
    seqan::Dna5String seq;

    if (readRecord(id, seq, reader, seqan::Fasta()) != 0)
            return 1;  // Could not record from file.

    std::cout << seq << "\n";

    // Open input stream, BamStream can read SAM and BAM files.
    seqan::BamStream bamStreamIn("../data/test/SOLiD/hiv-part-sim.sam");
    seqan::BamAlignmentRecord record;
    while (!atEnd(bamStreamIn))
    {
        readRecord(record, bamStreamIn);
        Align<Dna5String> align;

        bamRecordToAlignment(align, seq, record);

        std::cout << align << "\n";
    }

    return 0;
}
*/
