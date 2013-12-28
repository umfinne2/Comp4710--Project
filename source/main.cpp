#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

#include "global_alignment.hpp"
#include "local_alignment.hpp"
#include "local_global_alignment.hpp"
#include "uncertain_local_global_alignment.hpp"
#include "align_lib.hpp"

using namespace std;
using namespace seqan;

void usage(void)
{
    cout << "To run the program you need to pass the following arguments in this order:" << endl;
    cout << "\tAlgorithm - the algorithm you would like to test [ options: ga | la | lga | ulga ]" << endl;
    cout << "\tReference - the path to the fasta (*.fa) file containing the reference sequence used to generate the reads." << endl;
    cout << "\tReads - the path to the fastq (*.fq) file containing the reads generated." << endl;
    cout << "\tSAM - the path to the sam file containing the correct mapping of the reads to the reference sequence." << endl;
    //cout << "\tOutput - the type of output to print [ options: debug | csv ]" << endl;
    exit(1);
}

int main(int argc, char **argv)
{
    if (argc != 5)
    {
        cout << "Error: Invalid number of arguments\n";
        usage();
    }

    fstream in_ref(argv[2], ios::binary | ios::in);
    fstream in_reads(argv[3], ios::binary | ios::in);
    RecordReader<fstream, SinglePass<> > reader_ref(in_ref);
    RecordReader<fstream, SinglePass<> > reader_reads(in_reads);
    BamStream bamStreamIn(argv[4]);
    BamAlignmentRecord bam_record;

    CharString ref_id;
    TSequence ref_seq;
    CharString read_id;
    CharString read_qual;
    TSequence read_seq;
    TAlign corr_align;
    TAlign test_align;

    if (!atEnd(reader_ref) && readRecord(ref_id, ref_seq, reader_ref, Fasta()) != 0)
    {
        cout << "Error: Couldn't read reference sequence\n";
        return 1;
    }

    //cout << ref_id << "\t" << ref_seq << endl;
    cout << "Read_ID\tPercent_Match\tMapped\n";

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

        if (strcmp(argv[1], "ga") == 0)
        {
            GlobalAlignment::needle(test_align, ref_seq, read_seq);
        }

        else if (strcmp(argv[1],"la") == 0)
        {
            LocalAlignment::smith_waterman(test_align, ref_seq, read_seq);
        }

        else if (strcmp(argv[1], "lga") == 0)
        {
            LocalGlobalAlignment::lga(test_align, ref_seq, read_seq);
        }

        else if (strcmp(argv[1], "ulga") == 0)
        {
            UncertainLocalGlobalAlignment::ulga(test_align, ref_seq, read_seq, read_qual);
        }

        else
        {
            cout << "Error: Invalid algorithm type must be ga | la | lga | ulga" << endl;
            usage();
        }

        float pm = AlignLib::percent_match(corr_align, test_align);

        //cout << "******************** Read ID: " << bam_record.qName << " ***************\n";
        //cout << "Correct Alignment:\n" << corr_align << endl;
        //cout << "Test Alignment: " << score << "\n" << test_align << endl;
        cout << bam_record.qName << "\t" << pm << "\t" << ceil(pm) << endl;
        //cout << endl;
    }

    return 0;
}
