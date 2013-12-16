#include <fstream>
#include <iostream>

#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/bam_io.h>

using namespace std;
using namespace seqan;

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

