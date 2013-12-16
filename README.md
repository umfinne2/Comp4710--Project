Comp4710--Project
=================

A comparison of a LGA (Local Global Alignment) & ULGA (Uncertain Local Global Alignment) pairwise sequence alignment algorithms to traditional Needleman-Wunsch & Smith-Waterman algs. for comparative genome assembly

### [Seqan](http://www.seqan.de/) ###
Provides a set of easy to use C++ libraries for reading and writing fastq & sam files.

### [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) ###
Is a read simulator we will be using to benchmark the different algorithms.


### Commands ###
All these commands should be executed from the root probject directory.
For the programs in the scripts folder just pipe the output to a file as needed.
```
#Multiplies the default SOLiD profile error values by x (2 in this case)
./scripts/error_profile_multiplier.pl bin/art/SOLiD_profile/profile_default 2

#Uses art_SOLiD to generate fastq reads and sam alignments from hiv-part.fa
#and stores the files in the test/SOLiD folder with the prefix hiv-part-sim.
#  - 34 represents the avg read length
#  - 2 represents the coverage
#  - the '-p data/test/profile_test' uses our error profile file
./bin/art/art_SOLiD -s -p data/test/profile_test data/test/hiv-part.fa data/test/SOLiD/hiv-part-sim 34 2

#Because the fastq reads produced write the sequence with the number 0-3 versus A,C,G,T, this will fix them.
./scripts/fastq_solid_converter.pl data/test/SOLiD/hiv-part-sim.fq

#For now compile .cpp files with -I ../lib [file]
g++ -I ../lib source/needle.cpp
```

### TODO ###
*   Create sw, lga and ulga alignment functions
*   Create a separate main which takes ref.fa, \*.fq and *.sam files and prints avg running time and avg % accuracy.
*   Set up a proper Makefile

