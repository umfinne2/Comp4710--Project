Comp4710--Project
=================

A comparison of a LGA (Local Global Alignment) & ULGA (Uncertain Local Global Alignment) pairwise sequence alignment algorithms to traditional Needleman-Wunsch & Smith-Waterman algs. for comparative genome assembly

### [Seqan](http://www.seqan.de/) ###
Provides a set of easy to use C++ libraries for reading and writing fastq & sam files.

### [ART](http://www.niehs.nih.gov/research/resources/software/biostatistics/art/) ###
Is a read simulator we will be using to benchmark the different algorithms.


### Compiling ###
From the project's root directory

```
make
```

### How to run ###
From the project root directory
```
#replace <algorithm> with one of ga, la, lga, ulga depending upon which
algorithm you want to run the test with

./a.out <algorithm> data/experimental/hiv-ref.fa data/experimental/454/hiv-sim.fq data/experimental/454/hiv-sim.sam

#Global Alignment (~2.5 mins)
./a.out ga data/experimental/hiv-ref.fa data/experimental/454/hiv-sim.fq data/experimental/454/hiv-sim.sam

#Local Alignment (~2.5 mins)
./a.out la data/experimental/hiv-ref.fa data/experimental/454/hiv-sim.fq data/experimental/454/hiv-sim.sam

#Local-Global Alignment (~2.5 mins)
./a.out lga data/experimental/hiv-ref.fa data/experimental/454/hiv-sim.fq data/experimental/454/hiv-sim.sam

#Uncertain Local-Global Alignment (~5 mins)
./a.out ulga data/experimental/hiv-ref.fa data/experimental/454/hiv-sim.fq data/experimental/454/hiv-sim.sam
```

### Running Tests ###
```
./runTests.sh
```

### Results ###
The results for each of the algorithm are located in
data/experimental/raw_results
