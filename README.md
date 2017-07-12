# SureMap
SureMap is a software package for mapping DNA sequences against a reference genome, such as the human genome. It consists of three algorithms: SureMap-short, SureMap-hard and SureMap-long. SureMap-short is designed for short sequence reads up to 200bp. SureMap-hard is designed for 'hard to map read', the reads which has a high error rate due to the sequencer machine error, or the reads come from the regions in the genome which has a high variation rate against reference genome. This algorithm is suitable for the reads that other aligners can not align them (aligners like bwa-MEM or Bowtie or ...) 
SureMap-long is designed for long and very long reads (possibly very noisy reads) like pacBio or Nanopore.
For all the algorithms, SureMap first needs to construct the FM-index for the reference genome. SureMap-Indexer is doing this job?

### Prerequisites

Download https://github.com/y-256/libdivsufsort and install. Make sure to compile libdivsufsort to static 64-bit libraries, i.e. set options in the main CMakeLists.txt to option(BUILD_SHARED_LIBS "Set to OFF to build static libraries" OFF) option(BUILD_DIVSUFSORT64 "Build libdivsufsort64" ON)


SureMap is using the suffix array algorithm implemented by Juha Karkkainen and Dominik Kempa, based on the paper: Engineering a Lightweight External Memory Suffix Array Construction Algorithm. In Proc. 2nd International Conference on Algorithms for Big Data (ICABD) 2014.

## install

* To make 'SureMap-Indexer' only:

        make indexer

* To make 'SureMap-Aligner' only:

        make aligner
        
* To compile all tools:

        make

## Availabality
SureMap is released under MIT/X11 license. The latest source code is freely available at https://github.com/MohammadJRS/SureMap/.
precompiled binary for x86_64-linux is available in binary folder.

## Authors

* **MohammadJavad Rezaei Seraji** 
* **Seyed Abolfazl Motahari**

Sharif University Of Technology

