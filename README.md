# SureMap
SureMap is a software package for mapping DNA sequences against a reference genome, such as the human genome. It consists of three algorithms: SureMap-short, SureMap-hard and SureMap-long. SureMap-short is designed for short sequence reads up to 200bp. SureMap-hard is designed for 'hard to map reads', the reads which has a high error rate due to the sequencer machine error, or the reads that come from the regions in the genome which has a high variation rate against reference genome. This algorithm is suitable for the reads that other aligners can not align them (aligners like bwa-MEM or Bowtie or ...).
SureMap-long is designed for long and very long reads (possibly very noisy reads) like pacBio or Nanopore.
For all the algorithms, SureMap first needs to construct the FM-index for the reference genome. SureMap-Indexer is made for doing this job.

### Prerequisites

Download https://github.com/y-256/libdivsufsort and install. Make sure to compile libdivsufsort to static 64-bit libraries, i.e. setoptions in the main CMakeLists.txt to:
        option(BUILD_SHARED_LIBS "Set to OFF to build static libraries" OFF) 
        option(BUILD_DIVSUFSORT64 "Build libdivsufsort64" ON)


SureMap is using the suffix array algorithm implemented by Juha Karkkainen and Dominik Kempa, based on the paper: Juha Kärkkäinen, Dominik Kempa, Simon J. Puglisi. 'Parallel External Memory Suffix Sorting'. In Proc. 26th Symposium on Combinatorial Pattern Matching (CPM 2015), Springer, 2015, pp. 329-342. 

## install

* To make 'SureMap-Indexer' only:

        make indexer

* To make 'SureMap-Aligner' only:

        make aligner
        
* To compile all tools:

        make

## Example Usage
this command runs SureMap-short algorithm in fast mode with 10 threads to align short reads with error rate less than 0.1 of their length against reference genome.

        SureMap-Aligner.o short -t 10 -e 0.1 -m fast -o output.sam /path/to/index/directory/hg19.fa illumina.fastq

this command, runs SureMap-hard algorithm in sensitive mode to align short and hard-to-map reads against reference genome. maximum number of reports per read is 10.

        SureMap-Aligner.o hard -v 4 -k 10 -g -m sensitive -o output.sam /path/to/index/directory/hg19.fa bwa-unmapped.fastq
        
this command, runs SureMap-long algorithm to align long reads against reference genome. 'L' denotes fragment lengths in our approximate global alignment algorithm.

        SureMap-Aligner.o long -t 30 -L 500 -o output.sam /path/to/index/directory/hg19.fa pacbio.fastq

this command build index files for the reference genome 'ref.fa' file

        SureMap-Indexer.o ref.fa <outputFolder>
        
        
for additional options, please run SureMap-Aligner.o without any parameter to see the full options.

## Availabality
SureMap is released under MIT/X11 license. The latest source code is freely available at https://github.com/MohammadJRS/SureMap/.
precompiled binary for x86_64-linux is available in binary folder.

## Authors

* **MohammadJavad Rezaei Seraji** 
* **Seyed Abolfazl Motahari**

Sharif University Of Technology

