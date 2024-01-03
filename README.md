# kISS: k-ordered Induced Suffix Sorting

k-ordered Induced Suffix Sorting (kISS) is an algorithm tailored to efficiently construct and leverage k-ordered FM-indexes, particularly for high-speed read mapping in extensive genome datasets. The innovation behind k-ordered FM-indexes aims to enhance the efficiency and speed of pattern matching tasks in genomics. kISS achieves this by integrating state-of-the-art induced-sorting-based suffix array construction, capitalizing on its superior time and speed complexities.

## Requirement
To use kISS, ensure the following components are available:
- GNU g++ version 11 or later
- cmake version 3.18.0 or higher
- boost version 1.74.0 or higher
- [Biovoltron](https://github.com/JHHLAB/Biovoltron)

## Installation Steps
Begin by installing necessary dependencies:
```
sudo apt install g++-11 cmake make libboost-all-dev
```

Clone the kISS repository and navigate to the build directory, and generate build files using CMake and compile the project:
```
git clone --recurse-submodules https://github.com/jhhung/kISS
cd kISS
mdkir build
cd build
cmake ..
make -j
```

## Usage Examples
### Help message
To view available options and commands, execute:
```
$ ./kISS -h
 _     ___  ____  ____
| | __|_ _|/ ___|/ ___|
| |/ / | | \___ \\___ \
|   <  | |  ___) |___) |
|_|\_\|___||____/|____/ 0.0.1-alpha

kISS [--generic-option ...] cmd [--cmd-specific-option ...]:

Generic options:
  -h [ --help ]                    produce help message
  -v [ --version ]                 print version string
  -t [ --num_threads ] NUM (=128)  number of thread
  --verbose                        print more information

kISS suffix_sort [--option ...] <FASTA/text filename>
	For FASTA files, the filename should end with .fasta/.fa;
	otherwise, it's treated as text.

Options:
  -k NUM (=256)                    a k-ordered value, where each suffix is
                                   sorted based on the first k characters.
                                   Using -1 indicates unbounded sorting.

kISS fmindex_build [--option ...] <FASTA/text filename>
	For FASTA files, the filename should end with .fasta/.fa;
	otherwise, it's treated as text.

Options:
  --sa_sample_rate NUM (=4)        Controls the frequency of data point
                                   sampling during suffix array construction,
                                   balancing between accuracy and resource
                                   efficiency.
  -k NUM (=256)                    Sets the maximum query length for the
                                   fmindex search as k - sa_sample_rate; using
                                   -1 indicates an unlimited query length.

kISS fmindex_query [--option ...] <FASTA/text filename>
	For FASTA files, the filename should end with .fasta/.fa;
	otherwise, it's treated as text.

Options:
  -q STR                           Content of the query string.
  -f STR                           File containing multiple query string
                                   contents, with each query string separated
                                   by a newline.
```
### k-ordered suffix array sorter
Sorts suffixes within a file using 256-ordered suffix sorting under 24 threads:
```
./kISS suffix_sort -k 256 -t 24 data/example.fa
```

### k-ordered fmindex
Constructs a k-ordered FMindex for example.fasta with specific parameters:
```
./kISS fmindex_build \
    genome.fasta \
    --sa_sample_rate 4 \
    -k 256 \
    -t 128

# Query examples
./kISS fmindex_query example.fa -q GCTAGCTCTAG
./kISS fmindex_query example.fa -f example.query
```

## Citation
Zheng-Dao Yang, Hsuan-Yu Kuo, Po-Wei Hsieh and Jui-Hung Hung, Efficient Construction and Utilization of k-ordered FM-indexes with kISS for Ultra Fast Read Mapping in Large Genomes.

## Contact
For any inquiries or additional information, contact Jui-Hung Hung via email jhh@cs.nycu.edu.tw
