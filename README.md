# kISS: k-ordered Induced Suffix Sorting

k-ordered Induced Suffix Sorting (kISS) is an algorithm tailored to efficiently construct and leverage k-ordered FM-indexes, particularly for high-speed read mapping in extensive genome datasets. The innovation behind k-ordered FM-indexes aims to enhance the efficiency and speed of pattern matching tasks in genomics. kISS achieves this by integrating state-of-the-art induced-sorting-based suffix array construction, capitalizing on its superior time and speed complexities.

## Requirement
To use kISS, ensure the following components are available:
- GNU g++ version 11 or later
- cmake version 3.18.0 or higher
- boost version 1.74.0 or higher

## Installation Steps
Begin by installing necessary dependencies:
```
sudo apt install g++-11 cmake make libboost-all-dev
```

Clone the kISS repository and navigate to the build directory, and generate build files using CMake and compile the project:
```
git clone --recurse-submodules https://github.com/jhhung/kISS
cd kISS
cmake -B build

cd build && cmake .. && make && cd -
```

## Usage Examples
### Help message
To view available options and commands, execute:
```
$ ./build/kISS -h
 _     ___  ____  ____
| | __|_ _|/ ___|/ ___|
| |/ / | | \___ \\___ \
|   <  | |  ___) |___) |
|_|\_\|___||____/|____/ 0.0.1-alpha

kISS [--generic-option ...] cmd [--cmd-specific-option ...]:

Generic options:
  -h [ --help ]                         produce help message
  -v [ --version ]                      print version string
  -g [ --generic ]                      (Under construction) Select this option
                                        if the input FASTA file contains bases
                                        other than ATCG.
                                        When turned on, some specific
                                        optimizations cannot be done and the
                                        performancemay be slightly worse.
  -t [ --num_threads ] NUM (=128)       number of thread
  --verbose                             print more information


./kISS suffix_sort [--option ...] <FASTA filename>

Options:
  -k [ --kordered ] NUM (=256)          a k-ordered value, where each suffix is
                                        sorted based on the first k characters.
                                        Using -1 indicates unbounded sorting.
  -s [ --sorting-algorithm ] ALGO (=PARALLEL_SORTING)
                                        The sorting strategy for the step
                                        "Parallel k-ordered Sorting of LMS
                                        Suffixes" in the kISS pipeline.
                                        Valid arguments are PARALLEL_SORTING
                                        and PREFIX_DOUBLING.


./kISS fmindex_build [--option ...] <FASTA filename>

Options:
  -k [ --kordered ] NUM (=256)          (Under construction) Sets the maximum
                                        query length for the fmindex search as
                                        k - sa_sample_rate; using -1 indicates
                                        an unlimited query length.

./kISS fmindex_query [--option ...] <FASTA filename>

Options:
  -q [ --query ] STR                    Content of the query string.
                                        (Under constrction) Currently, a
                                        maximum length of 32 is supported.
  -p [ --patterms ] STR                 patterms files generated by FMTree
  -n [ --headn ] NUM (=10)              Output the first n location
```
### k-ordered suffix array sorter
Sorts suffixes within a file using 256-ordered suffix sorting under 24 threads:
```
$ ./build/kISS suffix_sort example/drosophia_chr1_2.fa -k 256 -t 24
[2024-04-12 15:16:54.427] [stderr] [info] [suffix_sort.hpp:62] n = 48800648, k = 256, suffix sorting elapsed 0.48087947400000003
```
Download chm13v2.0.fa from remote and sorts suffixes within chm13v2.0 using 256-ordered suffix sorting under 24 threads with more information:
```
$ wget https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz -O example/chm13v2.0.fa.gz
$ gunzip example/chm13v2.0.fa.gz
$ ./build/kISS suffix_sort example/chm13v2.0.fa -k 256 -t 24 --verbose
[2024-04-12 15:21:56.770] [stderr] [debug] [kiss1_core.hpp:196] Preparation elapsed 0.037508534
[2024-04-12 15:22:01.361] [stderr] [debug] [kiss1_core.hpp:200] SA.resize(n + 1) elapsed 4.591083359000001
[2024-04-12 15:22:02.291] [stderr] [debug] [kiss1_core.hpp:204] get_lms elapsed 0.929645511
[2024-04-12 15:22:16.974] [stderr] [debug] [kiss1_core.hpp:210] lms_suffix_direct_sort elapsed 14.683336318
[2024-04-12 15:22:18.501] [stderr] [debug] [kiss1_core.hpp:215] put_lms_suffix elapsed 1.526712206
[2024-04-12 15:22:31.447] [stderr] [debug] [kiss1_core.hpp:219] induced_sort elapsed 12.94582755
[2024-04-12 15:22:31.447] [stderr] [info] [suffix_sort.hpp:62] n = 3117292070, k = 256, suffix sorting elapsed 34.714330689
```

### k-ordered fmindex with FMtree validation
Please visit https://github.com/JHHLAB/FMtree for verification of FMtree combined with kISS in the results.

### k-ordered fmindex
> [!note]
> This is a version implemented by ourselves. Currently, we have discovered performance issues and are investigating.
> For the faster FMIndex implementation, please visit https://github.com/JHHLAB/FMtree

Constructs a k-ordered FMindex for example/drosophia_chr1_2.fa with specific parameters:
#### Build FMIndex
```
$ ./build/kISS fmindex_build example/drosophia_chr1_2.fa -k 256 -t 24
```

#### Query examples
```
$ ./build/kISS fmindex_query example/drosophia_chr1_2.fa -q GCTAGCTCTAG -n 5
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:56] query = GCTAGCTCTAG found 6 times
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:59] The 1-st position is 11165052, content of substring is GCTAGCTCTAG
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:59] The 2-nd position is 7392748, content of substring is GCTAGCTCTAG
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:59] The 3-rd position is 27638949, content of substring is GCTAGCTCTAG
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:59] The 4-th position is 27179738, content of substring is GCTAGCTCTAG
[2024-04-12 15:29:04.084] [stderr] [info] [fmindex_query.hpp:59] The 5-th position is 31027614, content of substring is GCTAGCTCTAG

$ ./build/kISS fmindex_query example/drosophia_chr1_2.fa -q TGCTTAGCTAG -n 8
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:56] query = TGCTTAGCTAG found 11 times
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 1-st position is 21997768, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 2-nd position is 11508644, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 3-rd position is 11626997, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 4-th position is 43082665, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 5-th position is 31234229, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 6-th position is 25159681, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 7-th position is 46565569, content of substring is TGCTTAGCTAG
[2024-04-12 15:29:36.001] [stderr] [info] [fmindex_query.hpp:59] The 8-th position is 44749501, content of substring is TGCTTAGCTAG
```

## Benchmark
See [Benchmark Section](experiment/README.md#benchmark)

## Citation
Zheng-Dao Yang, Hsuan-Yu Kuo, Po-Wei Hsieh and Jui-Hung Hung, Efficient Construction and Utilization of k-ordered FM-indexes with kISS for Ultra Fast Read Mapping in Large Genomes.

## Contact
For any inquiries or additional information, contact Jui-Hung Hung via email jhh@cs.nycu.edu.tw
