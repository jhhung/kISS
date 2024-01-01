# kISS

k-ordered Induced Suffix Sorting (kISS) is an innovative algorithm designed for the efficient construction and utilization of k-ordered FM-indexes, particularly for ultra-fast read mapping in large genomes. The k-ordered concept for FM-indexes aims to improve the speed and efficiency of pattern matching tasks. kISS achieves this by incorporating state-of-the-art induced-sorting-based suffix array construction, leveraging its superior time and speed complexities.

## Requirement
- GNU g++ >= 11
- cmake >= 3.18.0
- boost >= 1.74.0
- [Biovoltron](https://github.com/JHHLAB/Biovoltron)

## Installation
```
sudo apt install g++-11 cmake make libboost-all-dev
git clone --recurse-submodules https://github.com/jhhung/kISS

cd kISS
mdkir build
cd build
cmake ..
make -j
```

## Example
### k-ordered suffix array sorter
```
./kISS suffix_sort -k 256 -t 24 data/example.fa
```

### k-ordered fmindex
```
./kISS fmindex_build \
    genome.fasta \
    --sa_sample_rate 4 \
    -k 256 \
    -t 128

./kISS fmindex_query data/example.fa -q "GCTAGCTCTAG"
./kISS fmindex_query data/example.fa -f data/example.query
```

## Citation
Zheng-Dao Yang, Hsuan-Yu Kuo, Po-Wei Hsieh and Jui-Hung Hung, Efficient Construction and Utilization of k-ordered FM-indexes with kISS for Ultra Fast Read Mapping in Large Genomes.

## Contact
Jui-Hung Hung, jhhung at jhh[at]cs.nycu.edu.tw
