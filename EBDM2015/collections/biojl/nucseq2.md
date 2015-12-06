---
name: nucseq2
slide_author: Ben J. Ward
---
### Counting Kmers
```julia
@time DNAKmerCounts{10}(chr1);
```
1.901 seconds (8 allocations: 4096 KB)
Versus 2.724 seconds in Biostrings (R) and 5.954 seconds in SeqAn (C++)

### Counting Nucleotides
```julia
@time NucleotideCounts(chr1)
```
67.946 milliseconds (196 k allocations: 3058 KB)
Versus 1193 milliseconds in Biostrings (R) 3165 milliseconds in SeqAn (C++)

```julia
count_c(x::Uint64) = count_ones((((~x) >>> 1) & x) & 0x5555555555555555)
```
