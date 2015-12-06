---
name: nucseq3
slide_author: Ben J. Ward
---
### Translation
```julia
const seqs = collect(read("/Users/dcjones/human-coding-sequences.fa", FASTA))
@time for seq in seqs
    translate(convert(RNASequence, seq.seq[1:end-3]), Seq.standard_genetic_code, true)
end
```
467.524 milliseconds (765 k allocations: 65283 KB)
632 milliseconds in Biostrings (R) and 938 milliseconds in SeqAn (C++)

### Subsequences are cheap
```julia
@time [chr1[i:end] for i in 1:1000000];
```
48.829 milliseconds (1000 k allocations: 54688 KB)

Immutable by convention.
