---
name: nucseq
slide_author: Ben J. Ward
---
# Nucleotide Sequences
### 2 bit encoding:
```
A => 0b00
C => 0b01
T => 0b10
G => 0b11
```

### Reverse complement
```julia
using Bio.Seq
const chr1 = first(read("/Users/dcjones/Homo_sapiens.GRCh38.dna.primary_assembly.fa", FASTA)).seq
@time reverse_complement(chr1);
```
151.301 milliseconds (11 allocations: 91171 KB)
Versus 956 milliseconds in Biostrings (R) and 445 milliseconds in SeqAn (C++)
