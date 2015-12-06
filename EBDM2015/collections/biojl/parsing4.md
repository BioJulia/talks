---
name: parsing4
slide_author: Ben J. Ward
---
# Parsers are fast
```julia
using Bio.Seq

@time collect(read("/Users/dcjones/Homo_sapiens.GRCh38.dna.primary_assembly.fa", FASTA));
```
26.561 seconds (8951 allocations: 1423 MB, 4.74% gc time)
30.8 seconds in Biostrings (R), 39.7 seconds in SeqAn(C++)

Other cool stuff is possible too!
```
using Bio.Seq

# Memory mapping
read("/Users/dcjones/Homo_sapiens.GRCh38.dna.primary_assembly.fa", FASTA, memory_map=true)

# Reading from commands
for seqrec in take(drop(read(`quip -cd /Users/dcjones/SRR094166.fastq.qp`, FASTQ), 10000), 3)
    println(seqrec)
end
```
