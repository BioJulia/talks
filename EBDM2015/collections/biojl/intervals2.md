---
name: intervals2
slide_author: Ben J. Ward
---
### Intersection
```julia
genes = read("genes.bed", BED)
reads = IntervalCollection(read("reads.bed", BED))
@time for (a,b) in intersect(reads, genes)
end
```
3.190 seconds (17592 k allocations: 612 MB, 59.76% gc time)
5.9 seconds for `bedtools intersect -loj`.

.center[![intervals image](img/interval-intersection.svg)]
