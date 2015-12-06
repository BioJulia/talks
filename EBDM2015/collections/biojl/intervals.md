---
name: intervals
slide_author: Ben J. Ward
---
# Intervals

```julia
immutable Interval{T}
    seqname::ASCIIString
    first::Int64
    last::Int64
    strand::Strand
    metadata::T
end
```
B+-trees are provided by IntervalTrees:

.center[![B+-tree image](img/btree.svg)]

Intersection is a fundamental operation for sets of intervals.
Find all pairs of intervals from A and B that intersect.
