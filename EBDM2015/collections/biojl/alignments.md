---
name: alignments
slide_author: Ben J. Ward
---
# Sequence alignments

Use an anchor based system.

```julia
immutable AlignmentAnchor
    seqpos::Int
    refpos::Int
    op::Operation
end

immutable Alignment
    anchors::Vector{AlignmentAnchor}
    firstref::Int
    lastref::Int
end

immutable AlignedSequence{S}
    seq::S
    aln::Alignment
end
```

Operations supported covers all CIGAR ops, encoded as single bytes.
