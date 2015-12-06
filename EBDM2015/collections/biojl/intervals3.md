---
name: intervals3
slide_author: Ben J. Ward
---
### Coverage
```julia
@time coverage(read("bodymap-heart-2-chr1.bed", BED))

IntervalCollection with 3208603 intervals:
  chr1:10586-10633    .    1
  ⋮
```
34.209 seconds      (419 M allocations: 16216 MB, 31.36% gc time)
`bedtools genomecov` takes 54.4 seconds to do the same.

.center[![coveragesvg](img/coverage.svg)]
