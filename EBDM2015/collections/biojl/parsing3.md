---
name: parsing3
slide_author: Ben J. Ward
---
# Enter State machine specifications
With Ragel we can define a FASTA parser in 9 lines that's faster and more accurate that most hand-written C/C++ parsers.
```
newline     = '\r'? '\n';
hspace      = [ \t\v];
whitespace  = space | newline;

identifier  = (any - space)+;
description = ((any - hspace) [^\r\n]*);
letters     = (any - space - '>')+;
sequence    = whitespace* letters? (whitespace+ letters)*;
fasta_entry = '>' identifier (hspace+ description)? newline sequence whitespace*;

main := whitespace* (fasta_entry)*;
```
The same specification can be used to generate equivalent parsers C, D, Go, Java, Ruby, C#, OCaml.
