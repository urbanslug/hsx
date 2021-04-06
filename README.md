# HSX
C++ implementation of the lastz fasta indexing [HSX format][1].

To be merged into https://github.com/urbanslug/mashz

### CPP

Compile with `make`.

#### Usage


Reads paths to fasta files from standard input

```
./hsx-cpp hsxexA.fa hsxexB.fa hsxexC.fa > abc.hsx
```

Copy/move the hsx files to be in the same dir as the fasta files.
(Currently required but will be fixed).


### Python

#### Usage

```
py/build_fasta_hsx.py hsxexA.fa hsxexB.fa hsxexC.fa > abc.hsx
```

[1]: http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/hsx_format.html
