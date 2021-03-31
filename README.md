# HSX


C++ implementation of the lastz fasta indexing [HSX format][1].

C++ 17 (for now).

To be merged into https://github.com/urbanslug/mashz

### CPP

Compile with `make`.

#### Usage
hsx expects the fasta file and index in the same diretory.

Reads paths to fasta files from standard input

```
./hsx-cpp hsxexA.fa hsxexB.fa hsxexC.fa
```

Writes to `out.hsx`.


### Python

#### Usage

```
py/build_fasta_hsx.py hsxexA.fa hsxexB.fa hsxexC.fa > abc.hsx
```

[1]: http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/hsx_format.html
