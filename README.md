# HSX

C++ implementation of the lastz fasta indexing [HSX format][1].

To be merged into https://github.com/urbanslug/mashz

```
g++ -g -Wall hsx.cpp
```

Reads paths to fasta files from standard input
```
echo hsxexA.fa hsxexB.fa hsxexC.fa | ./a.out
```

Python

```
py/build_fasta_hsx.py hsxexA.fa hsxexB.fa hsxexC.fa > abc.hsx
```

[1]: http://www.bx.psu.edu/miller_lab/dist/README.lastz-1.02.00/hsx_format.html
