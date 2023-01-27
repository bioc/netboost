# Netboost
Boosting supported network analysis for high-dimensional omics applications.
The corresponding publication has been published in IEEE TCBB:

Pascal Schlosser, Jochen Knaus, Maximilian Schmutz, Konstanze Döhner, Christoph Plass, Lars Bullinger, Rainer Claus, Harald Binder, Michael Lübbert and Martin Schumacher, Netboost: Boosting-supported network analysis improves high-dimensional omics prediction in acute myeloid leukemia and Huntington's disease, IEEE/ACM Trans Comput Biol Bioinform (2021): 10.1109/TCBB.2020.2983010.

This package comes bundled with the MC-UPGMA clustering package by Yaniv Loewenstein.

# Requirements
This package is only working on MacOS and Linux (no Windows atm).

Required for building are C/C++ compilers, GNU make, GZIP, Perl.

# Installation
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("netboost", version = "3.12")
```

# Example
```R
browseVignettes("netboost")
```

# Contact
If you have any issues using the package then please get in touch with Pascal Schlosser (pascal.schlosser at uniklinik-freiburg.de).
Bug reports etc are most welcome, we want the package to be easy to use for everyone!