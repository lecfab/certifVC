# Quality certification for Vertex Cover
CertifVC is a C++ program to quickly find a small vertex cover and a certificate of its quality.
The associated paper presents the methods and algorithms, and discusses the experimental results.
The detailed results are in [this table](datasets.md).

## Installation & compilation
`$ git clone https://github.com/lecfab/certifVC.git`

`$ cd certifVC`

`$ make`

## Running

`$ ./certifVC DATASET -s SOLHEUR -b BOUNDHEUR` (more information with `$ ./certifVC --help`)

- `$ ./certifVC toygraph.txt` to obtain a certification with default heuristics (see below)
- `$ ./certifVC toygraph.txt -c` to output the small vertex cover found by SOLHEUR
- `$ ./certifVC toygraph.txt -s comparison -b matchLB` to use specific solution-heuristic and bounding-heuristic



#### Parameters
-   `DATASET`: an _undirected_ graph representation for nodes [0 to N-1]. It must be a text
file where each line corresponds to an edge of the form `a b`
(i.e. a SPACE b, with a and b long unsigned integers). Warning: if DATASET contains loops (`a a`)
or double edges (`a b`, `b a`), remove them with `$ ./undirect DATASET NEWDATASET`.

-   `-s SOLHEUR`: the following solution-heuristics are available:
    -   `greedyVC` (default): nodes of high degree are added to the cover until all edges are covered
    -   `greedyMatching`: matching of small size obtained by selecting edges incident to high-degree nodes
    -   `comparison`: run greedy vertex cover algorithms with various parameters to find the smallest cover (see details in the code)


-   `-b BOUNDHEUR`: the following bounding-heuristics are available:
    -   `cliqueLB` (default): greedy clique cover starting with nodes of low degree
    -   `matchLB`: bound obtained from a large matching obtained  by selecting edges incident to low-degree nodes
    -   `comparison`: run greedy clique cover algorithms with various parameters to find the highest bound (see details in the code)

#### Options
-   `-c`: write the smallest vertex cover in standard output so that they can be saved in a FILE (`> FILE`) or piped to another PROGRAM (`| PROGRAM`)


## Contributors

Copyright © 2023  Fabrice Lécuyer (fabrice.lecuyer@lip6.fr)

This program is free software: you can redistribute it and/or modify it with the same open-source licence (or a later version).
It is distributed in the hope that it will be useful, but without any warranty. See the [GNU General Public License](LICENCE.md) for more details.
