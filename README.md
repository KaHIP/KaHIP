KaHIP v2.10  [![Codacy Badge](https://api.codacy.com/project/badge/Grade/b4dc07be37694dc4841beba6ba038c19)](https://app.codacy.com/manual/schulzchristian/KaHIP?utm_source=github.com&utm_medium=referral&utm_content=schulzchristian/KaHIP&utm_campaign=Badge_Grade_Dashboard)
[![Build Status](https://travis-ci.org/schulzchristian/KaHIP.svg?branch=master)](https://travis-ci.org/schulzchristian/KaHIP) 
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![FOSSA Status](https://app.fossa.com/api/projects/git%2Bgithub.com%2Fschulzchristian%2FKaHIP.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2Fschulzchristian%2FKaHIP?ref=badge_shield)
=====

The graph partitioning framework KaHIP -- Karlsruhe High Quality Partitioning.

The graph partitioning problem asks for a division of a graph's node set into k equally sized blocks such that the number of edges that run between the blocks is minimized. KaHIP is a family of graph partitioning programs. It includes KaFFPa (Karlsruhe Fast Flow Partitioner), which is a multilevel graph partitioning algorithm, in its variants Strong, Eco and Fast, KaFFPaE (KaFFPaEvolutionary) which is a parallel evolutionary algorithm that uses KaFFPa to provide combine and mutation operations, as well as KaBaPE which extends the evolutionary algorithm. Moreover, specialized techniques are included to partition road networks (Buffoon), to output a vertex separator from a given partition as well as techniques geared towards the efficient partitioning of social networks. Here is an overview of our framework:

<p align="center">
<img src="./img/MGPall_en_new.png"
  alt="framework overview"
  width="601" height="558">
</p>

## NEW in v2.10: 



*ParHIP (Parallel High Quality Partitioning):* Our distributed memory parallel partitioning techniques designed to partition hierarchically structured networks such as web graphs or social networks.

*Mapping Algorithms:* Our new algorithms to map the blocks onto processors to minimize overall communication time based on hierarchical partitionings of the task graph and fast local search algorithms.

*Edge Partitioning Algorithms:* Our new algorithms to compute edge partitionings of graphs. 


## Main project site:
https://algo2.iti.kit.edu/documents/kahip/index.html

Installation Notes
=====

Before you can start you need to install the following software packages:

- Scons (https://www.scons.org/)
- OpenMPI (https://www.open-mpi.org/). Note: due to removed progress threads in OpenMPI > 1.8, please use an OpenMPI version < 1.8 or Intel MPI to obtain a scalable parallel algorithm.

Once you installed the packages, just type ./compile.sh. Once you did that you can try to run the following command:

./deploy/kaffpa examples/delaunay_n15.graph --k 2 --preconfiguration=strong

For a description of the graph format please have a look into the manual.



Licence
=====
The program is licenced under MIT licence.
If you publish results using our algorithms, please acknowledge our work by quoting the following paper:

```
@inproceedings{sandersschulz2013,
             AUTHOR = {Sanders, Peter and Schulz, Christian},
             TITLE = {{Think Locally, Act Globally: Highly Balanced Graph Partitioning}},
             BOOKTITLE = {Proceedings of the 12th International Symposium on Experimental Algorithms (SEA'13)},
             SERIES = {LNCS},
             PUBLISHER = {Springer},
             YEAR = {2013},
             VOLUME = {7933},
             PAGES = {164--175}
}
```

If you use our parallel partitioner ParHIP please also cite the following paper:

```
@inproceedings{meyerhenkesandersschulz2017,
             AUTHOR = {Meyerhenke, Henning and Sanders, Peter and Schulz, Christian},
             TITLE = {{Parallel Graph Partitioning for Complex Networks}},
             JOURNAL = {IEEE Transactions on Parallel and Distributed Systems (TPDS)},
             VOLUME = {28},
             NUMBER = {9},
             PAGES = {2625--2638},
             YEAR = {2017}
}
```

If you use mapping algorithm please also cite the following paper:

```
@inproceedings{schulztraeff2017,
             AUTHOR = {Schulz, Christian and Träff, Jesper Larsson},
             TITLE = {{Better Process Mapping and Sparse Quadratic Assignment}},
             BOOKTITLE = {Proceedings of the 16th International Symposium on Experimental Algorithms (SEA'17)},
             PUBLISHER = {Schloss Dagstuhl - Leibniz-Zentrum fuer Informatik},
             VOLUME = {75},
             SERIES = {LIPIcs},
             PAGES = {4:1--4:15},
             YEAR = {2017}
}
```

If you use edge partitioning algorithms please also cite the following paper:

```
@inproceedings{edgepartitioning2019,
             AUTHOR = {Schlag, Sebastian and Schulz, Christian and Seemaier, Daniel and Strash, Darren},
             TITLE = {{Scalable Edge Partitioning}},
             BOOKTITLE = {Proceedings of the 21th Workshop on Algorithm Engineering and Experimentation (ALENEX)},
             PUBLISHER = {SIAM},
             PAGES = {211--225},
             YEAR = {2019}
}
```

Project Contributors (sorted by last name)
=====
Yaroslav Akhremtsev

Roland Glantz
 
Dennis Luxen

Henning Meyerhenke

Ilya Safro

Peter Sanders

Sebastian Schlag

Christian Schulz (maintainer)

Daniel Seemaier

Darren Strash

Jesper Larsson Träff
