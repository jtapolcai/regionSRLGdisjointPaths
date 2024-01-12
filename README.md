# Source code for paper Erika Bérczi-Kovács, Péter Gyimesi, Balázs Vass, János Tapolcai, *"Efficient Algorithm for Region-Disjoint Survivable Routing in Backbone Networks"*, Infocom 2024

The code was written by Péter Gyimesi

## compile
```
g++ kod.cpp -o srlg-path -std=c++11
```

To find only two disjoint path:
```
g++ kod.cpp -o srlg-path2 -std=c++11 -D TWOPATH
```

To find k with binary search
```
g++ kod.cpp -o srlg-path3 -std=c++11 -D BINSEARCHCUT -D FAST
```

An implementation of Dervish algorithm from *B. Vass, E. Bérczi-Kovács, A. Barabás, Z. L. Hajdú, and J. Tapolcai,
“Polynomial-time algorithm for the regional SRLG-disjoint paths prob-
lem,” in Proc. IEEE INFOCOM, London, United Kingdom, May 2022.*
```
g++ kod.cpp -o srlg-path4 -std=c++11 -D GREEDY -D FAST
```

## Input format:

n - number of vertices

Following n lines:

x_i, y_i, id_i - coordinates and ID of vertices

Next line: s, t - start and end vertices

m - number of edges

Next m lines:

a_i, b_i, id_i - endpoints and ID of each edge

k - number of SRLGs (Shared Risk Link Groups)

Next k lines:

x, a_1, a_2 ... a_x, - edges traversed by each SRLG 

in the same line we define y the capacity of the SRLG

## Output format:

Writes ut the results in an xml format.

Determines if there are infinitely many paths; otherwise, it lists k - the number of maximum SRLG-independent paths (if k is not very large)

Next k lines: 
paths length in hops, the nodes along the path x_1=s, x_2, ... x_p=t
