# Combinatorial discrepancy algorithms

Experimental implementation of the Lovett-Meka combinatorial discrepancy algorithm (see the [original paper](https://epubs.siam.org/doi/abs/10.1137/130929400) for reference) and its derandomization by Levy, Ramadas and Rothvoss (see the [original paper](https://link.springer.com/chapter/10.1007/978-3-319-59250-3_31) for reference).

## Lovett-Meka Algorithm

```c++
lm(SetSystem ss, vector<float> cst);
```

This function defined in [discrepancy.cpp](discrepancy.cpp) executes the Lovett-Meka algorithm on the set system `ss` with constraints `cst` (`cst` and `ss.sets` must have the same size).

### Implementation choices

- The algorithm is implemented so that it doesn't fail, that is, if a coordinate of the coloring exceeds 1 or -1, it is simply snapped back to 1 or -1

## Levy-Ramadas-Rothvoss algorithm

```c++
Coloring lrr(SetSystem ss, vector<float> constraints, vector<float> (*incrf)(vector<vector<float>>))
```

This function defined in [discrepancy.cpp](discrepancy.cpp) executes the Levy-Ramadas-Rothvoss algorithm on the set system `ss` with constraints `cst` (`cst` and `ss.sets` must have the same size). 

The function `incrf` is used to compute the increment z(t) at each iteration from a list containing a basis of U(t).

We provide 2 functions `simple_sum` that sums the basis vectors and `random_sum` that sums the basis vector with a random gaussian coefficient.

### Implementation choices

- If the dimension from which the increment has to be chosen is 0 (which might happen when there are less that 16 remaining colors to chose), we stop the execution of the algorithm and chose the remaining colors to be of the sign of the color in the partial coloring.
- The arbitrary increment is chosen to be the sum of the basis vector of the space from which it is selected (more experiment coming on different ways to select this vector)

The implementation is experimental and not documented at this stage. Feel free to reach out for explanations.
