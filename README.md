# Combinatorial discrepancy algorithms

Experimental implementation of the Lovett-Meka combinatorial discrepancy algorithm (see the [original paper](https://epubs.siam.org/doi/abs/10.1137/130929400) for reference) and its derandomization by Levy, Ramadas and Rothvoss (see the [original paper](https://link.springer.com/chapter/10.1007/978-3-319-59250-3_31) for reference).

## Lovett-Meka algorithm implementation choices

- The algorithm is implemented so that it doesn't fail, that is, if a coordinate of the coloring exceeds 1 or -1, it is simply snapped back to 1 or -1

## Levy-Ramadas-Rothvoss algorithm implementation choices

- If the dimension from which the increment has to be chosen is 0 (which might happen when there are less that 16 remaining colors to chose), we stop the execution of the algorithm and chose the remaining colors to be of the sign of the color in the partial coloring.
- The arbitrary increment is chosen to be the sum of the basis vector of the space from which it is selected (more experiment coming on different ways to select this vector)

The implementation is experimental and not documented at this stage. Feel free to reach out for explanations.
