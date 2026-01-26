::: {#tbl-homoscedastic-sampling}

| **target acceptance rate** | **non centered ESS** | **non centered #evals** | **non centered #divergent** | **partially centered ESS** | **partially centered #evals** | **partially centered #divergent** |
|---------------------------:|---------------------:|------------------------:|----------------------------:|---------------------------:|------------------------------:|----------------------------------:|
|                        0.8 |                 6.2k |                    1.9M |  42 (27x0+4x1+4x2+4x3+1x18) |                        15k |                          930k |                  4 (37x0+2x1+1x2) |
|                        0.9 |                 8.1k |                    2.0M |           10 (33x0+4x1+3x2) |                        14k |                          980k |                      2 (38x0+2x1) |
|                       0.99 |                 7.5k |                    4.1M |                    0 (40x0) |                        11k |                          2.3M |                          0 (40x0) |

**Sampling statistics for the homoscedastic model**. 40 independent chains with 1000 post-warmup draws each are used to sample from the posterior for varying target acceptance rates (left-most column, 0.8, 0.9 and 0.99). We report the minimal ESS (2nd and 5th columns), the total number of gradient evaluations (3rd and 6th columns), and the number of divergences (as "total (histogram over chains)", 4th and 7th columns) for the non-centered parametrizations (columns 2, 3, and 4) and the partially centered parametrizations (columns 5, 6, and 7). While the individual chains are independent, all partially centered chains use the same optimized parametrization due to the initial long chain from the main text (if apllicable) with 10,000 post warm-up draws.


:::
