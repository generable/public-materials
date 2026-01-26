::: {#tbl-heteroscedastic-sampling}

| **target acceptance rate** | **non centered ESS** | **non centered #evals** |                                                           **non centered #divergent** | **partially centered ESS** | **partially centered #evals** |                    **partially centered #divergent** |
|---------------------------:|---------------------:|------------------------:|--------------------------------------------------------------------------------------:|---------------------------:|------------------------------:|-----------------------------------------------------:|
|                        0.8 |                  281 |                    9.3M | 2922 (1x0+5x1+...) |                       9.2k |                          4.1M | 98 (12x0+12x1+...) |
|                        0.9 |                 9.6k |                     14M |                                                125 (7x0+7x1+...) |                        10k |                          4.6M |            63 (18x0+6x1+...) |
|                       0.99 |                  10k |                     26M |                                                                          3 (37x0+3x1) |                       9.6k |                          9.3M |                                         3 (37x0+3x1) |

**Sampling statistics for the heteroscedastic model**. 40 independent chains with 1000 post-warmup draws each are used to sample from the posterior for varying target acceptance rates (left-most column, 0.8, 0.9 and 0.99). We report the minimal ESS (2nd and 5th columns), the total number of gradient evaluations (3rd and 6th columns), and the number of divergences (as "total (histogram over chains)", 4th and 7th columns) for the non-centered parametrizations (columns 2, 3, and 4) and the partially centered parametrizations (columns 5, 6, and 7). While the individual chains are independent, all partially centered chains use the same optimized parametrization due to the initial long chain from the main text (if apllicable) with 10,000 post warm-up draws.


:::
