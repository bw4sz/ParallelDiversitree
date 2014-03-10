DEoptim Methods
========================================================

he R implementation of Differential Evolution (DE), DEoptim, was first published on the Comprehensive R Archive Network (CRAN) in 2005 by David Ardia. Early versions were written in pure R. Since version 2.0-0 (published to CRAN in 2009) the package has relied on an interface to a C implementation of DE, which is significantly faster on most problems as compared to the implementation in pure R. The C interface is in many respects similar to the MS Visual C++ v5.0 implementation of the Differential Evolution algorithm distributed with the book Differential Evolution -- A Practical Approach to Global Optimization by Price, K.V., Storn, R.M., Lampinen J.A, Springer-Verlag, 2006, and found on-line at http://www.icsi.berkeley.edu/~storn/. Since version 2.0-3 the C implementation dynamically allocates the memory required to store the population, removing limitations on the number of members in the population and length of the parameter vectors that may be optimized. 


```r
require(DEoptim)
```

```
## Loading required package: DEoptim
```

```
## Warning: package 'DEoptim' was built under R version 3.0.2
```

```
## DEoptim package Differential Evolution algorithm in R Authors: D. Ardia,
## K. Mullen, B. Peterson and J. Ulrich
```

```r
Genrose <- function(x) {
    ## One generalization of the Rosenbrock banana valley function (n parameters)
    n <- length(x)
    
    ## make it take some time ...
    Sys.sleep(0.001)
    1 + sum(100 * (x[-n]^2 - x[-1])^2 + (x[-1] - 1)^2)
}

# get some run-time on simple problems
maxIt <- 250
n <- 5

oneCore <- system.time(DEoptim(fn = Genrose, lower = rep(-25, n), upper = rep(25, 
    n), control = list(NP = 10 * n, itermax = maxIt)))
```

```
## Iteration: 1 bestvalit: 1272736.421049 bestmemit:   -7.295027   -4.004996   -8.467699   -8.068786   15.255972
## Iteration: 2 bestvalit: 251068.519973 bestmemit:   -3.301509   -1.636951   -2.658846   -8.188993   21.376103
## Iteration: 3 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 4 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 5 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 6 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 7 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 8 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 9 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 10 bestvalit: 9806.187542 bestmemit:    3.025087    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 11 bestvalit: 7228.280817 bestmemit:    0.765661    2.827341    0.427748    0.056221   -3.122956
## Iteration: 12 bestvalit: 7228.280817 bestmemit:    0.765661    2.827341    0.427748    0.056221   -3.122956
## Iteration: 13 bestvalit: 4940.524364 bestmemit:   -0.952127    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 14 bestvalit: 4940.524364 bestmemit:   -0.952127    2.078030   -1.410220   -1.881906    3.613283
## Iteration: 15 bestvalit: 3851.588404 bestmemit:   -1.138887   -0.338844    1.206451   -0.955146    6.250324
## Iteration: 16 bestvalit: 2561.646205 bestmemit:   -1.138887   -0.338844    1.206451    3.242903    6.250324
## Iteration: 17 bestvalit: 1329.210413 bestmemit:    0.008180    0.917078    0.694447    2.372435    2.662058
## Iteration: 18 bestvalit: 1214.281499 bestmemit:    0.343335    0.116887   -1.167184   -0.142652   -2.853183
## Iteration: 19 bestvalit: 270.465473 bestmemit:    0.008180    0.917078    0.694447    1.755610    2.662058
## Iteration: 20 bestvalit: 270.465473 bestmemit:    0.008180    0.917078    0.694447    1.755610    2.662058
## Iteration: 21 bestvalit: 209.204555 bestmemit:   -0.662775    0.917078    0.694447    1.755610    2.662058
## Iteration: 22 bestvalit: 209.204555 bestmemit:   -0.662775    0.917078    0.694447    1.755610    2.662058
## Iteration: 23 bestvalit: 198.355098 bestmemit:    0.487926   -1.114999    1.139594    1.003935    0.993390
## Iteration: 24 bestvalit: 198.355098 bestmemit:    0.487926   -1.114999    1.139594    1.003935    0.993390
## Iteration: 25 bestvalit: 198.355098 bestmemit:    0.487926   -1.114999    1.139594    1.003935    0.993390
## Iteration: 26 bestvalit: 198.355098 bestmemit:    0.487926   -1.114999    1.139594    1.003935    0.993390
## Iteration: 27 bestvalit: 198.355098 bestmemit:    0.487926   -1.114999    1.139594    1.003935    0.993390
## Iteration: 28 bestvalit: 112.511716 bestmemit:   -0.238135   -0.685245    1.139594    1.003935    0.993390
## Iteration: 29 bestvalit: 112.511716 bestmemit:   -0.238135   -0.685245    1.139594    1.003935    0.993390
## Iteration: 30 bestvalit: 112.511716 bestmemit:   -0.238135   -0.685245    1.139594    1.003935    0.993390
## Iteration: 31 bestvalit: 112.511716 bestmemit:   -0.238135   -0.685245    1.139594    1.003935    0.993390
## Iteration: 32 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 33 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 34 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 35 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 36 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 37 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 38 bestvalit: 38.869029 bestmemit:   -0.667418    0.887384    0.432166    0.128765   -0.155186
## Iteration: 39 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 40 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 41 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 42 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 43 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 44 bestvalit: 9.704141 bestmemit:   -0.046216   -0.123185    0.044327   -0.157421    0.057153
## Iteration: 45 bestvalit: 8.157060 bestmemit:   -0.046216   -0.032110    0.044327   -0.157421    0.057153
## Iteration: 46 bestvalit: 8.157060 bestmemit:   -0.046216   -0.032110    0.044327   -0.157421    0.057153
## Iteration: 47 bestvalit: 8.099767 bestmemit:   -0.046216   -0.024636    0.044327   -0.157421    0.057153
## Iteration: 48 bestvalit: 8.099767 bestmemit:   -0.046216   -0.024636    0.044327   -0.157421    0.057153
## Iteration: 49 bestvalit: 8.099767 bestmemit:   -0.046216   -0.024636    0.044327   -0.157421    0.057153
## Iteration: 50 bestvalit: 8.099767 bestmemit:   -0.046216   -0.024636    0.044327   -0.157421    0.057153
## Iteration: 51 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 52 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 53 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 54 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 55 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 56 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 57 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 58 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 59 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 60 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 61 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 62 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 63 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 64 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 65 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 66 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 67 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 68 bestvalit: 5.492287 bestmemit:   -0.046216   -0.024636    0.044327    0.040504    0.057153
## Iteration: 69 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 70 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 71 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 72 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 73 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 74 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 75 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 76 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 77 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 78 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 79 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 80 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 81 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 82 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 83 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 84 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 85 bestvalit: 4.801499 bestmemit:   -0.920817    0.802652    0.737639    0.497879    0.129827
## Iteration: 86 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 87 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 88 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 89 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 90 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 91 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 92 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 93 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 94 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 95 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 96 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 97 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 98 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 99 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 100 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 101 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 102 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 103 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 104 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 105 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 106 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 107 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 108 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 109 bestvalit: 3.082584 bestmemit:   -0.799749    0.656133    0.451218    0.235693    0.079872
## Iteration: 110 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 111 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 112 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 113 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 114 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 115 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 116 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 117 bestvalit: 3.075001 bestmemit:   -0.799749    0.656133    0.451218    0.232229    0.079872
## Iteration: 118 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 119 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 120 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 121 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 122 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 123 bestvalit: 2.789663 bestmemit:   -0.920817    0.802652    0.634732    0.425722    0.129827
## Iteration: 124 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 125 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 126 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 127 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 128 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 129 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 130 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 131 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 132 bestvalit: 2.536916 bestmemit:   -0.917123    0.946751    0.895647    0.857424    0.726748
## Iteration: 133 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 134 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 135 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 136 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 137 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 138 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 139 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 140 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 141 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 142 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 143 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 144 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 145 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 146 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 147 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 148 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 149 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 150 bestvalit: 2.321592 bestmemit:   -0.883050    0.823816    0.680498    0.501955    0.279508
## Iteration: 151 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 152 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 153 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 154 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 155 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 156 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 157 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 158 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 159 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 160 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 161 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 162 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 163 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 164 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 165 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 166 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 167 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 168 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 169 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 170 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 171 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 172 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 173 bestvalit: 1.866582 bestmemit:   -0.941098    0.883212    0.753568    0.549229    0.311057
## Iteration: 174 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 175 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 176 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 177 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 178 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 179 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 180 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 181 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 182 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 183 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 184 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 185 bestvalit: 1.792734 bestmemit:   -1.010118    0.946751    0.932174    0.857424    0.726748
## Iteration: 186 bestvalit: 1.761770 bestmemit:   -0.929247    0.882465    0.781989    0.587837    0.363568
## Iteration: 187 bestvalit: 1.761770 bestmemit:   -0.929247    0.882465    0.781989    0.587837    0.363568
## Iteration: 188 bestvalit: 1.761770 bestmemit:   -0.929247    0.882465    0.781989    0.587837    0.363568
## Iteration: 189 bestvalit: 1.761770 bestmemit:   -0.929247    0.882465    0.781989    0.587837    0.363568
## Iteration: 190 bestvalit: 1.761770 bestmemit:   -0.929247    0.882465    0.781989    0.587837    0.363568
## Iteration: 191 bestvalit: 1.623000 bestmemit:   -0.992584    0.944225    0.916440    0.794660    0.639549
## Iteration: 192 bestvalit: 1.623000 bestmemit:   -0.992584    0.944225    0.916440    0.794660    0.639549
## Iteration: 193 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 194 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 195 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 196 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 197 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 198 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 199 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 200 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 201 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 202 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 203 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 204 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 205 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 206 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 207 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 208 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 209 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 210 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 211 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 212 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 213 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 214 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 215 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 216 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 217 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 218 bestvalit: 1.204277 bestmemit:   -0.963160    0.948346    0.911523    0.829704    0.681058
## Iteration: 219 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 220 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 221 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 222 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 223 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 224 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 225 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 226 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 227 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 228 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 229 bestvalit: 1.204155 bestmemit:   -0.963160    0.949355    0.911523    0.829704    0.681058
## Iteration: 230 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 231 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 232 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 233 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 234 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 235 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 236 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 237 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 238 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 239 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 240 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 241 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 242 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 243 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 244 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 245 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 246 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 247 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 248 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 249 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
## Iteration: 250 bestvalit: 1.176574 bestmemit:   -1.021937    1.037669    1.042176    1.072108    1.152761
```

```r

# using parallel package
withParallel <- system.time(DEoptim(fn = Genrose, lower = rep(-25, n), upper = rep(25, 
    n), control = list(NP = 10 * n, itermax = maxIt, parallelType = 1)))
```

```
## Iteration: 1 bestvalit: 3054430.465629 bestmemit:   -0.177166   -5.680058  -12.026124   -9.598344   23.162236
## Iteration: 2 bestvalit: 805681.478506 bestmemit:   -9.487350    4.620064   -0.313921   -5.718627   16.596290
## Iteration: 3 bestvalit: 445255.743130 bestmemit:   -0.177166   -5.680058    7.518580   -1.690689   23.162236
## Iteration: 4 bestvalit: 391035.401857 bestmemit:   -4.785344    0.880788    7.791014    3.557962   22.943242
## Iteration: 5 bestvalit: 165196.981731 bestmemit:    5.864749    4.620064   -0.313921   -5.718627   16.596290
## Iteration: 6 bestvalit: 148370.884628 bestmemit:    5.864749    4.620064   -0.313921   -5.145297   16.596290
## Iteration: 7 bestvalit: 24282.609772 bestmemit:    3.030234   -2.692641    2.039524   -4.208976   16.360002
## Iteration: 8 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 9 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 10 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 11 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 12 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 13 bestvalit: 10453.982741 bestmemit:    1.432584    3.338995    2.088700    0.398491    2.391072
## Iteration: 14 bestvalit: 6889.683283 bestmemit:    0.561857   -0.824392    2.579261    0.109800    4.591928
## Iteration: 15 bestvalit: 4038.882926 bestmemit:    0.374560   -2.168164   -0.976984    2.391856    5.072608
## Iteration: 16 bestvalit: 4038.882926 bestmemit:    0.374560   -2.168164   -0.976984    2.391856    5.072608
## Iteration: 17 bestvalit: 790.929924 bestmemit:   -1.164850    1.054662    0.838595    3.096571    8.399980
## Iteration: 18 bestvalit: 790.929924 bestmemit:   -1.164850    1.054662    0.838595    3.096571    8.399980
## Iteration: 19 bestvalit: 790.929924 bestmemit:   -1.164850    1.054662    0.838595    3.096571    8.399980
## Iteration: 20 bestvalit: 790.929924 bestmemit:   -1.164850    1.054662    0.838595    3.096571    8.399980
## Iteration: 21 bestvalit: 790.929924 bestmemit:   -1.164850    1.054662    0.838595    3.096571    8.399980
## Iteration: 22 bestvalit: 365.927854 bestmemit:    0.219615   -0.758252    1.332862    0.239223   -0.053373
## Iteration: 23 bestvalit: 365.927854 bestmemit:    0.219615   -0.758252    1.332862    0.239223   -0.053373
## Iteration: 24 bestvalit: 162.669597 bestmemit:    0.975722    0.056769    0.311856    0.286280    0.897666
## Iteration: 25 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 26 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 27 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 28 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 29 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 30 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 31 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 32 bestvalit: 82.538514 bestmemit:   -0.265382    0.056769    0.311856    0.286280    0.897666
## Iteration: 33 bestvalit: 68.299176 bestmemit:   -0.265382    0.056769    0.311856    0.543827    0.897666
## Iteration: 34 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 35 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 36 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 37 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 38 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 39 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 40 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 41 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 42 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 43 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 44 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 45 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 46 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 47 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 48 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 49 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 50 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 51 bestvalit: 40.538998 bestmemit:    0.088822   -0.081721   -0.145327    0.203296   -0.485762
## Iteration: 52 bestvalit: 39.588586 bestmemit:    0.132498    0.465634   -0.078463    0.239223   -0.053373
## Iteration: 53 bestvalit: 39.588586 bestmemit:    0.132498    0.465634   -0.078463    0.239223   -0.053373
## Iteration: 54 bestvalit: 29.133309 bestmemit:    0.132498    0.465634    0.419974    0.239223   -0.053373
## Iteration: 55 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 56 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 57 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 58 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 59 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 60 bestvalit: 7.575871 bestmemit:   -0.799078    0.705191    0.407368   -0.003043    0.027888
## Iteration: 61 bestvalit: 6.367517 bestmemit:   -0.799078    0.705191    0.407368    0.301767    0.027888
## Iteration: 62 bestvalit: 6.367517 bestmemit:   -0.799078    0.705191    0.407368    0.301767    0.027888
## Iteration: 63 bestvalit: 6.367517 bestmemit:   -0.799078    0.705191    0.407368    0.301767    0.027888
## Iteration: 64 bestvalit: 6.367517 bestmemit:   -0.799078    0.705191    0.407368    0.301767    0.027888
## Iteration: 65 bestvalit: 6.367517 bestmemit:   -0.799078    0.705191    0.407368    0.301767    0.027888
## Iteration: 66 bestvalit: 5.425288 bestmemit:    0.793768    0.705785    0.536960    0.423233    0.111445
## Iteration: 67 bestvalit: 5.425288 bestmemit:    0.793768    0.705785    0.536960    0.423233    0.111445
## Iteration: 68 bestvalit: 5.425288 bestmemit:    0.793768    0.705785    0.536960    0.423233    0.111445
## Iteration: 69 bestvalit: 4.957194 bestmemit:    0.633852    0.506979    0.291838    0.064396    0.042863
## Iteration: 70 bestvalit: 4.957194 bestmemit:    0.633852    0.506979    0.291838    0.064396    0.042863
## Iteration: 71 bestvalit: 4.957194 bestmemit:    0.633852    0.506979    0.291838    0.064396    0.042863
## Iteration: 72 bestvalit: 4.957194 bestmemit:    0.633852    0.506979    0.291838    0.064396    0.042863
## Iteration: 73 bestvalit: 4.957194 bestmemit:    0.633852    0.506979    0.291838    0.064396    0.042863
## Iteration: 74 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 75 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 76 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 77 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 78 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 79 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 80 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 81 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 82 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 83 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 84 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 85 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 86 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 87 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 88 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 89 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 90 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 91 bestvalit: 3.851106 bestmemit:    0.714077    0.506979    0.291838    0.064396    0.042863
## Iteration: 92 bestvalit: 3.732403 bestmemit:   -0.837165    0.698560    0.545162    0.276539    0.168670
## Iteration: 93 bestvalit: 3.732403 bestmemit:   -0.837165    0.698560    0.545162    0.276539    0.168670
## Iteration: 94 bestvalit: 3.732403 bestmemit:   -0.837165    0.698560    0.545162    0.276539    0.168670
## Iteration: 95 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 96 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 97 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 98 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 99 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 100 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 101 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 102 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 103 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 104 bestvalit: 3.624888 bestmemit:   -0.828022    0.705191    0.441299    0.245266    0.027888
## Iteration: 105 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 106 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 107 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 108 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 109 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 110 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 111 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 112 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 113 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 114 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 115 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 116 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 117 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 118 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 119 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 120 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 121 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 122 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 123 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 124 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 125 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 126 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 127 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 128 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 129 bestvalit: 2.386553 bestmemit:    1.033095    1.022789    1.103292    1.130357    1.270745
## Iteration: 130 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 131 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 132 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 133 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 134 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 135 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 136 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 137 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 138 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 139 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 140 bestvalit: 1.843690 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.142935
## Iteration: 141 bestvalit: 1.787904 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.138368
## Iteration: 142 bestvalit: 1.787904 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.138368
## Iteration: 143 bestvalit: 1.787904 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.138368
## Iteration: 144 bestvalit: 1.787904 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.138368
## Iteration: 145 bestvalit: 1.787904 bestmemit:    0.970189    1.007052    1.020660    1.039704    1.138368
## Iteration: 146 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 147 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 148 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 149 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 150 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 151 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 152 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 153 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 154 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 155 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 156 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 157 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 158 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 159 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 160 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 161 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 162 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 163 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 164 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 165 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 166 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 167 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 168 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 169 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 170 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 171 bestvalit: 1.431920 bestmemit:    1.033095    1.060997    1.076439    1.130357    1.270745
## Iteration: 172 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 173 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 174 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 175 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 176 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 177 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 178 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 179 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 180 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 181 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 182 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 183 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 184 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 185 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 186 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 187 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 188 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 189 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 190 bestvalit: 1.285466 bestmemit:    0.994816    1.019667    1.022701    1.060107    1.158008
## Iteration: 191 bestvalit: 1.237375 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.138730
## Iteration: 192 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 193 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 194 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 195 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 196 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 197 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 198 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 199 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 200 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 201 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 202 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 203 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 204 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 205 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 206 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 207 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 208 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 209 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 210 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 211 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 212 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 213 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 214 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 215 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 216 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 217 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 218 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 219 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 220 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 221 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 222 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 223 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 224 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 225 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 226 bestvalit: 1.191382 bestmemit:    0.992777    1.016625    1.026571    1.077809    1.160326
## Iteration: 227 bestvalit: 1.130713 bestmemit:    0.986348    1.003578    0.998724    0.990096    0.995726
## Iteration: 228 bestvalit: 1.130713 bestmemit:    0.986348    1.003578    0.998724    0.990096    0.995726
## Iteration: 229 bestvalit: 1.130713 bestmemit:    0.986348    1.003578    0.998724    0.990096    0.995726
## Iteration: 230 bestvalit: 1.130713 bestmemit:    0.986348    1.003578    0.998724    0.990096    0.995726
## Iteration: 231 bestvalit: 1.130713 bestmemit:    0.986348    1.003578    0.998724    0.990096    0.995726
## Iteration: 232 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 233 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 234 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 235 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 236 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 237 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 238 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 239 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 240 bestvalit: 1.082913 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.095374
## Iteration: 241 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 242 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 243 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 244 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 245 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 246 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 247 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 248 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 249 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
## Iteration: 250 bestvalit: 1.081560 bestmemit:    1.005701    1.012094    1.023313    1.035111    1.046396
```

```r

require(iterators)
```

```
## Loading required package: iterators
```

```
## Warning: package 'iterators' was built under R version 3.0.2
```

```r

# using foreach
require(doSNOW)
```

```
## Loading required package: doSNOW Loading required package: foreach Loading
## required package: snow
## 
## Attaching package: 'snow'
## 
## The following objects are masked from 'package:parallel':
## 
## clusterApply, clusterApplyLB, clusterCall, clusterEvalQ, clusterExport,
## clusterMap, clusterSplit, makeCluster, parApply, parCapply, parLapply,
## parRapply, parSapply, splitIndices, stopCluster
```

```r
cl <- makeCluster(8, "SOCK")
registerDoSNOW(cl)

withParallelforeach <- system.time(DEoptim(fn = Genrose, lower = rep(-25, n), 
    upper = rep(25, n), control = list(NP = 10 * n, itermax = maxIt, parallelType = 2)))
```

```
## Iteration: 1 bestvalit: 131439.969794 bestmemit:   -3.370363   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 2 bestvalit: 131439.969794 bestmemit:   -3.370363   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 3 bestvalit: 131439.969794 bestmemit:   -3.370363   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 4 bestvalit: 131439.969794 bestmemit:   -3.370363   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 5 bestvalit: 123706.651922 bestmemit:    2.924959   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 6 bestvalit: 123706.651922 bestmemit:    2.924959   -3.832626   -4.383272   -3.684877   -0.377452
## Iteration: 7 bestvalit: 98907.872555 bestmemit:    4.031989   -4.308610    3.327579   -1.864102   -9.374932
## Iteration: 8 bestvalit: 5180.539025 bestmemit:    2.589026    2.746886    1.865523    2.016382    2.777115
## Iteration: 9 bestvalit: 5180.539025 bestmemit:    2.589026    2.746886    1.865523    2.016382    2.777115
## Iteration: 10 bestvalit: 5025.996394 bestmemit:    2.589026    2.746886    1.865523    2.016382    4.266230
## Iteration: 11 bestvalit: 5025.996394 bestmemit:    2.589026    2.746886    1.865523    2.016382    4.266230
## Iteration: 12 bestvalit: 2872.553837 bestmemit:    2.589026    1.959580    1.865523    2.016382    4.266230
## Iteration: 13 bestvalit: 1305.921930 bestmemit:   -1.402662   -1.129974   -0.089081   -0.582282    1.414707
## Iteration: 14 bestvalit: 1305.921930 bestmemit:   -1.402662   -1.129974   -0.089081   -0.582282    1.414707
## Iteration: 15 bestvalit: 889.763047 bestmemit:   -0.437789    1.161754    1.126037   -1.527466    2.335239
## Iteration: 16 bestvalit: 889.763047 bestmemit:   -0.437789    1.161754    1.126037   -1.527466    2.335239
## Iteration: 17 bestvalit: 889.763047 bestmemit:   -0.437789    1.161754    1.126037   -1.527466    2.335239
## Iteration: 18 bestvalit: 889.763047 bestmemit:   -0.437789    1.161754    1.126037   -1.527466    2.335239
## Iteration: 19 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 20 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 21 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 22 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 23 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 24 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 25 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 26 bestvalit: 77.544604 bestmemit:    1.076674    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 27 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 28 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 29 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 30 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 31 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 32 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 33 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 34 bestvalit: 77.357903 bestmemit:    0.191847    0.597184    0.020070   -0.185756   -0.470200
## Iteration: 35 bestvalit: 62.368075 bestmemit:   -0.733497    0.009045    0.422805    0.499323    0.432452
## Iteration: 36 bestvalit: 62.368075 bestmemit:   -0.733497    0.009045    0.422805    0.499323    0.432452
## Iteration: 37 bestvalit: 54.232168 bestmemit:    0.191847    0.465862    0.157895   -0.185756   -0.470200
## Iteration: 38 bestvalit: 14.087212 bestmemit:    0.660891    0.166645   -0.030116    0.141635    0.012141
## Iteration: 39 bestvalit: 14.087212 bestmemit:    0.660891    0.166645   -0.030116    0.141635    0.012141
## Iteration: 40 bestvalit: 14.087212 bestmemit:    0.660891    0.166645   -0.030116    0.141635    0.012141
## Iteration: 41 bestvalit: 14.087212 bestmemit:    0.660891    0.166645   -0.030116    0.141635    0.012141
## Iteration: 42 bestvalit: 14.087212 bestmemit:    0.660891    0.166645   -0.030116    0.141635    0.012141
## Iteration: 43 bestvalit: 13.207969 bestmemit:   -0.784743    0.364045    0.248965    0.057806    0.141861
## Iteration: 44 bestvalit: 12.514084 bestmemit:    0.660891    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 45 bestvalit: 12.514084 bestmemit:    0.660891    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 46 bestvalit: 9.187158 bestmemit:    0.402582    0.275078    0.211368   -0.074550    0.075596
## Iteration: 47 bestvalit: 8.672048 bestmemit:    0.032561    0.126292    0.168573   -0.008563    0.051698
## Iteration: 48 bestvalit: 8.672048 bestmemit:    0.032561    0.126292    0.168573   -0.008563    0.051698
## Iteration: 49 bestvalit: 6.935187 bestmemit:   -0.581629    0.364045    0.248965    0.057806    0.141861
## Iteration: 50 bestvalit: 6.935187 bestmemit:   -0.581629    0.364045    0.248965    0.057806    0.141861
## Iteration: 51 bestvalit: 5.597497 bestmemit:   -0.477844    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 52 bestvalit: 5.597497 bestmemit:   -0.477844    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 53 bestvalit: 5.597497 bestmemit:   -0.477844    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 54 bestvalit: 5.597497 bestmemit:   -0.477844    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 55 bestvalit: 5.597497 bestmemit:   -0.477844    0.166645   -0.030116   -0.027624    0.012141
## Iteration: 56 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 57 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 58 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 59 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 60 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 61 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 62 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 63 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 64 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 65 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 66 bestvalit: 5.521360 bestmemit:   -0.477844    0.175005   -0.030116   -0.027624    0.012141
## Iteration: 67 bestvalit: 5.503235 bestmemit:    0.628684    0.472420    0.267015    0.118932   -0.072278
## Iteration: 68 bestvalit: 5.085138 bestmemit:   -0.506791    0.275078    0.131514   -0.000777    0.075596
## Iteration: 69 bestvalit: 5.085138 bestmemit:   -0.506791    0.275078    0.131514   -0.000777    0.075596
## Iteration: 70 bestvalit: 5.085138 bestmemit:   -0.506791    0.275078    0.131514   -0.000777    0.075596
## Iteration: 71 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 72 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 73 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 74 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 75 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 76 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 77 bestvalit: 4.624710 bestmemit:    0.581337    0.392678    0.120397    0.054624    0.029938
## Iteration: 78 bestvalit: 4.390422 bestmemit:    0.774588    0.575319    0.404318    0.096759    0.003716
## Iteration: 79 bestvalit: 4.390422 bestmemit:    0.774588    0.575319    0.404318    0.096759    0.003716
## Iteration: 80 bestvalit: 4.390422 bestmemit:    0.774588    0.575319    0.404318    0.096759    0.003716
## Iteration: 81 bestvalit: 4.390422 bestmemit:    0.774588    0.575319    0.404318    0.096759    0.003716
## Iteration: 82 bestvalit: 4.200306 bestmemit:    0.774588    0.575319    0.404318    0.110922    0.003716
## Iteration: 83 bestvalit: 4.200306 bestmemit:    0.774588    0.575319    0.404318    0.110922    0.003716
## Iteration: 84 bestvalit: 4.200306 bestmemit:    0.774588    0.575319    0.404318    0.110922    0.003716
## Iteration: 85 bestvalit: 4.200306 bestmemit:    0.774588    0.575319    0.404318    0.110922    0.003716
## Iteration: 86 bestvalit: 3.648384 bestmemit:    0.750229    0.576403    0.303000    0.080066    0.036578
## Iteration: 87 bestvalit: 3.648384 bestmemit:    0.750229    0.576403    0.303000    0.080066    0.036578
## Iteration: 88 bestvalit: 3.648384 bestmemit:    0.750229    0.576403    0.303000    0.080066    0.036578
## Iteration: 89 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 90 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 91 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 92 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 93 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 94 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 95 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 96 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 97 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 98 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 99 bestvalit: 3.482778 bestmemit:    0.774588    0.575319    0.343274    0.110922    0.003716
## Iteration: 100 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 101 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 102 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 103 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 104 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 105 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 106 bestvalit: 3.220269 bestmemit:    0.828312    0.681536    0.427594    0.167351    0.035984
## Iteration: 107 bestvalit: 3.158843 bestmemit:    0.833896    0.716731    0.506398    0.206699    0.052333
## Iteration: 108 bestvalit: 3.158843 bestmemit:    0.833896    0.716731    0.506398    0.206699    0.052333
## Iteration: 109 bestvalit: 3.158843 bestmemit:    0.833896    0.716731    0.506398    0.206699    0.052333
## Iteration: 110 bestvalit: 3.158843 bestmemit:    0.833896    0.716731    0.506398    0.206699    0.052333
## Iteration: 111 bestvalit: 3.158843 bestmemit:    0.833896    0.716731    0.506398    0.206699    0.052333
## Iteration: 112 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 113 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 114 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 115 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 116 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 117 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 118 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 119 bestvalit: 3.151670 bestmemit:    0.826647    0.656248    0.420564    0.190466    0.031485
## Iteration: 120 bestvalit: 3.145809 bestmemit:    0.828312    0.681536    0.436935    0.222067    0.035984
## Iteration: 121 bestvalit: 3.145809 bestmemit:    0.828312    0.681536    0.436935    0.222067    0.035984
## Iteration: 122 bestvalit: 3.145809 bestmemit:    0.828312    0.681536    0.436935    0.222067    0.035984
## Iteration: 123 bestvalit: 3.022368 bestmemit:    0.828312    0.681536    0.491135    0.222067    0.035984
## Iteration: 124 bestvalit: 3.022368 bestmemit:    0.828312    0.681536    0.491135    0.222067    0.035984
## Iteration: 125 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 126 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 127 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 128 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 129 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 130 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 131 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 132 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 133 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 134 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 135 bestvalit: 2.917360 bestmemit:    0.926011    0.817305    0.655604    0.383191    0.103581
## Iteration: 136 bestvalit: 2.916202 bestmemit:    0.828312    0.700008    0.491135    0.222067    0.048243
## Iteration: 137 bestvalit: 2.912317 bestmemit:    0.848467    0.698900    0.484185    0.224635    0.052039
## Iteration: 138 bestvalit: 2.912317 bestmemit:    0.848467    0.698900    0.484185    0.224635    0.052039
## Iteration: 139 bestvalit: 2.912317 bestmemit:    0.848467    0.698900    0.484185    0.224635    0.052039
## Iteration: 140 bestvalit: 2.905269 bestmemit:    0.828312    0.700008    0.491135    0.222067    0.057439
## Iteration: 141 bestvalit: 2.767048 bestmemit:    0.862210    0.743409    0.532083    0.274880    0.060090
## Iteration: 142 bestvalit: 2.767048 bestmemit:    0.862210    0.743409    0.532083    0.274880    0.060090
## Iteration: 143 bestvalit: 2.767048 bestmemit:    0.862210    0.743409    0.532083    0.274880    0.060090
## Iteration: 144 bestvalit: 2.755162 bestmemit:    0.862210    0.743409    0.537962    0.274880    0.060090
## Iteration: 145 bestvalit: 2.755162 bestmemit:    0.862210    0.743409    0.537962    0.274880    0.060090
## Iteration: 146 bestvalit: 2.755162 bestmemit:    0.862210    0.743409    0.537962    0.274880    0.060090
## Iteration: 147 bestvalit: 2.755162 bestmemit:    0.862210    0.743409    0.537962    0.274880    0.060090
## Iteration: 148 bestvalit: 2.755162 bestmemit:    0.862210    0.743409    0.537962    0.274880    0.060090
## Iteration: 149 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 150 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 151 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 152 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 153 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 154 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 155 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 156 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 157 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 158 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 159 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 160 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 161 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 162 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 163 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 164 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 165 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 166 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 167 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 168 bestvalit: 2.599676 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.082492
## Iteration: 169 bestvalit: 2.531017 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.110022
## Iteration: 170 bestvalit: 2.531017 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.110022
## Iteration: 171 bestvalit: 2.531017 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.110022
## Iteration: 172 bestvalit: 2.531017 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.110022
## Iteration: 173 bestvalit: 2.531017 bestmemit:    0.861263    0.749863    0.564184    0.315736    0.110022
## Iteration: 174 bestvalit: 2.429707 bestmemit:    0.889243    0.807875    0.643964    0.401207    0.130048
## Iteration: 175 bestvalit: 2.429707 bestmemit:    0.889243    0.807875    0.643964    0.401207    0.130048
## Iteration: 176 bestvalit: 2.429707 bestmemit:    0.889243    0.807875    0.643964    0.401207    0.130048
## Iteration: 177 bestvalit: 2.428087 bestmemit:    0.889243    0.807875    0.643964    0.397557    0.130048
## Iteration: 178 bestvalit: 2.428087 bestmemit:    0.889243    0.807875    0.643964    0.397557    0.130048
## Iteration: 179 bestvalit: 2.428087 bestmemit:    0.889243    0.807875    0.643964    0.397557    0.130048
## Iteration: 180 bestvalit: 2.428087 bestmemit:    0.889243    0.807875    0.643964    0.397557    0.130048
## Iteration: 181 bestvalit: 2.428087 bestmemit:    0.889243    0.807875    0.643964    0.397557    0.130048
## Iteration: 182 bestvalit: 2.275551 bestmemit:    0.882395    0.810730    0.659146    0.457906    0.194615
## Iteration: 183 bestvalit: 2.275551 bestmemit:    0.882395    0.810730    0.659146    0.457906    0.194615
## Iteration: 184 bestvalit: 2.275551 bestmemit:    0.882395    0.810730    0.659146    0.457906    0.194615
## Iteration: 185 bestvalit: 1.963190 bestmemit:    0.939010    0.852584    0.718292    0.518021    0.267306
## Iteration: 186 bestvalit: 1.963190 bestmemit:    0.939010    0.852584    0.718292    0.518021    0.267306
## Iteration: 187 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 188 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 189 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 190 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 191 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 192 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 193 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 194 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 195 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 196 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 197 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 198 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 199 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 200 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 201 bestvalit: 1.914432 bestmemit:    0.910732    0.845835    0.718242    0.522223    0.263393
## Iteration: 202 bestvalit: 1.912022 bestmemit:    0.910732    0.845835    0.720213    0.522223    0.263393
## Iteration: 203 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 204 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 205 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 206 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 207 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 208 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 209 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 210 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 211 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 212 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 213 bestvalit: 1.772089 bestmemit:    0.946137    0.900008    0.817952    0.625886    0.381178
## Iteration: 214 bestvalit: 1.683917 bestmemit:    0.954685    0.897420    0.811903    0.627188    0.408479
## Iteration: 215 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 216 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 217 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 218 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 219 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 220 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 221 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 222 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 223 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 224 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 225 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 226 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 227 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 228 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 229 bestvalit: 1.580584 bestmemit:    0.954685    0.897420    0.811903    0.641453    0.408479
## Iteration: 230 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 231 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 232 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 233 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 234 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 235 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 236 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 237 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 238 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 239 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 240 bestvalit: 1.562449 bestmemit:    0.954685    0.897420    0.808591    0.641453    0.408479
## Iteration: 241 bestvalit: 1.518312 bestmemit:    0.968325    0.910629    0.832766    0.692845    0.506526
## Iteration: 242 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 243 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 244 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 245 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 246 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 247 bestvalit: 1.380838 bestmemit:    0.982451    0.931715    0.868790    0.768326    0.603591
## Iteration: 248 bestvalit: 1.336056 bestmemit:    0.988075    0.940730    0.893045    0.799570    0.629266
## Iteration: 249 bestvalit: 1.313139 bestmemit:    0.979993    0.933419    0.886815    0.784203    0.611635
## Iteration: 250 bestvalit: 1.313139 bestmemit:    0.979993    0.933419    0.886815    0.784203    0.611635
```

```r
stopCluster(cl)

# compare speeds
spe <- list(oneCore, withParallel, withParallelforeach)
names(spe) <- c("oneCore", "withParallel", "withParallelforeach")
spe
```

```
## $oneCore
##    user  system elapsed 
##    0.01    0.00   12.57 
## 
## $withParallel
##    user  system elapsed 
##    0.84    0.33    6.18 
## 
## $withParallelforeach
##    user  system elapsed 
##    3.46    0.36    5.99
```



