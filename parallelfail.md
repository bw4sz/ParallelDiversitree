Diagnosing the paralell ability of the current quasse function
========================================================



```r
require(diversitree)
```

```
## Loading required package: diversitree
## Loading required package: deSolve
## Loading required package: ape
## Loading required package: subplex
## Loading required package: Rcpp
```

```r
lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
```




```r
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa = 15, x0 = 0, single.lineage = FALSE)
plot(phy)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2.png) 




```r
states <- phy$tip.state
states.sd <- 1/200
lik <- make.quasse(phy, states, states.sd, sigmoid.x, constant.x)
lik.nodrift <- constrain(lik, drift ~ 0)
argnames(lik.nodrift)
```

```
## [1] "l.y0"      "l.y1"      "l.xmid"    "l.r"       "m.c"       "diffusion"
```




```r
# provide pars
pars <- c(0.1, 0.2, 0, 2.5, 0.03, 0.01)

lik.nodrift(pars)
```

[1] -62.04

```r

require(foreach)
```

```
## Loading required package: foreach
```

```r
require(doSNOW)
```

```
## Loading required package: doSNOW
## Loading required package: iterators
## Loading required package: snow
```

```r

cl <- makeCluster(2, "SOCK")
registerDoSNOW(cl)
foreach(i = 1:5) %dopar% {
    lik.nodrift(pars)
}
```

```
## Error: task 1 failed - "Corrupt QuaSSE integrator: ptr is NULL (are you
## using multicore?)"
```


