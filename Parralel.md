Paralellizing QuaSSE - Options and Approaches
========================================================

From the tutorial, begin with an unparallelized version. 

Well start with a simulated tree. The tree simulation differs slightly from the other methods, because
there is no longer a canonical argument list (speciation and extinction rates are arbitrary functions of the
character state). Here is a set of functions; speciation rate is a sigmoidal function which ranges from 0.1
to 0.2 with an inflection point at x = 0, extinction is constant at rate 0.03, and the model of character
evolution is Brownian motion with diffusion parameter 0.025.



```r
# set knitr options
opts_chunk$set(message = FALSE, warning = FALSE)
```



```r
require(diversitree)
lambda <- function(x) sigmoid.x(x, 0.1, 0.2, 0, 2.5)
mu <- function(x) constant.x(x, 0.03)
char <- make.brownian.with.drift(0, 0.025)
```



```r
set.seed(1)
phy <- tree.quasse(c(lambda, mu, char), max.taxa = 15, x0 = 0, single.lineage = FALSE)
```


We need to specify the standard deviation for the states; here I will just assume that all taxa have a state
standard deviation of 1=200


```r
states <- phy$tip.state
states.sd <- 1/200
```


Then, build the likelihood as usual. The difference compared with other models is that we have to specify
the speciation and extinction functions (here, sigmoid.x and constant.x, respectively). There are a
number of other provided functions (see ?constant.x for a list), but any function that takes x as the first
argument may be used.


```r
lik <- make.quasse(phy, states, states.sd, sigmoid.x, constant.x)
```


This can be used in ML calculations as usual. There is a starting.point.quasse function that may
be useful in selecting sensible starting points, but some effort is still required to convert this into a full
vector as it just returns constant rate speciation, extinction, and diffusion rates.


```r
p <- starting.point.quasse(phy, states)
p
```

```
##    lambda        mu diffusion 
##   0.16108   0.02569   0.03164
```


Lets ignore drift: the argument list we need is:


```r
lik.nodrift <- constrain(lik, drift ~ 0)
argnames(lik.nodrift)
```

```
## [1] "l.y0"      "l.y1"      "l.xmid"    "l.r"       "m.c"       "diffusion"
```


A sensible starting point here might be


```r
p.start <- c(p[1], p[1], mean(states), 1, p[2:3])
names(p.start) <- argnames(lik.nodrift)
p.start
```

```
##      l.y0      l.y1    l.xmid       l.r       m.c diffusion 
##   0.16108   0.16108   0.57097   1.00000   0.02569   0.03164
```

```r

# set lower bound
lower <- c(0, 0, min(states), -Inf, 0, 0)
```


Then run find.mle, as usual. The control argument here just tells the subplex algorithm to use an
initial step size of 0.1 (rather than 1), which reduces the number of function evaluations somewhat.


```r
time.f <- system.time(fit <- find.mle(lik.nodrift, p.start, control = list(parscale = 0.1), 
    lower = lower, verbose = 0))
```


*Time of the find.mle function*
Function took 585.82 seconds to run
