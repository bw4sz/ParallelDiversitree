Writing likelihood functions in R
========================================================

I'm having trouble with the syntax for writing likelihood functions in R. This is important since we want to create lots of dummy examples of small phylogenies, run them both on diversitree and any parallel approach we come up with.

Creating own function, let's start simple
Lets take a look with some data found in ?mle2


```r
x = 0:10
y = c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8)
d = data.frame(x, y)
```


Here, were fitting a possion model that depends on an intercept term plus a linear term. The exp() is mainly there to make sure that the value of lambda doesnt go negative, which isnt allowed (it would imply a negative number of occurrences for our outcome of interest).


```r
require(bbmle)
```

```
## Loading required package: bbmle
## Loading required package: stats4
```

```r
# using R's native optim function()
fit0 = mle2(y ~ dpois(lambda = exp(intercept + slope * x)), start = list(intercept = mean(y), 
    slope = 0), data = d)

# now the same exact code written as a function, how do i do this?
simpleF <- function(intercept = mean(y), slope = 0) {
    -sum(dpois(y, lambda = exp(intercept + slope * x)), log = TRUE)
}

# test
simpleF()
```

```
## [1] -1
```

```r

fit1 = mle2(simpleF, start = list(intercept = mean(y), slope = 0))

fit1
```

```
## 
## Call:
## mle2(minuslogl = simpleF, start = list(intercept = mean(y), slope = 0))
## 
## Coefficients:
## intercept     slope 
##     11.55      0.00 
## 
## Log-likelihood: 1
```

```r
fit0
```

```
## 
## Call:
## mle2(minuslogl = y ~ dpois(lambda = exp(intercept + slope * x)), 
##     start = list(intercept = mean(y), slope = 0), data = d)
## 
## Coefficients:
## intercept     slope 
##    3.0949   -0.1524 
## 
## Log-likelihood: -28.95
```


Clearly we have a syntax issue. Thoughts?
===============================
