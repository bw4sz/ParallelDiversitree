Parallel Optomization Test using PPSO
========================================================


Introduction to ppso - Particle Swarm Optomization
====

The package provides a few dummy functions to show it works - here is the rastrigin function that apparently is commonly used to test optimzation (new to me!)

![](http://upload.wikimedia.org/wikipedia/commons/8/8b/Rastrigin_function.png)



```r
require(ppso)
```

```
## Loading required package: ppso
```

```r
# simple application (all file I/O disabled)
result <- optim_pso(objective_function = rastrigin_function, projectfile = NULL, 
    logfile = NULL)
print(result)
```

```
## $value
## [1] -1.988
## 
## $par
## [1] -0.004496  0.023956
## 
## $function_calls
## [1] 200
## 
## $break_flag
## [1] "max iterations reached"
```

```r

# actual minimum -2 at (0,0)
result = optim_dds(objective_function = rastrigin_function, projectfile = NULL, 
    logfile = NULL)
print(result)
```

```
## $value
## [1] -2
## 
## $par
## [1] 0.001772 0.002502
## 
## $function_calls
## [1] 500
## 
## $break_flag
## [1] "max number of function calls reached"
```


The arguments take in a function that is evaluate at x, and a list of parameters.

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


fit0 and fit1, why do i get diff results?
