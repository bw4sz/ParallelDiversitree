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
```

```
## Warning: Package lhs not installed, lhs_init disabled
```

```r
print(result)
```

```
## $value
## [1] -1.991
## 
## $par
## [1]  0.019235 -0.007651
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
```

```
## Warning: Package lhs not installed, lhs_init disabled
```

```r
print(result)
```

```
## $value
## [1] -1.998
## 
## $par
## [1] -0.008495  0.004030
## 
## $function_calls
## [1] 500
## 
## $break_flag
## [1] "max number of function calls reached"
```


The arguments take in a function that is evaluate at x, and a list of parameters.

