---
title: "Maximum Entropy: Claims Vignette"
author: "Tucker McElroy"
date: '2024-09-10'
output: html_document
---



# Claims

We examine the "National Unemployment Insurance Weekly Claims" data from the
[Department of Labor](https://oui.doleta.gov/unemploy/claims.asp), referred to
as *claims* for short. The data covers a span from January 7, 1967 through 
April 22, 2023.  
  
## Load Code

Load up the libraries needed.


```r
library(seasonal)
```

```
## Warning: package 'seasonal' was built under R version 4.1.2
```

```r
library(devtools)
```

```
## Warning: package 'devtools' was built under R version 4.1.1
```

```
## Loading required package: usethis
```

```
## Warning: package 'usethis' was built under R version 4.1.1
```

```r
library(Rcpp)
```

Load up code from Maximum Entropy repository.






