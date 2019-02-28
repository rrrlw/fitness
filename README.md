# fitness

[![Travis-CI Build Status](https://travis-ci.org/rrrlw/fitness.svg?branch=master)](https://travis-ci.org/rrrlw/fitness)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/rrrlw/fitness?branch=master&svg=true)](https://ci.appveyor.com/project/rrrlw/fitness)
[![Coverage Status](https://img.shields.io/codecov/c/github/rrrlw/fitness/master.svg)](https://codecov.io/github/rrrlw/fitness?branch=master)

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN version](http://www.r-pkg.org/badges/version/fitness)](https://CRAN.R-project.org/package=fitness)
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/fitness)](https://CRAN.R-project.org/package=fitness)

## Overview

The fitness package permits working with fitness landscapes in R.

## Installation

The development version of fitness can be installed using the devtools package.

```r
# install from GitHub
devtools::install_github("rrrlw/fitness")
```

## Sample code

## Functionality

The following models of fitness landscapes are implemented:

1. NK model
1. Rough Mt. Fuji model

In addition, the following functionality for generic fitness landscapes
is implemented:

* Item 1
* Item 2

## Contribute

To contribute to fitness, you can create issues for any bugs/suggestions on the [issues page](https://github.com/rrrlw/fitness/issues). You can also fork the fitness repository and create pull requests to add features you think will be useful for users.

## Plan

read_* functions to read in each type of model (or general fitness landscape?)
write_* functions to write each type of model (or general fitness landscape?)
method to find optimal peak, to find all peaks (same for valleys?)
functions to generate NK landscapes, RMF landscapes
utility functions: calculate lattice distance, list all neighbors (values or coordinates?) for given spot
