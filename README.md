## LCMAP Change Detection

This project contains application source code for Change Detection C library
and related scripts.

## Dependencies

### CentOS

TBD

### Ubuntu

On Ubuntu, you will need the following packages installed in order to compile the
CCDC C library:

```bash
$ sudo apt-get install build-essential
$ sudo apt-get install libgsl0-dev libgsl0ldbl gsl-bin
$ sudo apt-get install r-base r-recommended
```

Currently, the following R library is also required:

```bash
$ sudo R
```
```r
> install.packages("glmnet", repos = "http://cran.us.r-project.org")
```

The classification tool requires the additional installs:

```bash
sudo apt-get install libmatio-dev libmatio2
sudo apt-get install gfortran
```


## Installation

To install, simply run the top-level ``make`` target:

```bash
$ make
```

The executables and scripts will be installed into ``./bin`` by default. This
can be overridden by setting a ``BIN`` environment variable or using a ``BIN``
variable when running the target:

```bash
$ BIN=/my/path/bin make
```

Similarly, you may override the include and lib paths, should they be different
on your system:
 * ``XML2INC``
 * ``ESPANIC``
 * ``GSL_SCI_INC``
 * ``GSL_SCI_LIB``
 * ``FORTRAN``
 * ``HDF5INC``
 * ``HDF5LIB``
 * ``MATIOLIB``

## Usage

TBD

## Implementation

### CCDC - Continuous Change Detection and Classification (Algorithm)

* <b>NOTE:</b> This algorithm is not validated and considered prototype.
* See [CCDC ADD](http://landsat.usgs.gov/documents/ccdc_add.pdf) for the detailed description.

## More Information

This project is hosted by the US Geological Survey (USGS) Earth Resources Observation and
Science (EROS) Land Change Monitoring, Assessment, and Projection (LCMAP) Project.
For questions regarding this source code, please contact the Landsat Contact Us page and
specify USGS LCMAP in the "Regarding" section. https://landsat.usgs.gov/contactus.php
