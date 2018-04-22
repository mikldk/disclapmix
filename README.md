# disclapmix

An R package (<https://www.r-project.org/>) to perform inference in a mixture of discrete Laplace distributions using the EM algorithm.
Intended for forensic Y chromosomal STR (Y-STR) haplotype analyses. 

See documentation included in package (vignettes and manual) at <https://mikldk.github.io/disclapmix/>.


## Installation

You first need `R` (<https://www.r-project.org/>). 
Then you can install `disclapmix` from CRAN using

```r
install.packages("disclapmix")
```

You can also install the development version from GitHub by using the `remotes` package (<https://CRAN.R-project.org/package=remotes>):

```r
# install.packages("remotes")
remotes::install_github("mikldk/disclapmix")
```

## Getting started

Refer to the included vignettes. You can get an overview of the included vignettes by the following `R` command:

```r
vignette(package = "disclapmix")
```

To read a vignette, type:

```r
vignette("introduction", package = "disclapmix")
```

## References

Andersen MM, PS Eriksen, N Morling (2013). 
*The discrete Laplace exponential family and estimation of Y-STR haplotype frequencies*.
Journal of Theoretical Biology 329. 
<https://doi.org/10.1016/j.jtbi.2013.03.009>

Andersen MM, PS Eriksen, N Morling (2013). 
*A gentle introduction to the discrete Laplace method for estimating Y-STR haplotype frequencies*.
arXiv:1304.2129. 
<https://arxiv.org/abs/1304.2129>

## Disclaimer

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

## License

License: GPL (>= 2).

## Build status

Travis CI:

[![Travis-CI Build Status](https://travis-ci.org/mikldk/disclapmix.svg?branch=master)](https://travis-ci.org/mikldk/disclapmix)

