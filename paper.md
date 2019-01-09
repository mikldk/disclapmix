---
title: 'Discrete Laplace mixture model with applications in forensic genetics'
tags:
  - Y-chromosome
  - forensic genetics
  - population genetics
  - weight of evidence
authors:
 - name: Mikkel Meyer Andersen
   orcid: 0000-0002-0234-0266
   affiliation: 1
affiliations:
 - name: Department of Mathematical Sciences, Aalborg University, Denmark
   index: 1
date: 23 April 2018
bibliography: paper.bib
---

# Summary

This R package implements a model based on a mixture of multivariate discrete Laplace distributions that has applications in forensic genetics. The implementation consists of parameter estimation and various functionalities. 
The method and this package were (and still are) used by multiple groups for e.g. frequency estimation  [@AndersenDisclap2013; @YHRD15; @Willuweit2018; @Roewer2018; @Familias; @CeredaDIP; @CeredaRare], cluster analysis [@AndersenDisclapCluster2014], and mixture interpretation [@AndersenDisclapMixture2015; Taylor2018]. 
Below, background for the method and package is described.

Estimating haplotype frequencies is important in e.g. forensic genetics, where the frequencies are used to calculate the likelihood ratio for the evidential weight of a DNA profile found at a crime scene [@AndersenDisclap2013; @DJBwoe2]. 
Estimation is naturally based on a population model, motivating the investigation of the Fisher-Wright model [@Fisher1930; @Wright1931; @Ewens1972; @OhtaKimura1973] of evolution for haploid lineage DNA markers.

An exponential family (a class of probability distributions that is well understood in probability theory) called the 'discrete Laplace distribution' was described in [@AndersenDisclap2013] that also illustrates how well the discrete Laplace distribution approximates a more complicated distribution that arises by investigating the well-known population genetic Fisher-Wright model of evolution by a single-step mutation process [@Caliebe2010].

In [@AndersenDisclap2013], it was also shown that the discrete Laplace distribution can be used to estimate haplotype frequencies for haploid lineage DNA markers (such as Y-chromosomal short tandem repeats), which in turn can be used to assess the evidential weight of a DNA profile found at a crime scene. This was done by making inference in a mixture of multivariate discrete Laplace distributions using the EM algorithm to estimate the probabilities of membership of a set of unobserved subpopulations and also estimate the central haplotypes of the subpopulations. The implementation was made efficient by avoiding to construct the model matrix explicitely which made it possible to perform the calculations on a normal computer.

This package implements the method described in [@AndersenDisclap2013] as a freely available open source software R package using both R [@R] and C++ [@Rcpp].
The documentation of `disclapmix` consists of manual pages for the various available functions, articles describing how to perform contiguous analyses (*vignettes*), and unit tests.

I would like to thank Poul Svante Eriksen (Aalborg University, Denmark) for tips, hints, helpful discussions and help with implementation and debugging.

# References
