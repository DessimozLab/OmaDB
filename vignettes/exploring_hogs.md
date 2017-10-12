---
title: "Exploring Hierarchical orthologous groups"
author: "Klara Kaleb"
date: "2017-10-12"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Exploring Hierarchical orthologous groups}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This little vignette shows you how to get started with the data available for HOGs in the `roma` package.

## The HOGs

Hierarchical orthologous groups (also known as HOGs) are sets of genes that are defined with respect to particular taxonomic ranges of interest[1]. They group genes that have descended from a single common ancestral genes in that taxonomic range. 

HOGs hold a lot of useful information and have many applications in various contexts, including inference of gene function, study of gene evolution dynamics and comparative genomics. Each HOG has a taxonomic range - within it, a given HOG can branch into constructs known as subHOGs which arise in an event of gene duplication.

HOGs can be retrived either by their hog id or by one of their members. Let's say we are interested in a gene that goes by the name of HUMAN22168 - below is the way to access its HOG.





