# gSegMulti

This repository provides the R functions for the paper:
[Graph-based multiple change-point detection](https://arxiv.org/abs/2110.01170)

**Authors:** Yuxuan Zhang and Hao Chen from University of California, Davis.

## Paper Overview
This paper propose a new multiple change-point detection framework for multivariate and non-Euclidean data. First, we combine graph-based statistics with wild binary segmentation or seeded binary segmentation to search for a pool of candidate change-points. 
We then prune the candidate change-points through a novel goodness-of-fit statistic.
Numerical studies show that this new framework outperforms existing methods under a wide range of settings. The resulting change-points can further be arranged hierarchically based on the goodness-of-fit statistic. The new framework is illustrated on a Neuropixels recording of an awake mouse.

## R Functions
In the R folder, you can find all main functions from the paper:

- `gsegwbs.R`
- `gsegsbs.R`
- `gbe.R`
