# COP-5725: Database Systems - Fall 2023

## Term Project: Jackie Ye (ZY22B), Yichao Wang ()

### Implementation of paper: ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph

### Paper Link: https://dl.acm.org/doi/abs/10.14778/3603581.3603601

***

## Environment Requirement: Clion 2023.2.2

### Please open the project's root directory in Clion and then click 'Run'

***

## Datasets

### All datasets raw files should be placed under ***cmake-build-debug/data/***

### DEEP10M:

### BigGraph:

***

## Additional Information

### Since the original approach has leveraged parallels computing (OpenMP) to increase the efficiency, they are able to evaluate the approach using the whole dataset. However, due to the limitation of computing resources, it would be very likely to take a prohibitively long time to index the whole dataset on a commercial PC. Therefore, we recommend to test our code by indexing less than 1000 of N.