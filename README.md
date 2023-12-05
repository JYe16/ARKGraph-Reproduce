# COP-5725: Database Systems - Fall 2023

## Term Project: Jackie Ye (ZY22B), Yichao Wang (YW20BP)

### Implementation of paper: ARKGraph: All-Range Approximate K-Nearest-Neighbor Graph

### Paper Link: https://dl.acm.org/doi/abs/10.14778/3603581.3603601

***

## Environment Requirement: Clion 2023.2.2

### Please open the project's root directory in Clion and then click 'Run'

***

## Datasets

### All datasets raw files should be placed under ***cmake-build-debug/data/***

### DEEP10M: Click [here](https://fsu-my.sharepoint.com/:u:/g/personal/zy22b_fsu_edu/Eekl0KQuG-pNh66Sa7SWmhoBJMYzO73MI8Q28GXcNcPQQg?e=OOgEs0) to download.

### BigGraph: Click [here](https://fsu-my.sharepoint.com/:u:/g/personal/zy22b_fsu_edu/EV9Vpk2aYzpEv7wZ7K3JlJwBfyE5waEhvWFIVRVtWZLUBQ?e=hRpVLh) to download.

***

## Additional Information

### Since the original approach has leveraged parallels computing (OpenMP) to increase the efficiency, they are able to evaluate the approach using the whole dataset. However, due to the limitation of computing resources, it would be very likely to take a prohibitively long time to index the whole dataset on a commercial PC. Therefore, we recommend to test our code by indexing less than 100 of N.
