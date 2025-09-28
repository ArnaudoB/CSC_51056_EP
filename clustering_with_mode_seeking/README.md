# Hill Climbing Lab

## Goal
The goal of this lab is to implement a hill climbing clustering algorithm with density estimation and persistence analysis in C++.

## Implemented
- Reading 2D point data from a file (`.xy` format).
- Construction of the k-NN graph (naive quadratic implementation).
- Density estimation as the inverse RMS distance to k neighbors.
- Hill climbing forest construction and label assignment.
- Persistence computation with ToMATo rule (`compute_persistence_tau`).

## To do
- Enable higher dimensional points
- Add **visualization** of the clusters and persistence values.
- Optimize the **k-NN search** using a kd-tree (currently quadratic).
