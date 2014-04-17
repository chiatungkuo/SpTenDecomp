## Nonnegative Tensor Decomposition with Direct Factor-level Sparsity Control
This directory contains MATLAB code for nonnegative tensor CP (Parafac) decomposition with direct control over factor level sparsity.
This code builds upon and utilizes routines from existing public software packages as listed below.
The code runs on MATLAB version 6.0 or higher.

File(s): 

+ cpNonnegSp.m: decompose a tensor into nonnegative factors with specified sparseness on each mode 
+ tuckerNonneg: decompose a tensor into Tucker model with nonnegative core and factors

Dependence: 

+ [Tensor toolbox](http://www.sandia.gov/~tgkolda/TensorToolbox/) (for the tensor constructs in general)
+ [N-way toolbox](http://www.models.life.ku.dk/nwaytoolbox/) (for fast NNLS solver)
+ [nmfpack](http://www.cs.helsinki.fi/u/phoyer/software.html) (for the sparse projection implementation)

This method is described in the paper [Directed Interpretable Discovery in Tensors with Sparse Projection](http://kuo.idav.ucdavis.edu/pubs/sdm2014).
The data sets used in the experiments in the paper above consist of [ORL face images from AT&T Laboratories Cambridge](http://www.cl.cam.ac.uk/research/dtg/attarchive/facedatabase.html) and another private fMRI scan data set. 
If you use the software please cite the reference above. 
