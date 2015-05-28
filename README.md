fmmtl: FMM Template Library
=====

fmmtl is a structured dense matrix algorithms library which aids in the research, development, and use of advanced methods for systems of equations of the form:

![equation](http://latex.codecogs.com/gif.latex?r_i%3D%5Csum_jK%28t_i%2Cs_j%29%5C%2Cc_j)<br/>
where<br/>
![equation](http://latex.codecogs.com/gif.latex?K) is the _kernel_ generating the elements of the matrix,<br/>
![equation](http://latex.codecogs.com/gif.latex?s_j) are the _sources_ of the kernel,<br/>
![equation](http://latex.codecogs.com/gif.latex?c_j) are the _charges_ of the sources,<br/>
![equation](http://latex.codecogs.com/gif.latex?t_i) are the _targets_ of the kernel (which may be equivalent to the sources),<br/>
![equation](http://latex.codecogs.com/gif.latex?r_i) are the _results_.<br/>

This is a kernel-matrix equation. Matrices of this form can be found in a wide variety of fields include physics, statistics, and machine learning. Methods for accelerating matrix-vector products and direct solutions of systems of this form often take advantage of a (often heirarchically) low-rank representation of the kernel, K.

This library intends to collect kernels and their low-rank expansions and provide fast, abstracted algorithms for working with any of them.

Primary Authors:
* Cris Cecka (ccecka@seas.harvard.edu)
* Simon Layton

Contributors:
* Christopher Cooper
* Aparna Chandramowlishwaran
* Esmail Fadae, Brian Bresnahan

Dependencies:
* [g++ 4.7](http://gcc.gnu.org/) or higher.
* [Thrust Header Library](https://github.com/thrust/thrust). (Tested with version 1.8)
* Experimental: [FLENS Header Library](https://github.com/ccecka/FLENS). (Extended with CUBLAS/CUSOLVER support)

Optional:
* [CUDA 5.5+](https://developer.nvidia.com/cuda-downloads) for GPU acceleration.


Building:
* Edit Makefile.inc to edit paths
* Enter unit_tests/ or examples/
* 'make'
