# stls_sdp

This is a [CVX](http://cvxr.com/cvx/) code implementing a semidefinite programming (SDP) relaxation for Structured Total Least Squares (STLS).
The STLS problem is the following:
given an affine space of matrices *L*, and a matrix *θ* in *L*,
find the closest rank deficient matrix to *θ* in *L*.
Some common choices for the affine space *L* include the spaces of Hankel, Toeplitz and Sylvester matrices.
The SDP relaxation implemented here is introduced in the paper:

- [Cifuentes (2019)](http://www.mit.edu/~diegcif/), "A convex relaxation for structured total least squares", [arXiv:1904.09661](https://arxiv.org/abs/1904.09661)
