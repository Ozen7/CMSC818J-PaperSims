# CMSC818J-PaperSims
A set of Python simulators meant to be used in a paper written for my class "CMSC 818J"

# How To Run:
Each .py file has a "run_*****" method that can be imported to run that particular simulator. The only things you have to provide are the two matrices to be multiplied (matrix and 1-d vector, in Sparse-TPU's case). Make sure that the two matrices are in a scipy compressed format (Scipy CSC, Scipy CSR, etc.). For an example, take a look at the Runner.py file.

# Libraries needed:
- Scipy
- Numpy
- bisect (used in an optimization for ExTensor sim)
