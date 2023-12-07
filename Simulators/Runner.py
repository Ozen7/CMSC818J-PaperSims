from MatRaptorSim import run_matraptor
from ExTensorSim import run_extensor
import numpy as np
import scipy.io as sio
from scipy.sparse import csc_array
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix
from scipy.sparse import random

def generate_random_sparse_matrix(rows, cols, density):
    """
    Generate a random sparse matrix with the given number of rows, columns, and density.

    Parameters:
    - rows: Number of rows in the matrix.
    - cols: Number of columns in the matrix.
    - density: Density of non-zero elements in the matrix (between 0 and 1).

    Returns:
    - A random sparse matrix.
    """
    # Ensure density is within the valid range
    density = max(0.0, min(1.0, density))

    # Generate a random sparse matrix
    sparse_matrix = random(rows, cols, density=density, format='csr')

    return sparse_matrix


sparse_matrix_1 = generate_random_sparse_matrix(1000, 1000, 0.01)   
sparse_matrix_2 = generate_random_sparse_matrix(1000, 1000, 0.1)
sparse_matrix_3 = generate_random_sparse_matrix(1000, 500, 0.1)
sparse_matrix_4 = generate_random_sparse_matrix(500, 1000, 1)



print("Random Sparse Matrices:")

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 1 (1000,1000,0.01)")

print("MatRaptor")
run_matraptor(sparse_matrix_1,sparse_matrix_1)

print("ExTensor")
run_extensor(sparse_matrix_1,sparse_matrix_1,30,10)

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 2 (1000,1000,0.1)")

print("MatRaptor")
run_matraptor(sparse_matrix_1,sparse_matrix_2)

print("ExTensor")
run_extensor(sparse_matrix_1,sparse_matrix_2,30,10)

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 3 (1000,500,1)")

print("MatRaptor")
run_matraptor(sparse_matrix_1,sparse_matrix_3)

print("ExTensor")
run_extensor(sparse_matrix_1,sparse_matrix_3,30,10)

print("\n\nsparse matrix 2 (1000,1000,0.1) x sparse matrix 3 (1000,500,1)")

print("MatRaptor")
run_matraptor(sparse_matrix_2,sparse_matrix_3)

print("ExTensor")
run_extensor(sparse_matrix_2,sparse_matrix_3,30,10)

print("\n\nsparse matrix 3 (1000,500,1) x sparse matrix 4 (500,1000,10)")

print("MatRaptor")
run_matraptor(sparse_matrix_3,sparse_matrix_4)

print("ExTensor")
run_extensor(sparse_matrix_3,sparse_matrix_4,30,10)



mbeacxc = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
p2pgnutella = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\p2p-Gnutella31.mtx')
pesa = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\pesa.mtx')
route = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route.mtx')
routehi = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_hi.mtx')
routeb = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_b.mtx')
routec = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_c.mtx')

print("----mbeacxc----")
print("MatRaptor")
run_matraptor(mbeacxc, mbeacxc)

print("ExTensor")
run_extensor(mbeacxc,mbeacxc,30,10)

print("----gnutella----")
print("MatRaptor")
run_matraptor(p2pgnutella,p2pgnutella)

print("ExTensor")
run_extensor(p2pgnutella,p2pgnutella,30,10)


print("----pesa----")
print("MatRaptor")
run_matraptor(pesa,pesa.T)

print("ExTensor")
run_extensor(pesa,pesa.T,30,10)


# SPMV
#print("routehi")
#run_matraptor(routehi, route)

#print("routeb")
#run_matraptor(routeb, route)