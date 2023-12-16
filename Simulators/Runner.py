from MatRaptorSim import run_matraptor
from ExTensorSim import run_extensor
from SparseTPU import run_sparsetpu
from SpArch import run_sparch
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
sparse_matrix_2 = generate_random_sparse_matrix(1000, 1000, 0.01)
sparse_matrix_3 = generate_random_sparse_matrix(1000, 500, 0.01)
sparse_matrix_4 = generate_random_sparse_matrix(500, 1000, 0.01)

sparse_matrix_5 = generate_random_sparse_matrix(10000, 10000, 0.0001)
sparse_matrix_6 = generate_random_sparse_matrix(10000, 10000, 0.0001)

sparse_matrix_9 = generate_random_sparse_matrix(1000, 1000, 0.001)

gen = np.random.default_rng()
vec1 = gen.integers(low=1, high = 10, size = (1000,1))
vec2 = gen.integers(low=1, high = 10, size = (500,1))
vec3 = gen.integers(low=1, high = 10, size = (1000,1))


sparse_matrix_10 = generate_random_sparse_matrix(128,10000)
vec4 = gen.integers(low=1, high=10, size=(10000,1) )



'''

Note for anybody who ends up here:

The route matrix and routehi, routeb, routec vectors do NOT work on sparse-tpu. For some reason, a certain weight just blows up and instantly hits the integer limit.
I suspect it has something to do with the fact that they are floating-point values, as integer-valued workloads work completely fine. 
print("SpMV - Real World Datasets:")

route = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route.mtx')
routehi = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_hi.mtx')
routeb = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_b.mtx')
routec = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_c.mtx')

print("\n---------SM1 (1000 x 1000) @0.01 & V1(1000)---------")
print("\nMatRaptor")
run_matraptor(sparse_matrix_1, vec1)

print("\nExTensor")
run_extensor(sparse_matrix_1,vec1,100,50)

print("\nSpArch")
run_sparch(sparse_matrix_1,vec1)

print("\nSparse-TPU")
run_sparsetpu(sparse_matrix_1,vec1)

print("\n---------SM3 (1000 x 500) @ 0.01 & V2 (500)---------")

print("\nMatRaptor")
run_matraptor(sparse_matrix_3,vec2)

print("\nExTensor")
run_extensor(sparse_matrix_3,vec2,100,50)

print("\nSpArch")
run_sparch(sparse_matrix_3,vec2)

print("\nSparse-TPU")
run_sparsetpu(sparse_matrix_3,vec2)

print("\n---------SM9 (1000 x 1000) @ 0.001 & V3(1000)---------")

print("\nMatRaptor")
run_matraptor(sparse_matrix_9,vec3)

print("\nExTensor")
run_extensor(sparse_matrix_9,vec3,100,50)

print("\nSpArch")
run_sparch(sparse_matrix_9,vec3)

print("\nSparse-TPU")
run_sparsetpu(sparse_matrix_9,vec3)

print("\n\n SPGEMM - Random Sparse Matrices:")

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 1 (1000,1000,0.01)")

print("\n---------MatRaptor---------")
run_matraptor(sparse_matrix_1,sparse_matrix_1)

print("\n---------ExTensor---------")
run_extensor(sparse_matrix_1,sparse_matrix_1,100,50)

print("\n---------SpArch---------")
run_sparch(sparse_matrix_1,sparse_matrix_1)

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 2 (1000,1000,0.1)")

print("\n---------MatRaptor---------")
run_matraptor(sparse_matrix_1,sparse_matrix_2)

print("\n---------ExTensor---------")
run_extensor(sparse_matrix_1,sparse_matrix_2,100,50)

print("\n---------SpArch---------")
run_sparch(sparse_matrix_1,sparse_matrix_2)

print("\n\nsparse matrix 1 (1000,1000,0.01) x sparse matrix 3 (1000,500,0.01)")

print("\n---------MatRaptor---------")
run_matraptor(sparse_matrix_1,sparse_matrix_3)

print("\n---------ExTensor---------")
run_extensor(sparse_matrix_1,sparse_matrix_3,100,50)

print("\n---------SpArch---------")
run_sparch(sparse_matrix_1,sparse_matrix_3)

print("\n\nsparse matrix 3 (1000,500,0.01) x sparse matrix 4 (500,1000,0.01)")

print("\n---------MatRaptor---------")
run_matraptor(sparse_matrix_3,sparse_matrix_4)

print("\n---------ExTensor---------")
run_extensor(sparse_matrix_3,sparse_matrix_4,100,50)

print("\n---------SpArch---------")
run_sparch(sparse_matrix_3,sparse_matrix_4)

print("\n\nsparse matrix 5 (10000,10000,0.0001) x sparse matrix 6 (10000,10000,0.0001)")

print("\n---------MatRaptor---------")
run_matraptor(sparse_matrix_5,sparse_matrix_6)

print("\n---------ExTensor---------")
run_extensor(sparse_matrix_5,sparse_matrix_6,1000,250)

print("\n---------SpArch---------")
run_sparch(sparse_matrix_5,sparse_matrix_6)


print("SPGEMM - Random Real-World Matrices:")

mbeacxc = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
p2pgnutella = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\p2p-Gnutella31.mtx')
pesa = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\pesa.mtx')
route = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route.mtx')
routehi = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_hi.mtx')
routeb = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_b.mtx')
routec = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_c.mtx')

print("\nmbeacxc")
print("\n---------MatRaptor---------")
run_matraptor(mbeacxc, mbeacxc)

print("\n---------ExTensor---------")
run_extensor(mbeacxc,mbeacxc,30,10)

print("\n---------SpArch---------")
run_sparch(mbeacxc,mbeacxc)

print("\ngnutella")
print("\n---------MatRaptor---------")
run_matraptor(p2pgnutella,p2pgnutella)

print("\n---------ExTensor---------")
run_extensor(p2pgnutella,p2pgnutella,30,10)

print("\n---------SpArch---------")
run_sparch(p2pgnutella,p2pgnutella)


print("\npesa")
print("\n---------MatRaptor---------")
run_matraptor(pesa,pesa.T)

print("\n---------ExTensor---------")
run_extensor(pesa,pesa.T,30,10)

print("\n---------SpArch---------")
run_sparch(pesa,pesa.T)

'''
