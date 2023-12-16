from typing import Callable
import numpy as np
import scipy.io as sio
from scipy.sparse import csc_array
from scipy.sparse import csr_matrix
from scipy.sparse import coo_matrix

'''
mbeacxc = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\mbeacxc.mtx')
p2pgnutella = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\p2p-Gnutella31.mtx')
pesa = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets\pesa.mtx')
route = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route.mtx')
routepath = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_hi.mtx')
routeb = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_b.mtx')
routec = sio.mmread('C:\Workspace\CMSC818J\PaperSims\CMSC818J-PaperSims\Simulators\Datasets/route_c.mtx')
print(mbeacxc.shape)
print(p2pgnutella.shape)
print(pesa.shape)
print(route.shape)
print(routepath.shape)
print(routeb.shape)
print(routec.shape)
'''

data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
indices = [0, 2, 4, 1, 3, 4, 0, 2, 3]
indptr = [0, 2, 4, 7, 9]
datacsr = csr_matrix((data,indices,indptr))
print(datacsr.data)
print(datacsr.indices)
print(datacsr.indptr)

"""sumary_line


def csr_to_c2sr(data, indices, indptr, num_channels):
    # Create arrays to store C2SR format data
    c2sr_values = [[] for x in range(num_channels)]
    c2sr_column_ids = [[] for x in range(num_channels)]
    c2sr_row_length = [[] for x in range(num_channels)]
    c2sr_row_pointer = [[] for x in range(num_channels)]

    # Loop through each row in CSR format
    for i in range(len(indptr) - 1):
        # Get the start and end indices for the current row in CSR format
        start_idx = indptr[i]
        end_idx = indptr[i + 1]

        # Get the channel assigned to the current row
        channel = i % num_channels
        
        # Append the row length and row pointer for the current row in C2SR format
        c2sr_row_length[channel].append(end_idx - start_idx)
        c2sr_row_pointer[channel].append(len(c2sr_values[channel]))

        # Loop through the non-zero elements of the current row
        for j in range(start_idx, end_idx):
            # Append the value and column id for the current element in C2SR format
            c2sr_values[channel].append(data[j])
            c2sr_column_ids[channel].append(indices[j])
        
        
        

    return c2sr_values, c2sr_column_ids, c2sr_row_length, c2sr_row_pointer

# Example usage
# data = [1, 2, 3, 4, 5, 6, 7, 8, 9]
# indices = [0, 2, 4, 1, 3, 4, 0, 2, 3]
# indptr = [0, 2, 4, 7, 9]
datacsr = csr_matrix(sparsematrix)
data = datacsr.data
indices = datacsr.indices
indptr = datacsr.indptr

num_channels = 3

c2sr_values, c2sr_column_ids, c2sr_row_length, c2sr_row_pointer = csr_to_c2sr(data, indices, indptr, num_channels)

# Print the results
assert(len(c2sr_values) == len(c2sr_column_ids))

print("C2SR Values:", c2sr_values)
print("C2SR Column IDs:", c2sr_column_ids)
print("C2SR Row Lengths:", c2sr_row_length)
print("C2SR Row Pointers:", c2sr_row_pointer)


x = np.array([1,2,3])
y = np.array([[1,2,3],[3,2,1],[4,5,6]])
print(np.matmul(x,y))
print(np.matmul(y.T,x.T))

x = np.array([[1,2,3,4],[1,2,3,4],[1,2,3,4],[1,2,3,4]])
y = np.array([[1,2,3,4],[3,2,1,4],[4,5,6,4],[1,2,3,4]])

C = np.zeros((4,4))
for chunk_start_A in range(0, 4, 2):
    for chunk_start_B in range(0, 4, 2):
        # Extract chunks from the single column of A and single row of B
        A_chunk = x[chunk_start_A:chunk_start_A+2, :]
        B_chunk = y[:, chunk_start_B:chunk_start_B+2]
        
        print(A_chunk)
        print(B_chunk)
        print(C)
        print(A_chunk @ B_chunk)

        # Multiply and accumulate
        C[chunk_start_A:chunk_start_A+2, chunk_start_B:chunk_start_B+2] += A_chunk @ B_chunk

print(C)
print(np.matmul(x,y))


"""