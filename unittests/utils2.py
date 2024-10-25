"""
This file contains some helper functions for testing
Provided by CS 61C staff
"""

import numc as nc
import numpy as np
import hashlib, struct
from typing import Union, List
import operator
import time

"""
Global vars
"""
num_samples = 1000
decimal_places = 5
func_mapping = {
    "add": operator.add,
    "sub": operator.sub,
    "mul": operator.mul,
    "neg": operator.neg,
    "abs": operator.abs,
    "pow": operator.pow
}

func_mapping_np = {
    "add": operator.add,
    "sub": operator.sub,
    "mul": np.matmul,
    "neg": operator.neg,
    "abs": operator.abs,
    "pow": np.linalg.matrix_power
}
# """
# Returns a dumbpy matrix and a numc matrix with the same data
# """
# def dp_nc_matrix(*args, **kwargs):
#     if len(kwargs) > 0:
#         return dp.Matrix(*args, **kwargs), nc.Matrix(*args, **kwargs)
#     else:
#         return dp.Matrix(*args), nc.Matrix(*args)

"""
Returns a random dumbpy matrix and a random numc matrix with the same data
seed, low, and high are optional
"""
def rand_np_nc_matrix(rows, cols, seed=0, low=0.0, high=1.0):
    mat_nc = nc.Matrix(rows, cols, rand=True, seed=seed, low=low, high=high)
    mat_np = np.array(nc.to_list(mat_nc))
    return mat_np, mat_nc

"""
Element-wise compare of np and nc matrix, within certain decimal places
"""
def cmp_np_nc_matrix(np_mat: np.ndarray, nc_mat: nc.Matrix, precision=decimal_places):
    if np_mat.shape != nc_mat.shape:
        print("matrix shape diff")
        print("np:", np_mat.shape)
        print("nc:", nc_mat.shape)
        return False
    if len(nc_mat.shape) == 1:
        # 1D matrix
        for i in range(nc_mat.shape[0]):
            if round(float(np_mat[i]), precision) != round(float(nc_mat[i]), precision):
                print("matrix value diff")
                print(i)
                print("np:", np_mat[i])
                print("nc:", nc_mat[i])
                return False
    else:
        # 2D matrix
        for i in range(nc_mat.shape[0]):
            for j in range(nc_mat.shape[1]):
                if round(float(np_mat[i][j]), precision) != round(float(nc_mat[i][j]), precision):
                    print("matrix value diff")
                    print(i, j)
                    print("np:", np_mat[i][j])
                    print("nc:", nc_mat[i][j])
                    return False     
    return True

# """
# Returns whether the given dumbpy matrix dp_mat is equal to the numc matrix nc_mat
# This function allows a reasonable margin of( floating point errors
# """
# def cmp_np_nc_matrix(dp_mat: dp.Matrix, nc_mat: nc.Matrix):
#     return rand_md5(dp_mat) == rand_md5(nc_mat)

"""
Test if numc returns the correct result given an operation and some matrices.
If speed_up is set to True, returns the speedup as well
"""
def compute(np_mat_lst: List[Union[np.ndarray, int]],
    nc_mat_lst: List[Union[nc.Matrix, int]], op: str, precision=decimal_places):
    f = func_mapping[op]
    fnp = func_mapping_np[op]
    nc_start, nc_end, np_start, np_end = None, None, None, None
    nc_result, np_result = None, None
    assert(op in list(func_mapping.keys()))
    assert(op in list(func_mapping_np.keys()))
    if op == "neg" or op == "abs":
        assert(len(np_mat_lst) == 1)
        assert(len(nc_mat_lst) == 1)
        nc_start = time.time()
        nc_result = f(nc_mat_lst[0])
        nc_end = time.time()

        np_start = time.time()
        np_result = fnp(np_mat_lst[0])
        np_end = time.time()
    else:
        assert(len(np_mat_lst) > 1)
        assert(len(nc_mat_lst) > 1)
        nc_start = time.time()
        nc_result = nc_mat_lst[0]
        for mat in nc_mat_lst[1:]:
            nc_result = f(nc_result, mat)
        nc_end = time.time()

        np_start = time.time()
        np_result = np_mat_lst[0]
        for mat in np_mat_lst[1:]:
            np_result = fnp(np_result, mat)
        np_end = time.time()
    # Check for correctness
    is_correct = cmp_np_nc_matrix(np_result, nc_result, precision=precision)
    return is_correct, (np_end - np_start) / (nc_end - nc_start)

"""
Print speedup
"""
def print_speedup(speed_up):
    print("Speed up is:", speed_up)

# """
# Generate a md5 hash by sampling random elements in nc_mat
# """
# def rand_md5(mat: Union[dp.Matrix, nc.Matrix]):
#     np.random.seed(1)
#     m = hashlib.md5()
#     if len(mat.shape) > 1:
#         rows, cols = mat.shape
#         total_cnt = mat.shape[0] * mat.shape[1]
#         if total_cnt < num_samples:
#             for i in range(rows):
#                 for j in range(cols):
#                     m.update(struct.pack("f", round(mat[i][j], decimal_places)))
#         else:
#             for _ in range(num_samples):
#                 i = np.random.randint(rows)
#                 j = np.random.randint(cols)
#                 m.update(struct.pack("f", round(mat[i][j], decimal_places)))
#     else:
#         total_cnt = mat.shape[0]
#         if total_cnt < num_samples:
#             for i in range(total_cnt):
#                 m.update(struct.pack("f", round(mat[i], decimal_places)))
#         else:
#             for _ in range(num_samples):
#                 i = np.random.randint(total_cnt)
#                 m.update(struct.pack("f", round(mat[i], decimal_places)))
#     return m.digest()
