import numc as nc
import random
import time

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""

version = "raw"
repeat = 20

with open(f"timing_{version}.txt", "w") as fout:
    fout.write(f"repeat:{repeat}\n")
    test = "medium add"
    row, col = 1000, 1000 
    print(f"timing {test}, ({row}, {col}):")
    total_time = 0.0
    result = nc.Matrix(row, col, rand=True, seed=0, low=-1)
    for i in range(repeat):
        mat = nc.Matrix(row, col, rand=True, seed=i+1, low=-1)
        start = time.time()
        result += mat
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{col}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "medium sub"
    row, col = 1000, 1000
    print(f"timing {test}, ({row}, {col}):")
    total_time = 0.0
    result = nc.Matrix(row, col, rand=True, seed=0, low=-1)
    for i in range(repeat):
        mat = nc.Matrix(row, col, rand=True, seed=i+1, low=-1)
        start = time.time()
        result -= mat
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{col}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "medium abs"
    row, col = 1000, 1000
    print(f"timing {test}, ({row}, {col}):")
    total_time = 0.0
    for i in range(repeat):
        mat = nc.Matrix(row, col, rand=True, seed=i, low=-1)
        start = time.time()
        abs(mat)
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{col}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "medium neg"
    row, col = 1000, 1000
    print(f"timing {test}, ({row}, {col}):")
    total_time = 0.0
    for i in range(repeat):
        mat = nc.Matrix(row, col, rand=True, seed=i, low=-1)
        start = time.time()
        -mat
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{col}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "medium mul"
    row, col1, col2 = 1001, 999, 1000
    print(f"timing {test}, ({row}, {col1}, {col2}):")
    total_time = 0.0
    for i in range(repeat):
        mat1 = nc.Matrix(row, col1, rand=True, seed=i, low=-1)
        mat2 = nc.Matrix(col1, col2, rand=True, seed=i+1, low=-1)
        start = time.time()
        mat1 * mat2
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{col1},{col2}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "medium pow"
    row = 1000
    print(f"timing {test}, ({row}, {row}):")
    total_time = 0.0
    for i in range(repeat):
        mat = nc.Matrix(row, row, rand=True, seed=i, low=-1)
        start = time.time()
        mat**3
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{row}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "pow one"
    row = 1000
    print(f"timing {test}, ({row}, {row}):")
    total_time = 0.0
    for i in range(repeat):
        mat = nc.Matrix(row, row, rand=True, seed=i, low=-1)
        start = time.time()
        mat**1
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{row}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")

    test = "pow zero"
    row = 1000
    print(f"timing {test}, ({row}, {row}):")
    total_time = 0.0
    for i in range(repeat):
        mat = nc.Matrix(row, row, rand=True, seed=i, low=-1)
        start = time.time()
        mat**0
        total_time += time.time() - start
    total_time = round(total_time, 10)
    print(f"{total_time}s\n")
    line_start = f"{test},({row},{row}):"
    fout.write(f"{line_start:<30}\t{total_time}s\t{total_time/repeat*1000}ms\n")
