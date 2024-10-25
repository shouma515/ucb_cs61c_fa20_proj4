from utils2 import *
from unittest import TestCase

import numc as nc
import random

"""
For each operation, you should write tests to test  on matrices of different sizes.
Hint: use dp_mc_matrix to generate dumbpy and numc matrices with the same data and use
      cmp_dp_nc_matrix to compare the results
"""
class TestBasic(TestCase):
    def test_utils(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(2, 2, seed=0)
        self.assertTrue(cmp_np_nc_matrix(np_mat, nc_mat))

class TestAdd(TestCase):
    def test_small_add(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(2, 2, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(2, 2, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_add(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(234, 789, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(234, 789, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_long_add(self):
        # TODO: YOUR CODE HERE
        mat_pairs = [rand_np_nc_matrix(123, 789, seed=i) for i in range(100)]
        np_mats, nc_mats = zip(*mat_pairs)
        is_correct, speed_up = compute(list(np_mats), list(nc_mats), "add")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_add(self):
    #     # TODO: YOUR CODE HERE
    #     np_mat1, nc_mat1 = rand_np_nc_matrix(5000, 10000, seed=0)
    #     np_mat2, nc_mat2 = rand_np_nc_matrix(5000, 10000, seed=1)
    #     is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "add")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

class TestSub(TestCase):
    def test_small_sub(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(2, 2, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(2, 2, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_sub(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(234, 789, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(234, 789, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_long_sub(self):
        # TODO: YOUR CODE HERE
        mat_pairs = [rand_np_nc_matrix(123, 789, seed=i) for i in range(100)]
        np_mats, nc_mats = zip(*mat_pairs)
        is_correct, speed_up = compute(list(np_mats), list(nc_mats), "sub")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_sub(self):
    #     # TODO: YOUR CODE HERE
    #     np_mat1, nc_mat1 = rand_np_nc_matrix(5000, 5839, seed=0)
    #     np_mat2, nc_mat2 = rand_np_nc_matrix(5000, 5839, seed=1)
    #     is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "sub")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

class TestAbs(TestCase):
    def test_small_abs(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(2, 2, seed=0, low=-1.0)
        is_correct, speed_up = compute([np_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_abs(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(234, 789, seed=0, low=-1.0)
        is_correct, speed_up = compute([np_mat], [nc_mat], "abs")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_abs(self):
    #     # TODO: YOUR CODE HERE
    #     np_mat, nc_mat = rand_np_nc_matrix(5000, 5839, seed=0, low=-1.0)
    #     is_correct, speed_up = compute([np_mat], [nc_mat], "abs")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

class TestNeg(TestCase):
    def test_small_neg(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(2, 2, seed=0, low=-1.0)
        is_correct, speed_up = compute([np_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_neg(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(234, 789, seed=0, low=-1.0)
        is_correct, speed_up = compute([np_mat], [nc_mat], "neg")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_neg(self):
    #     # TODO: YOUR CODE HERE
    #     np_mat, nc_mat = rand_np_nc_matrix(5000, 5839, seed=0, low=-1.0)
    #     is_correct, speed_up = compute([np_mat], [nc_mat], "neg")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

class TestMul(TestCase):
    def test_small_mul(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(2, 2, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(2, 2, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_mul(self):
        # TODO: YOUR CODE HERE
        np_mat1, nc_mat1 = rand_np_nc_matrix(765, 789, seed=0)
        np_mat2, nc_mat2 = rand_np_nc_matrix(789, 567, seed=1)
        is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "mul")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    # def test_long_mul(self):
    #     # TODO: YOUR CODE HERE
    #     random.seed(42)
    #     np_mats, nc_mats = [], []
    #     prev_col = 123
    #     for i in range(10):
    #         # TODO: 1-d array cause problem as nc.to_list(1-d-matrix) will give a list
    #         # and does not provide any info about actual col and row, thus a (n, 1) matrix
    #         # and a (1, n) matrix will result in a same np.ndarray, thus cause dimension
    #         # mismatch problem when multiplying matrix
    #         # how to reproduce: python random seed 42, prev_col 12, range(3)
    #         # TODO: multiply has precision problem when matrix is larger and chain is longer
    #         crnt_col = random.randint(100, 500)
    #         np_mat, nc_mat = rand_np_nc_matrix(prev_col, crnt_col, seed=i, low=-1)
    #         np_mats.append(np_mat)
    #         nc_mats.append(nc_mat)
    #         prev_col = crnt_col
    #     is_correct, speed_up = compute(np_mats, nc_mats, "mul", precision=2)
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)

    # def test_large_mul(self):
    #     # TODO: YOUR CODE HERE
    #     np_mat1, nc_mat1 = rand_np_nc_matrix(5000, 5839, seed=0)
    #     np_mat2, nc_mat2 = rand_np_nc_matrix(5000, 5839, seed=1)
    #     is_correct, speed_up = compute([np_mat1, np_mat2], [nc_mat1, nc_mat2], "mul")
    #     self.assertTrue(is_correct)
    #     print_speedup(speed_up)


class TestPow(TestCase):
    def test_small_pow(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(2, 2, seed=0)
        is_correct, speed_up = compute([np_mat, 3], [nc_mat, 3], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    def test_medium_pow(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(234, 234, seed=0)
        is_correct, speed_up = compute([np_mat, 3], [nc_mat, 3], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_pow_one(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(234, 234, seed=0)
        is_correct, speed_up = compute([np_mat, 1], [nc_mat, 1], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)
    
    def test_pow_zero(self):
        # TODO: YOUR CODE HERE
        np_mat, nc_mat = rand_np_nc_matrix(234, 234, seed=0)
        is_correct, speed_up = compute([np_mat, 0], [nc_mat, 0], "pow")
        self.assertTrue(is_correct)
        print_speedup(speed_up)

    # def test_large_pow(self):
    #     # TODO: YOUR CODE HERE
    #     pass

class TestGet(TestCase):
    def test_get(self):
        # TODO: YOUR CODE HERE
        # 2-d matrix
        a = nc.Matrix([[1,2,3],[4,5,6],[7,8,9]])
        self.assertEqual(a[0:1, 0:1], 1)
        self.assertListEqual(nc.to_list(a[0]), [1, 2, 3])
        self.assertListEqual(nc.to_list(a[0:2]), [[1, 2, 3], [4, 5, 6]])
        self.assertListEqual(nc.to_list(a[0:2, 0:2]), [[1, 2], [4, 5]])
        self.assertListEqual(nc.to_list(a[0:2, 0]), [1,4])
        self.assertListEqual(nc.to_list(a[0, 0:2]), [1,2])
        self.assertEqual(a[0, 0], 1)
        # negative index
        self.assertEqual(a[-3:-2, -3:-2], 1)
        self.assertListEqual(nc.to_list(a[-1]), [7, 8, 9])
        self.assertListEqual(nc.to_list(a[0:-1]), [[1, 2, 3], [4, 5, 6]])
        self.assertListEqual(nc.to_list(a[-2:, 0:-1]), [[4, 5], [7, 8]])
        self.assertListEqual(nc.to_list(a[0:2, -1]), [3,6])
        self.assertListEqual(nc.to_list(a[0, -2:3]), [2,3])
        self.assertEqual(a[-2, -2], 5)
        # 1-d matrix, 1 row
        b = nc.Matrix(1, 3, [1, 2, 3])
        self.assertEqual(b[0], 1)
        self.assertListEqual(nc.to_list(b[0:2]), [1,2])
        with self.assertRaises(TypeError):
            b[0:1, 0:1]
        # 1-d matrix, 1 col
        c = nc.Matrix(3, 1, [1, 2, 3])
        self.assertEqual(c[-1], 3)
        self.assertListEqual(nc.to_list(c[-3:3]), [1,2,3])
        with self.assertRaises(TypeError):
            c[0:1, 0:1]
        # invalid slices for matrix61c
        d = nc.Matrix(4, 4)
        with self.assertRaises(ValueError):
            d[0:4:2]
        with self.assertRaises(ValueError):
            d[0:0]

class TestSet(TestCase):
    def test_set(self):
        # TODO: YOUR CODE HERE
        a = nc.Matrix(3, 3)
        a[0:1, 0:1] = 1.0
        self.assertEqual(a[0, 0], 1)
        a[:, 0] = [1, 1, 1]
        self.assertListEqual(nc.to_list(a), [[1,0,0],[1,0,0],[1,0,0]])
        a[0, :] = [2, 2, 2]
        self.assertListEqual(nc.to_list(a), [[2,2,2],[1,0,0],[1,0,0]])
        a[0:-1, 0:2] = [[1, 2], [3, 4]]
        self.assertListEqual(nc.to_list(a), [[1.0, 2.0, 2.0], [3.0, 4.0, 0.0], [1.0, 0.0, 0.0]])
        # test slices
        a = nc.Matrix(2, 2)
        a[0:1, 0:-1] = 1.0
        self.assertListEqual(nc.to_list(a), [[1.0, 0.0], [0.0, 0.0]])
        a[1] = [2, 2]
        self.assertListEqual(nc.to_list(a), [[1.0, 0.0], [2.0, 2.0]])
        b = a[1]
        b[-1] = 3
        self.assertListEqual(nc.to_list(a), [[1.0, 0.0], [2.0, 3.0]])
        # test nested slices
        a = nc.Matrix(4, 4)
        b = a[-4:-1, 0:3]
        c = b[1:3, 1:3]
        c[0] = [2, 2]
        self.assertListEqual(nc.to_list(a), [[0.0, 0.0, 0.0, 0.0], [0.0, 2.0, 2.0, 0.0], [0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0]])
        self.assertListEqual(nc.to_list(b), [[0.0, 0.0, 0.0], [0.0, 2.0, 2.0], [0.0, 0.0, 0.0]])
        self.assertListEqual(nc.to_list(c), [[2.0, 2.0], [0.0, 0.0]])

class TestShape(TestCase):
    def test_shape(self):
        np_mat, nc_mat = rand_np_nc_matrix(12, 23, seed=0)
        self.assertTrue(np_mat.shape == nc_mat.shape)
        self.assertEqual(nc_mat.shape, (12, 23))
