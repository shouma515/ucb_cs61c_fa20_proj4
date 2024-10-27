#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

// private function declarations
int set_identity(matrix*);
int cp_matrix(matrix *, matrix*);

/*
 * Generates a random double between `low` and `high`.
 */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/*
 * Generates a random matrix with `seed`.
 */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Allocate space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. Remember to set all fieds of the matrix struct.
 * `parent` should be set to NULL to indicate that this matrix is not a slice.
 * You should return -1 if either `rows` or `cols` or both have invalid values, or if any
 * call to allocate memory in this function fails. If you don't set python error messages here upon
 * failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    // use **mat here so that all the allocation happen in this function
    // otherwise if use *mat, need to allocate space for a matrix struct outside of
    // this function and then pass its reference in
    if (rows <= 0 || cols <= 0) {
        // invalid argument value
        return -1;
    }

    *mat = (matrix *)malloc(sizeof(matrix));
    if (*mat == NULL) {
        // fail to allocate space
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;

    (*mat)->data = (double **)malloc(rows * sizeof(double *));;
    if ((*mat)->data == NULL) {
        free(*mat);
        return -1;
    }

    for (int i = 0; i < rows; i++) {
        (*mat)->data[i] = (double *)malloc(cols * sizeof(double));
        if ((*mat)->data[i] == NULL) {
            // free all previous allocated rows
            for (int j = 0; j < i; j++) {
                free((*mat)->data[j]);
            }
            free((*mat)->data);
            free(*mat);
            return -1;
        }
    }

    (*mat)->is_1d = rows == 1 || cols == 1 ? 1 : 0;
    (*mat)->ref_cnt = 1; // a reference to itself
    (*mat)->parent = NULL;

    // initialize all values to be zero
    fill_matrix(*mat, 0);

    return 0;
}

/*
 * Allocate space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * This is equivalent to setting the new matrix to be
 * from[row_offset:row_offset + rows, col_offset:col_offset + cols]
 * If you don't set python error messages here upon failure, then remember to set it in numc.c.
 * Return 0 upon success and non-zero upon failure.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int row_offset, int col_offset,
                        int rows, int cols) {
    /* TODO: YOUR CODE HERE */
    // check if arguments are valid
    if (from == NULL || rows <= 0 || cols <= 0 || row_offset >= from->rows || col_offset >= from->cols
        || row_offset + rows > from->rows || col_offset + cols > from->cols) {
        return -1;
    }

    *mat = (matrix *)malloc(sizeof(matrix));
    if (*mat == NULL) {
        // fail to allocate space
        return -1;
    }
    (*mat)->rows = rows;
    (*mat)->cols = cols;

    (*mat)->data = (double **)malloc(rows * sizeof(double *));;
    if ((*mat)->data == NULL) {
        free(*mat);
        return -1;
    }

    for (int i = 0; i < rows; i++) {
        int idx = i + row_offset;
        (*mat)->data[i] = from->data[idx] + col_offset;
    }

    (*mat)->is_1d = rows == 1 || cols == 1 ? 1 : 0;
    (*mat)->ref_cnt = 1;
    (*mat)->parent = from;

    from->ref_cnt += 1;

    return 0;

}

/*
 * This function will be called automatically by Python when a numc matrix loses all of its
 * reference pointers.
 * You need to make sure that you only free `mat->data` if no other existing matrices are also
 * referring this data array.
 * See the spec for more information.
 */
void deallocate_matrix(matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (mat == NULL) {
        return;
    }

    mat->ref_cnt -= 1;
    if (mat->ref_cnt > 0) {
        return;
    }

    // need to deallocate current matrix
    if (mat->parent == NULL) {
        // this is a top level metrics, need to free all memory in data field
        for (int i = 0; i < mat->rows; i++) {
            free(mat->data[i]);
        }
        free(mat->data);
        free(mat);
    } else {
        // if the parent has more than 1 reference, this will decrease its ref count by one;
        // else, this will also deallocate parent, and maybe its parent recursively.
        // All the children of this matrix, indirectly reference to this matrix's parent
        // through this matrix, so we only try to deallocate its parent when ref to this matrix
        // decrease to zero, i.e. when all its children are deallocated and itself is deallocated.
        deallocate_matrix(mat->parent);
        // matrix with parent only responsible for the double* array to hold the start of rows of the matrix.
        free(mat->data);
        free(mat);
    }
}

/*
 * Return the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid.
 */
double get(matrix *mat, int row, int col) {
    /* TODO: YOUR CODE HERE */
    return mat->data[row][col];
}

/*
 * Set the value at the given row and column to val. You may assume `row` and
 * `col` are valid
 */
void set(matrix *mat, int row, int col, double val) {
    /* TODO: YOUR CODE HERE */
    mat->data[row][col] = val;
}

/*
 * Set all entries in mat to val
 */
void fill_matrix(matrix *mat, double val) {
    /* TODO: YOUR CODE HERE */
    if (val == 0) {
        for (int i = 0; i < mat->rows; i++) {
            memset(mat->data[i], 0, mat->cols * sizeof(double));
        }
        return;
    }

    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] = val;
        }
    }
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (result == NULL || mat1 == NULL || mat2 == NULL ||
        result->rows != mat1->rows || result->rows != mat2->rows
        || result->cols != mat1->cols || result->cols != mat2->cols) {
            return -1;
        }
    
    // method 1, straightforward
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
        }
    }

    // // method 2, unrolling
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols / 8 * 8; j+=8) {
    //         result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //         result->data[i][j+1] = mat1->data[i][j+1] + mat2->data[i][j+1];
    //         result->data[i][j+2] = mat1->data[i][j+2] + mat2->data[i][j+2];
    //         result->data[i][j+3] = mat1->data[i][j+3] + mat2->data[i][j+3];
    //         result->data[i][j+4] = mat1->data[i][j+4] + mat2->data[i][j+4];
    //         result->data[i][j+5] = mat1->data[i][j+5] + mat2->data[i][j+5];
    //         result->data[i][j+6] = mat1->data[i][j+6] + mat2->data[i][j+6];
    //         result->data[i][j+7] = mat1->data[i][j+7] + mat2->data[i][j+7];
    //     }
        
    //     // deal with the last elements when cols is not a multiply of the factor
    //     for (int j = result->cols / 8 * 8; j < result->cols; j++) {
    //         result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //     }
    // }

    // // method 3, SIMD
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols/4*4; j+=4) {
    //         __m256d op1 = _mm256_loadu_pd(mat1->data[i]+j);
    //         __m256d op2 = _mm256_loadu_pd(mat2->data[i]+j);
    //         __m256d r = _mm256_add_pd(op1, op2);
    //         _mm256_storeu_pd(result->data[i]+j, r);
    //     }

    //     for (int j = result->cols / 4 * 4; j < result->cols; j++) {
    //         result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //     }
    // }

    // // method 4, OpenMP, simple parallel for
    // #pragma omp parallel for
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols/4*4; j+=4) {
    //         __m256d op1 = _mm256_loadu_pd(mat1->data[i]+j);
    //         __m256d op2 = _mm256_loadu_pd(mat2->data[i]+j);
    //         __m256d r = _mm256_add_pd(op1, op2);
    //         _mm256_storeu_pd(result->data[i]+j, r);
    //     }

    //     for (int j = result->cols / 4 * 4; j < result->cols; j++) {
    //         result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //     }
    // }

    // // method 4, OpenMP, manual job division, the simple method above does
    // // not work when there are very few rows, but a lot cols
    // // TODO: error when row is 1 and 2
    // // TODO: this method also has problem when number of element of the matrix
    // // is less then total number of threads
    // #pragma omp parallel
    // {
    //     int omp_tid = omp_get_thread_num();
    //     int total_thread = omp_get_num_threads();
    //     printf("total_thread: %d, crnt_thread: %d\n", total_thread, omp_tid);
    //     fflush(stdout);
    //     int row_start, row_end, col_start, col_end;
    //     if (result->rows >= total_thread) {
    //         // TODO: figure out a better way to distribute
    //         // workload more evenly among threads
    //         // e.g., first using the naive method below, then
    //         // move workload from highest to lowest until the diff
    //         // between highest and lowest is one
    //         int block_size = result->rows / total_thread;
    //         row_start = block_size * omp_tid;
    //         row_end = row_start + block_size;
    //         if (omp_tid == total_thread - 1) {
    //             // last thread need to handle all remaining nums
    //             row_end = result->rows;
    //         }
    //         col_start = 0;
    //         col_end = result->cols;
    //     } else {
    //         // very few rows
    //         int thread_per_row = total_thread / result->rows;
    //         int additional_thread = total_thread % result->rows;
    //         if (omp_tid < additional_thread * (thread_per_row + 1)) {
    //             // first additional_thread rows get 1 additional thread
    //             row_start = omp_tid / (thread_per_row + 1);
    //             row_end = row_start + 1;
    //             int col_block_size = result->cols / (thread_per_row + 1);
    //             int col_idx = omp_tid % (thread_per_row + 1);
    //             col_start = col_idx * col_block_size;
    //             col_end = col_start + col_block_size;
    //             if (col_idx == thread_per_row) {
    //                 col_end = result->cols;
    //             }
    //         } else {
    //             // last few rows only get thread_per_row thread
    //             int offset_idx = omp_tid - additional_thread * (thread_per_row + 1);
    //             row_start = offset_idx / thread_per_row + additional_thread;
    //             row_end = row_start + 1;
    //             int col_block_size = result->cols / thread_per_row;
    //             int col_idx = offset_idx % thread_per_row;
    //             col_start = col_idx * col_block_size;
    //             col_end = col_start + col_block_size;
    //             if (col_idx == thread_per_row - 1) {
    //                 col_end = result->cols;
    //             }
    //         }
    //     }
    //     for (int i = row_start; i < row_end; i++) {
    //         for (int j = col_start; j < col_end/4*4; j+=4) {
    //             __m256d op1 = _mm256_loadu_pd(mat1->data[i]+j);
    //             __m256d op2 = _mm256_loadu_pd(mat2->data[i]+j);
    //             __m256d r = _mm256_add_pd(op1, op2);
    //             _mm256_storeu_pd(result->data[i]+j, r);
    //         }

    //         for (int j = col_end/4*4; j < col_end; j++) {
    //             result->data[i][j] = mat1->data[i][j] + mat2->data[i][j];
    //         }
    //     }
    // }

    return 0;
}

/*
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (result == NULL || mat1 == NULL || mat2 == NULL ||
        result->rows != mat1->rows || result->rows != mat2->rows
        || result->cols != mat1->cols || result->cols != mat2->cols) {
            return -1;
        }
    
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = mat1->data[i][j] - mat2->data[i][j];
        }
    }

    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 */
// mutiply two 1000*1000 matrices
// 1. raw method one takes 15-20s for 10 times, i.e., 1.5-2s
// 2. change ijk to ikj takes 3.0s - 3.5s for 10 times to start with, i.e., 300ms
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    /* TODO: YOUR CODE HERE */
    if (result == NULL || mat1 == NULL || mat2 == NULL ||
        result->rows != mat1->rows || result->cols != mat2->cols
        || mat1->cols != mat2->rows) {
            return -1;
        }
    
    // // TODO: test the speed of the following two implementation
    // // method 1: takes 15.7s for 10 times, i.e., 1.57s
    // for (int i = 0; i < result->rows; i++) {
    //     for (int j = 0; j < result->cols; j++) {
    //         double sum = 0;
    //         for (int k = 0; k < mat1->cols; k++) {
    //             sum += mat1->data[i][k] * mat2->data[k][j];
    //         }
    //         result->data[i][j] = sum;
    //     }
    // }

    // method 2: more cache friendly as we store data row-wise, 3.25s for 10 times
    // TODO: maybe try to initialize result->data[i][j] to zero in the inner loop
    // using a test on k == 0? so that no additional functional is needed.
    // fill_matrix(result, 0);
    // for (int i = 0; i < result->rows; i++) {
    //     for (int k = 0; k < mat1->cols; k++) {
    //         for (int j = 0; j < result->cols; j++) {
    //             result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
    //         }
    //     }
    // }

    // // method 3: block
    // // // method 4: SIMD, ~5s for 10 times, 500ms
    // fill_matrix(result, 0);
    // for (int i = 0; i < result->rows; i++) {
    //     for (int k = 0; k < mat1->cols; k++) {
    //         __m256d a = _mm256_set1_pd(mat1->data[i][k]);
    //         for (int j = 0; j < result->cols / 4 * 4; j+=4) {
    //             __m256d v_b = _mm256_loadu_pd(mat2->data[k] + j);
    //             __m256d v_c = _mm256_loadu_pd(mat2->data[i] + j);
    //             __m256d v_r = _mm256_fmadd_pd(a, v_b, v_c);
    //             _mm256_storeu_pd(result->data[i] + j, v_r);
    //         }

    //         for (int j = result->cols / 4 * 4; j < result->cols; j++) {
    //             result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
    //         }
    //     }
    // }
    // // method 5: SIMD + unroll, ~2.5s for 10 times, 250ms
    // fill_matrix(result, 0);
    // for (int i = 0; i < result->rows; i++) {
    //     for (int k = 0; k < mat1->cols / 4 * 4; k+=4) {
    //         __m256d a1 = _mm256_set1_pd(mat1->data[i][k]);
    //         __m256d a2 = _mm256_set1_pd(mat1->data[i][k+1]);
    //         __m256d a3 = _mm256_set1_pd(mat1->data[i][k+2]);
    //         __m256d a4 = _mm256_set1_pd(mat1->data[i][k+3]);
    //         for (int j = 0; j < result->cols / 4 * 4; j+=4) {
    //             __m256d v_b1 = _mm256_loadu_pd(mat2->data[k] + j);
    //             __m256d v_b2 = _mm256_loadu_pd(mat2->data[k+1] + j);
    //             __m256d v_b3 = _mm256_loadu_pd(mat2->data[k+2] + j);
    //             __m256d v_b4 = _mm256_loadu_pd(mat2->data[k+3] + j);
    //             __m256d v_r = _mm256_loadu_pd(mat2->data[i] + j);
    //             v_r = _mm256_fmadd_pd(a1, v_b1, v_r);
    //             v_r = _mm256_fmadd_pd(a2, v_b2, v_r);
    //             v_r = _mm256_fmadd_pd(a3, v_b3, v_r);
    //             v_r = _mm256_fmadd_pd(a4, v_b4, v_r);
    //             _mm256_storeu_pd(result->data[i] + j, v_r);
    //         }

    //         for (int j = result->cols / 4 * 4; j < result->cols; j++) {
    //             result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
    //             result->data[i][j] += mat1->data[i][k+1] * mat2->data[k+1][j];
    //             result->data[i][j] += mat1->data[i][k+2] * mat2->data[k+2][j];
    //             result->data[i][j] += mat1->data[i][k+3] * mat2->data[k+3][j];
    //         }
    //     }

    //     for (int k = mat1->cols / 4 * 4; k < mat1->cols; k++) {
    //         __m256d a = _mm256_set1_pd(mat1->data[i][k]);
    //         for (int j = 0; j < result->cols / 4 * 4; j+=4) {
    //             __m256d v_b = _mm256_loadu_pd(mat2->data[k] + j);
    //             __m256d v_c = _mm256_loadu_pd(mat2->data[i] + j);
    //             __m256d v_r = _mm256_fmadd_pd(a, v_b, v_c);
    //             _mm256_storeu_pd(result->data[i] + j, v_r);
    //         }

    //         for (int j = result->cols / 4 * 4; j < result->cols; j++) {
    //             result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
    //         }
    //     }
    // }

    // method 6: SIMD + unroll + OpenMP, ~0.4s for 10 times, 40ms
    fill_matrix(result, 0);
    #pragma omp parallel for
    for (int i = 0; i < result->rows; i++) {
        for (int k = 0; k < mat1->cols / 4 * 4; k+=4) {
            __m256d a1 = _mm256_set1_pd(mat1->data[i][k]);
            __m256d a2 = _mm256_set1_pd(mat1->data[i][k+1]);
            __m256d a3 = _mm256_set1_pd(mat1->data[i][k+2]);
            __m256d a4 = _mm256_set1_pd(mat1->data[i][k+3]);
            for (int j = 0; j < result->cols / 4 * 4; j+=4) {
                __m256d v_b1 = _mm256_loadu_pd(mat2->data[k] + j);
                __m256d v_b2 = _mm256_loadu_pd(mat2->data[k+1] + j);
                __m256d v_b3 = _mm256_loadu_pd(mat2->data[k+2] + j);
                __m256d v_b4 = _mm256_loadu_pd(mat2->data[k+3] + j);
                __m256d v_r = _mm256_loadu_pd(mat2->data[i] + j);
                v_r = _mm256_fmadd_pd(a1, v_b1, v_r);
                v_r = _mm256_fmadd_pd(a2, v_b2, v_r);
                v_r = _mm256_fmadd_pd(a3, v_b3, v_r);
                v_r = _mm256_fmadd_pd(a4, v_b4, v_r);
                _mm256_storeu_pd(result->data[i] + j, v_r);
            }

            for (int j = result->cols / 4 * 4; j < result->cols; j++) {
                result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
                result->data[i][j] += mat1->data[i][k+1] * mat2->data[k+1][j];
                result->data[i][j] += mat1->data[i][k+2] * mat2->data[k+2][j];
                result->data[i][j] += mat1->data[i][k+3] * mat2->data[k+3][j];
            }
        }

        for (int k = mat1->cols / 4 * 4; k < mat1->cols; k++) {
            __m256d a = _mm256_set1_pd(mat1->data[i][k]);
            for (int j = 0; j < result->cols / 4 * 4; j+=4) {
                __m256d v_b = _mm256_loadu_pd(mat2->data[k] + j);
                __m256d v_c = _mm256_loadu_pd(mat2->data[i] + j);
                __m256d v_r = _mm256_fmadd_pd(a, v_b, v_c);
                _mm256_storeu_pd(result->data[i] + j, v_r);
            }

            for (int j = result->cols / 4 * 4; j < result->cols; j++) {
                result->data[i][j] += mat1->data[i][k] * mat2->data[k][j];
            }
        }
    }

    // TODO: optimization
    // 1. unroll
    // 2. block-wise multiply
    // 3. use memset to fill_matrix with zero

    return 0;
}


/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 */
// pow a 1000 * 1000 matrix to 3, start with 12.4s for 10 times, i.e. 1.24s
// 1. reorder base calculation to reduce 1 unnecessary multiplication, reduce time to 10.0s
// 2. with multiplication method 6, reduce pow to 1.36s for 10 times, 136ms
int pow_matrix(matrix *result, matrix *mat, int pow) {
    /* TODO: YOUR CODE HERE */
    // use fast pow: https://zhuanlan.zhihu.com/p/42639682
    if (result == NULL || mat == NULL || pow < 0 ||
        result->rows != mat->rows || result->cols != mat->cols
        || mat->rows != mat->cols) {
            return -1;
        }
    // code to specially treat pow = 0 or 1, so that those two cases
    // can be faster.
    if (pow == 1) {
        cp_matrix(result, mat);
        return 0;
    }
    set_identity(result);
    if (pow == 0) {
        return 0;
    }
    matrix *tmp = NULL; // use as space to hold multiplication result
    matrix *base = NULL; // base of pow calulation
    if (allocate_matrix(&tmp, mat->rows, mat->rows) != 0) {
        return -1;
    }
    if (allocate_matrix(&base, mat->rows, mat->rows) != 0) {
        return -1;
    }
    // ans of pow (here we use result mat so that we don't need to allocate another matrix)
    matrix *ans = result;
    int first_base = 1;
    while (pow)
    {
        if (first_base) {
            cp_matrix(base, mat);
            first_base = 0;
        } else {
            // calculate current base, base = previous_base ** 2
            // we do this complex if else, so that we can remove 1 unnecessary matrix multiplication,
            // comparing with do calculate base at the end of each loop(which do not require
            // test if this is the first base, i.e., mat**1)
            mul_matrix(tmp, base, base);
            matrix *t = base;
            base = tmp;
            tmp = t;
        }
        if (pow & 1) {
            // ans = ans * base
            mul_matrix(tmp, ans, base);
            matrix *t = ans;
            ans = tmp;
            tmp = t;
        }
        pow >>= 1;
    }

    if (ans != result) {
        cp_matrix(result, ans);
        deallocate_matrix(ans);
    }

    if (base != result) {
        deallocate_matrix(base);
    }

    if (tmp != result) {
        deallocate_matrix(tmp);
    }
    
    return 0;
}

/*
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int neg_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (result == NULL || mat == NULL ||
        result->rows != mat->rows || result->cols != mat->cols) {
            return -1;
        }
    
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = -mat->data[i][j];
        }
    }

    return 0;
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success and a nonzero value upon failure.
 */
int abs_matrix(matrix *result, matrix *mat) {
    /* TODO: YOUR CODE HERE */
    if (result == NULL || mat == NULL ||
        result->rows != mat->rows || result->cols != mat->cols) {
            return -1;
        }
    
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            result->data[i][j] = fabs(mat->data[i][j]);
        }
    }

    return 0;
}

int set_identity(matrix* mat) {
    if (mat == NULL || mat->rows != mat->cols) {
        return -1;
    }
    for (int i = 0; i < mat->rows; i++) {
        for (int j = 0; j < mat->cols; j++) {
            mat->data[i][j] = i == j ? 1 : 0;
        }
    }

    return 0;
}

int cp_matrix(matrix *dest, matrix* src) {
    if (dest == NULL || src == NULL ||
        dest->rows != src->rows || dest->cols != src->cols) {
            return -1;
        }
    
    for (int i = 0; i < dest->rows; i++) {
        for (int j = 0; j < dest->cols; j++) {
            dest->data[i][j] = src->data[i][j];
        }
    }

    return 0;
}