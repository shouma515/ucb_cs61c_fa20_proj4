#include "numc.h"
#include <structmember.h>

PyTypeObject Matrix61cType;

PyObject *get_shape(int rows, int cols);
int check_set_scalar(matrix* mat, int row, int col, PyObject *v);
int check_1D(int length, PyObject *v);
void set_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v);
int check_set_row_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v);
int check_set_col_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v);
int check_set_2D(matrix* mat, int row_offset, int col_offset, int rows, int cols, PyObject* v);

/* Helper functions for initalization of matrices and vectors */

/**
 * Create a Matrix61c instance with the given matrix, assume input matrix is valid.
 */
Matrix61c *create_matrix61c(matrix * mat) {
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = mat;
    ret->shape = get_shape(mat->rows, mat->cols);
    return ret;
} 


/*
 * Return a tuple given rows and cols
 */
PyObject *get_shape(int rows, int cols) {
  if (rows == 1 || cols == 1) {
    return PyTuple_Pack(1, PyLong_FromLong(rows * cols));
  } else {
    return PyTuple_Pack(2, PyLong_FromLong(rows), PyLong_FromLong(cols));
  }
}
/*
 * Matrix(rows, cols, low, high). Fill a matrix random double values
 */
int init_rand(PyObject *self, int rows, int cols, unsigned int seed, double low,
              double high) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    rand_matrix(new_mat, seed, low, high);
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(rows, cols, val). Fill a matrix of dimension rows * cols with val
 */
int init_fill(PyObject *self, int rows, int cols, double val) {
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed)
        return alloc_failed;
    else {
        fill_matrix(new_mat, val);
        ((Matrix61c *)self)->mat = new_mat;
        ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    }
    return 0;
}

/*
 * Matrix(rows, cols, 1d_list). Fill a matrix with dimension rows * cols with 1d_list values
 */
int init_1d(PyObject *self, int rows, int cols, PyObject *lst) {
    if (rows * cols != PyList_Size(lst)) {
        PyErr_SetString(PyExc_ValueError, "Incorrect number of elements in list");
        return -1;
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    int count = 0;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j, PyFloat_AsDouble(PyList_GetItem(lst, count)));
            count++;
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * Matrix(2d_list). Fill a matrix with dimension len(2d_list) * len(2d_list[0])
 */
int init_2d(PyObject *self, PyObject *lst) {
    int rows = PyList_Size(lst);
    if (rows == 0) {
        PyErr_SetString(PyExc_ValueError,
                        "Cannot initialize numc.Matrix with an empty list");
        return -1;
    }
    int cols;
    if (!PyList_Check(PyList_GetItem(lst, 0))) {
        PyErr_SetString(PyExc_ValueError, "List values not valid");
        return -1;
    } else {
        cols = PyList_Size(PyList_GetItem(lst, 0));
    }
    for (int i = 0; i < rows; i++) {
        if (!PyList_Check(PyList_GetItem(lst, i)) ||
                PyList_Size(PyList_GetItem(lst, i)) != cols) {
            PyErr_SetString(PyExc_ValueError, "List values not valid");
            return -1;
        }
    }
    matrix *new_mat;
    int alloc_failed = allocate_matrix(&new_mat, rows, cols);
    if (alloc_failed) return alloc_failed;
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            set(new_mat, i, j,
                PyFloat_AsDouble(PyList_GetItem(PyList_GetItem(lst, i), j)));
        }
    }
    ((Matrix61c *)self)->mat = new_mat;
    ((Matrix61c *)self)->shape = get_shape(new_mat->rows, new_mat->cols);
    return 0;
}

/*
 * This deallocation function is called when reference count is 0
 */
void Matrix61c_dealloc(Matrix61c *self) {
    deallocate_matrix(self->mat);
    Py_TYPE(self)->tp_free(self);
}

/* For immutable types all initializations should take place in tp_new */
PyObject *Matrix61c_new(PyTypeObject *type, PyObject *args,
                        PyObject *kwds) {
    /* size of allocated memory is tp_basicsize + nitems*tp_itemsize*/
    Matrix61c *self = (Matrix61c *)type->tp_alloc(type, 0);
    return (PyObject *)self;
}

/*
 * This matrix61c type is mutable, so needs init function. Return 0 on success otherwise -1
 */
int Matrix61c_init(PyObject *self, PyObject *args, PyObject *kwds) {
    /* Generate random matrices */
    if (kwds != NULL) {
        PyObject *rand = PyDict_GetItemString(kwds, "rand");
        if (!rand) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (!PyBool_Check(rand)) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
        if (rand != Py_True) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        PyObject *low = PyDict_GetItemString(kwds, "low");
        PyObject *high = PyDict_GetItemString(kwds, "high");
        PyObject *seed = PyDict_GetItemString(kwds, "seed");
        double double_low = 0;
        double double_high = 1;
        unsigned int unsigned_seed = 0;

        if (low) {
            if (PyFloat_Check(low)) {
                double_low = PyFloat_AsDouble(low);
            } else if (PyLong_Check(low)) {
                double_low = PyLong_AsLong(low);
            }
        }

        if (high) {
            if (PyFloat_Check(high)) {
                double_high = PyFloat_AsDouble(high);
            } else if (PyLong_Check(high)) {
                double_high = PyLong_AsLong(high);
            }
        }

        if (double_low >= double_high) {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }

        // Set seed if argument exists
        if (seed) {
            if (PyLong_Check(seed)) {
                unsigned_seed = PyLong_AsUnsignedLong(seed);
            }
        }

        PyObject *rows = NULL;
        PyObject *cols = NULL;
        if (PyArg_UnpackTuple(args, "args", 2, 2, &rows, &cols)) {
            if (rows && cols && PyLong_Check(rows) && PyLong_Check(cols)) {
                return init_rand(self, PyLong_AsLong(rows), PyLong_AsLong(cols), unsigned_seed, double_low,
                                 double_high);
            }
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    }
    PyObject *arg1 = NULL;
    PyObject *arg2 = NULL;
    PyObject *arg3 = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 3, &arg1, &arg2, &arg3)) {
        /* arguments are (rows, cols, val) */
        if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && (PyLong_Check(arg3)
                || PyFloat_Check(arg3))) {
            if (PyLong_Check(arg3)) {
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyLong_AsLong(arg3));
            } else
                return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), PyFloat_AsDouble(arg3));
        } else if (arg1 && arg2 && arg3 && PyLong_Check(arg1) && PyLong_Check(arg2) && PyList_Check(arg3)) {
            /* Matrix(rows, cols, 1D list) */
            return init_1d(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), arg3);
        } else if (arg1 && PyList_Check(arg1) && arg2 == NULL && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_2d(self, arg1);
        } else if (arg1 && arg2 && PyLong_Check(arg1) && PyLong_Check(arg2) && arg3 == NULL) {
            /* Matrix(rows, cols, 1D list) */
            return init_fill(self, PyLong_AsLong(arg1), PyLong_AsLong(arg2), 0);
        } else {
            PyErr_SetString(PyExc_TypeError, "Invalid arguments");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return -1;
    }
}

/*
 * List of lists representations for matrices
 */
PyObject *Matrix61c_to_list(Matrix61c *self) {
    int rows = self->mat->rows;
    int cols = self->mat->cols;
    PyObject *py_lst = NULL;
    if (self->mat->is_1d) {  // If 1D matrix, print as a single list
        py_lst = PyList_New(rows * cols);
        int count = 0;
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(py_lst, count, PyFloat_FromDouble(get(self->mat, i, j)));
                count++;
            }
        }
    } else {  // if 2D, print as nested list
        py_lst = PyList_New(rows);
        for (int i = 0; i < rows; i++) {
            PyList_SetItem(py_lst, i, PyList_New(cols));
            PyObject *curr_row = PyList_GetItem(py_lst, i);
            for (int j = 0; j < cols; j++) {
                PyList_SetItem(curr_row, j, PyFloat_FromDouble(get(self->mat, i, j)));
            }
        }
    }
    return py_lst;
}

PyObject *Matrix61c_class_to_list(Matrix61c *self, PyObject *args) {
    PyObject *mat = NULL;
    if (PyArg_UnpackTuple(args, "args", 1, 1, &mat)) {
        if (!PyObject_TypeCheck(mat, &Matrix61cType)) {
            PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
            return NULL;
        }
        Matrix61c* mat61c = (Matrix61c*)mat;
        return Matrix61c_to_list(mat61c);
    } else {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments");
        return NULL;
    }
}

/*
 * Add class methods
 */
PyMethodDef Matrix61c_class_methods[] = {
    {"to_list", (PyCFunction)Matrix61c_class_to_list, METH_VARARGS, "Returns a list representation of numc.Matrix"},
    {NULL, NULL, 0, NULL}
};

/*
 * Matrix61c string representation. For printing purposes.
 */
PyObject *Matrix61c_repr(PyObject *self) {
    PyObject *py_lst = Matrix61c_to_list((Matrix61c *)self);
    return PyObject_Repr(py_lst);
}

/* NUMBER METHODS */

/*
 * Add the second numc.Matrix (Matrix61c) object to the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_add(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c* mat61c = (Matrix61c*)args;
    if (self->mat->rows != mat61c->mat->rows || self->mat->cols != mat61c->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Matrice's dimensions must match when add.");
        return NULL;
    }
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, self->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix when add.");
        return NULL;
    }
    if (add_matrix(result, self->mat, mat61c->mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing add.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * Substract the second numc.Matrix (Matrix61c) object from the first one. The first operand is
 * self, and the second operand can be obtained by casting `args`.
 */
PyObject *Matrix61c_sub(Matrix61c* self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c* mat61c = (Matrix61c*)args;
    if (self->mat->rows != mat61c->mat->rows || self->mat->cols != mat61c->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Matrice's dimensions must match when sub.");
        return NULL;
    }
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, self->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix when sub.");
        return NULL;
    }
    if (sub_matrix(result, self->mat, mat61c->mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing sub.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * NOT element-wise multiplication. The first operand is self, and the second operand
 * can be obtained by casting `args`.
 */
PyObject *Matrix61c_multiply(Matrix61c* self, PyObject *args) {
    /* TODO: YOUR CODE HERE */
    if (!PyObject_TypeCheck(args, &Matrix61cType)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type numc.Matrix!");
        return NULL;
    }
    Matrix61c* mat61c = (Matrix61c*)args;
    if (self->mat->cols != mat61c->mat->rows) {
        PyErr_SetString(PyExc_ValueError, "Matrice's dimensions must match when multiply.");
        return NULL;
    }
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, mat61c->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix when multiply.");
        return NULL;
    }
    if (mul_matrix(result, self->mat, mat61c->mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing multiply.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * Negates the given numc.Matrix.
 */
PyObject *Matrix61c_neg(Matrix61c* self) {
    /* TODO: YOUR CODE HERE */
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, self->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix.");
        return NULL;
    }
    if (neg_matrix(result, self->mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing neg.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * Take the element-wise absolute value of this numc.Matrix.
 */
PyObject *Matrix61c_abs(Matrix61c *self) {
    /* TODO: YOUR CODE HERE */
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, self->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix.");
        return NULL;
    }
    if (abs_matrix(result, self->mat)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing abs.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * Raise numc.Matrix (Matrix61c) to the `pow`th power. You can ignore the argument `optional`.
 */
PyObject *Matrix61c_pow(Matrix61c *self, PyObject *pow, PyObject *optional) {
    /* TODO: YOUR CODE HERE */
    if (self->mat->rows != self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Power function only works for square matrices.");
        return NULL;
    }
    if (!PyLong_Check(pow)) {
        PyErr_SetString(PyExc_TypeError, "Argument must of type long.");
        return NULL;
    }
    long val_pow = PyLong_AsLong(pow);
    if (val_pow < 0) {
        PyErr_SetString(PyExc_ValueError, "Pow must be non-negative.");
        return NULL;
    }
    matrix* result = NULL;
    if (allocate_matrix(&result, self->mat->rows, self->mat->cols)) {
        PyErr_SetString(PyExc_RuntimeError, "Fail to allocate result matrix when pow.");
        return NULL;
    }
    if (pow_matrix(result, self->mat, val_pow)) {
        PyErr_SetString(PyExc_RuntimeError, "Error when executing pow.");
        return NULL;
    }
    Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
    ret->mat = result;
    ret->shape = get_shape(result->rows, result->cols);
    return ret;
}

/*
 * Create a PyNumberMethods struct for overloading operators with all the number methods you have
 * define. You might find this link helpful: https://docs.python.org/3.6/c-api/typeobj.html
 */
PyNumberMethods Matrix61c_as_number = {
    /* TODO: YOUR CODE HERE */
    .nb_add = Matrix61c_add,
    .nb_subtract = Matrix61c_sub,
    .nb_multiply = Matrix61c_multiply,
    .nb_negative = Matrix61c_neg,
    .nb_absolute = Matrix61c_abs,
    .nb_power = Matrix61c_pow
};


/* INSTANCE METHODS */

/*
 * Given a numc.Matrix self, parse `args` to (int) row, (int) col, and (double/int) val.
 * Return None in Python (this is different from returning null).
 */
PyObject *Matrix61c_set_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    int row, col;
    PyObject *val_obj = NULL;
    if (!PyArg_ParseTuple(args, "iiO", &row, &col, &val_obj)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments when set value");
        return NULL;
    }
    if (row < 0 || row >= self->mat->rows || col < 0 || col >= self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Index out of range");
        return NULL;
    }
    double val;
    if (PyLong_Check(val_obj)) {
        // It's an int
        val = PyLong_AsLong(val_obj); // cast to double here since the set method accept a double argument.
    } else if (PyFloat_Check(val_obj)) {
        // It's a double
        val = PyFloat_AsDouble(val_obj);
    } else {
        PyErr_SetString(PyExc_TypeError, "Value must be an int or a float");
        return NULL;
    }

    set(self->mat, row, col, val);

    // this is the macro to return Python None.
    Py_RETURN_NONE;
}

/*
 * Given a numc.Matrix `self`, parse `args` to (int) row and (int) col.
 * Return the value at the `row`th row and `col`th column, which is a Python
 * float/int.
 */
PyObject *Matrix61c_get_value(Matrix61c *self, PyObject* args) {
    /* TODO: YOUR CODE HERE */
    int row, col;
    if (!PyArg_ParseTuple(args, "ii", &row, &col)) {
        PyErr_SetString(PyExc_TypeError, "Invalid arguments when set value");
        return NULL;
    }
    if (row < 0 || row >= self->mat->rows || col < 0 || col >= self->mat->cols) {
        PyErr_SetString(PyExc_ValueError, "Index out of range");
        return NULL;
    }

    double val = get(self->mat, row, col);

    return PyFloat_FromDouble(val);
}

/*
 * Create an array of PyMethodDef structs to hold the instance methods.
 * Name the python function corresponding to Matrix61c_get_value as "get" and Matrix61c_set_value
 * as "set"
 * You might find this link helpful: https://docs.python.org/3.6/c-api/structures.html
 */
PyMethodDef Matrix61c_methods[] = {
    /* TODO: YOUR CODE HERE */
    // {"add", (PyCFunction)Matrix61c_add, METH_VARARGS, "Add a matrix with matching dimension to the current matrix, and produce a new numc.Matrix object."},
    // {"sub", (PyCFunction)Matrix61c_sub, METH_VARARGS, "Sub a matrix with matching dimension to the current matrix, and produce a new numc.Matrix object."},
    // {"multiply", (PyCFunction)Matrix61c_multiply, METH_VARARGS, "Multiply a matrix with matching dimension to the current matrix, and produce a new numc.Matrix object."},
    // {"neg", (PyCFunction)Matrix61c_neg, METH_NOARGS, "Negate the current matrix, and produce a new numc.Matrix object."},
    // {"abs", (PyCFunction)Matrix61c_abs, METH_NOARGS, "Abs the current matrix, and produce a new numc.Matrix object."},
    // {"pow", (PyCFunction)Matrix61c_pow, METH_KEYWORDS, "Raise the current matrix to pow th, and produce a new numc.Matrix object."},
    {"set", (PyCFunction)Matrix61c_set_value, METH_VARARGS, "Set the value at index (i, j) of the current matrix to val."},
    {"get", (PyCFunction)Matrix61c_get_value, METH_VARARGS, "Get the value of the current matrix at index (i, j)."},
    {NULL, NULL, 0, NULL}
};

/* INDEXING */

/*
 * Given a numc.Matrix `self`, index into it with `key`. Return the indexed result.
 */
PyObject *Matrix61c_subscript(Matrix61c* self, PyObject* key) {
    /* TODO: YOUR CODE HERE */
    // TODO:  is it ok to use PyIndex_Check here directly or it should be parsed from argument tuple first using PyArgs_ParseTuple?
    if (PyIndex_Check(key)) {
        // if argument is a int index
        Py_ssize_t i = PyNumber_AsSsize_t(key, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred()) {
            return NULL;
        }
        if (self->mat->is_1d) {
            if (i < 0) {
                i += self->mat->rows * self->mat->cols;
            }
            if (i < 0 || i >= self->mat->rows * self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "sequence index out of range");
                return NULL;
            }
            if (self->mat->rows == 1) {
                return PyFloat_FromDouble(get(self->mat, 0, i));
            } else {
                return PyFloat_FromDouble(get(self->mat, i, 0));
            }
        } else {
            if (i < 0) {
                i += self->mat->rows;
            }
            if (i < 0 || i >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "sequence index out of range");
                return NULL;
            }
            matrix* result = NULL;
            if (allocate_matrix_ref(&result, self->mat, i, 0, 1, self->mat->cols)) {
                PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                return NULL;
            }
            Matrix61c* ret = (Matrix61c*)Matrix61c_new(&Matrix61cType, NULL, NULL);
            ret->mat = result;
            ret->shape = get_shape(result->rows, result->cols);
            return ret;
        }
    }
    else if (PySlice_Check(key)) {
        // if key is a single slice
        Py_ssize_t start, stop, step, slicelength;
        if (PySlice_GetIndicesEx(key, self->mat->rows,
                                 &start, &stop, &step, &slicelength) < 0) {
            return NULL;
        }
        if (slicelength == 0 || step != 1) {
            PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
            return NULL;
        }
        matrix* result = NULL;
        if (self->mat->is_1d) {
            if (self->mat->rows == 1) {
                if (slicelength == 1) {
                    return PyFloat_FromDouble(get(self->mat, 0, start));
                }
                if (allocate_matrix_ref(&result, self->mat, 0, start, 1, slicelength)) {
                    PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                    return NULL;
                }
            } else {
                if (slicelength == 1) {
                    return PyFloat_FromDouble(get(self->mat, start, 0));
                }
                if (allocate_matrix_ref(&result, self->mat, start, 0, slicelength, 1)) {
                    PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                    return NULL;
                }
            }
        } else {
            if (allocate_matrix_ref(&result, self->mat, start, 0, slicelength, self->mat->cols)) {
                PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                return NULL;
            }
        }
        return create_matrix61c(result);
    } else if (PyTuple_Check(key) && PyTuple_GET_SIZE(key) == 2) {
        if (self->mat->is_1d) {
            PyErr_SetString(PyExc_TypeError, "for 1D matrix, key must be an int, a slice.");
            return NULL;
        }
        int row, col;
        PyObject* row_slice;
        PyObject* col_slice;
        // TODO: the following code has weird error when slice involving "-1", like [-1, :], [0, 0:-1], investigate it
        if (PyArg_ParseTuple(key, "ii", &row, &col)) {
            if (row < 0) {
                row += self->mat->rows;
            }
            if (col < 0) {
                col += self->mat->cols;
            }
            if (row < 0 || row >= self->mat->rows || col < 0 || col >= self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return NULL;
            }
            return PyFloat_FromDouble(get(self->mat, row, col));
        }
        else if (PyArg_ParseTuple(key, "iO!", &row, &PySlice_Type, &col_slice)) {
            if (row < 0) {
                row += self->mat->rows;
            }
            if (row < 0 || row >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return NULL;
            }
            Py_ssize_t start, stop, step, slicelength;
            if (PySlice_GetIndicesEx(col_slice, self->mat->cols,
                                    &start, &stop, &step, &slicelength) < 0) {
                return NULL;
            }
            if (slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return NULL;
            }
            if (slicelength == 1) {
                return PyFloat_FromDouble(get(self->mat, row, start));
            }
            matrix* result = NULL;
            if (allocate_matrix_ref(&result, self->mat, row, start, 1, slicelength)) {
                PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                return NULL;
            }
            return create_matrix61c(result);
        }
        else if (PyArg_ParseTuple(key, "O!i", &PySlice_Type, &row_slice, &col)) {
            if (col < 0) {
                col += self->mat->cols;
            }
            if (col < 0 || col >= self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return NULL;
            }
            Py_ssize_t start, stop, step, slicelength;
            if (PySlice_GetIndicesEx(row_slice, self->mat->rows,
                                    &start, &stop, &step, &slicelength) < 0) {
                return NULL;
            }
            if (slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return NULL;
            }
            if (slicelength == 1) {
                return PyFloat_FromDouble(get(self->mat, start, col));
            }
            matrix* result = NULL;
            if (allocate_matrix_ref(&result, self->mat, start, col, slicelength, 1)) {
                PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                return NULL;
            }
            return create_matrix61c(result);
        }
        else if (PyArg_ParseTuple(key, "O!O!", &PySlice_Type, &row_slice, &PySlice_Type, &col_slice)) {
            Py_ssize_t row_start, stop, step, row_slicelength;
            if (PySlice_GetIndicesEx(row_slice, self->mat->rows,
                                    &row_start, &stop, &step, &row_slicelength) < 0) {
                return NULL;
            }
            if (row_slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return NULL;
            }
            Py_ssize_t col_start, col_slicelength;
            if (PySlice_GetIndicesEx(col_slice, self->mat->cols,
                                    &col_start, &stop, &step, &col_slicelength) < 0) {
                return NULL;
            }
            if (col_slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return NULL;
            }

            if (row_slicelength == 1 && col_slicelength == 1) {
                return PyFloat_FromDouble(get(self->mat, row_start, col_start));
            }
            matrix* result = NULL;
            if (allocate_matrix_ref(&result, self->mat, row_start, col_start, row_slicelength, col_slicelength)) {
                PyErr_SetString(PyExc_RuntimeError, "Error when allocate matrix ref.");
                return NULL;
            }
            return create_matrix61c(result);
        } else {
            PyErr_SetString(PyExc_TypeError, "key of tuple type must be a slice/int tuple of size 2(for 2D matrix only)");
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "key must be an int, a slice, or a slice/int tuple of size 2(for 2D matrix only)");
        return NULL;
    }
}

/*
 * Given a numc.Matrix `self`, index into it with `key`, and set the indexed result to `v`.
 */
int Matrix61c_set_subscript(Matrix61c* self, PyObject *key, PyObject *v) {
    /* TODO: YOUR CODE HERE */
    if (PyIndex_Check(key)) {
        // if argument is a int index
        Py_ssize_t i = PyNumber_AsSsize_t(key, PyExc_IndexError);
        if (i == -1 && PyErr_Occurred()) {
            return -1;
        }
        if (self->mat->is_1d) {
            if (i < 0) {
                i += self->mat->rows * self->mat->cols;
            }
            if (i < 0 || i >= self->mat->rows * self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "sequence index out of range");
                return -1;
            }

            if (self->mat->rows == 1) {
                return check_set_scalar(self->mat, 0, i, v);
            }
            return check_set_scalar(self->mat, i, 0, v);
        } else {
            if (i < 0) {
                i += self->mat->rows;
            }
            if (i < 0 || i >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "sequence index out of range");
                return -1;
            }

            if (check_1D(self->mat->cols, v) != 0) {
                return -1;
            }

            set_1D(self->mat, i, 0, self->mat->cols, v);
            return 0;
        }
    }
    else if (PySlice_Check(key)) {
        // if key is a single slice
        Py_ssize_t start, stop, step, slicelength;
        if (PySlice_GetIndicesEx(key, self->mat->rows,
                                 &start, &stop, &step, &slicelength) < 0) {
            return -1;
        }
        if (slicelength == 0 || step != 1) {
            PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
            return -1;
        }
        if (self->mat->is_1d) {
            if (slicelength == 1) {
                return self->mat->rows == 1 ?
                    check_set_scalar(self->mat, 0, start, v) :
                    check_set_scalar(self->mat, start, 0, v);
            }
            return self->mat->rows == 1 ?
                check_set_row_1D(self->mat, 0, start, slicelength, v) :
                check_set_col_1D(self->mat, start, 0, slicelength, v);
        }
        return slicelength == 1 ?
            check_set_row_1D(self->mat, start, 0, self->mat->cols, v) :
            check_set_2D(self->mat, start, 0, slicelength, self->mat->cols, v);
    } else if (PyTuple_Check(key) && PyTuple_GET_SIZE(key) == 2) {
        if (self->mat->is_1d) {
            PyErr_SetString(PyExc_TypeError, "for 1D matrix, key must be an int, a slice.");
            return -1;
        }
        int row, col;
        PyObject* row_slice;
        PyObject* col_slice;
        if (PyArg_ParseTuple(key, "ii", &row, &col)) {
            if (row < 0) {
                row += self->mat->rows;
            }
            if (col < 0) {
                col += self->mat->rows;
            }
            if (row < 0 || row >= self->mat->rows || col < 0 || col >= self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return -1;
            }
            return check_set_scalar(self->mat, row, col, v);
        } else if (PyArg_ParseTuple(key, "iO!", &row, &PySlice_Type, &col_slice)) {
            if (row < 0) {
                row += self->mat->rows;
            }
            if (row < 0 || row >= self->mat->rows) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return -1;
            }
            Py_ssize_t start, stop, step, slicelength;
            if (PySlice_GetIndicesEx(col_slice, self->mat->cols,
                                    &start, &stop, &step, &slicelength) < 0) {
                return -1;
            }
            if (slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return -1;
            }
            if (slicelength == 1) {
                return check_set_scalar(self->mat, row, start, v);
            }
            return check_set_row_1D(self->mat, row, start, slicelength, v);
        }
        else if (PyArg_ParseTuple(key, "O!i", &PySlice_Type, &row_slice, &col)) {
            if (col < 0) {
                col += self->mat->cols;
            }
            if (col < 0 || col >= self->mat->cols) {
                PyErr_SetString(PyExc_IndexError, "Index out of range.");
                return -1;
            }
            Py_ssize_t start, stop, step, slicelength;
            if (PySlice_GetIndicesEx(row_slice, self->mat->rows,
                                    &start, &stop, &step, &slicelength) < 0) {
                return -1;
            }
            if (slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return -1;
            }
            if (slicelength == 1) {
                return check_set_scalar(self->mat, start, col, v);
            }
            return check_set_col_1D(self->mat, start, col, slicelength, v);
        }
        else if (PyArg_ParseTuple(key, "O!O!", &PySlice_Type, &row_slice, &PySlice_Type, &col_slice)) {
            Py_ssize_t row_start, stop, step, row_slicelength;
            if (PySlice_GetIndicesEx(row_slice, self->mat->rows,
                                    &row_start, &stop, &step, &row_slicelength) < 0) {
                return -1;
            }
            if (row_slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return -1;
            }
            Py_ssize_t col_start, col_slicelength;
            if (PySlice_GetIndicesEx(col_slice, self->mat->cols,
                                    &col_start, &stop, &step, &col_slicelength) < 0) {
                return -1;
            }
            if (col_slicelength == 0 || step != 1) {
                PyErr_SetString(PyExc_ValueError, "slice not valid for Matrix61c case.");
                return -1;
            }

            if (row_slicelength == 1 && col_slicelength == 1) {
                return check_set_scalar(self->mat, row_start, col_start, v);
            }
            if (row_slicelength == 1) {
                return check_set_row_1D(self->mat, row_start, col_start, col_slicelength, v);
            }
            if (col_slicelength == 1) {
                return check_set_col_1D(self->mat, row_start, col_start, row_slicelength, v);
            }
            return check_set_2D(self->mat, row_start, col_start, row_slicelength, col_slicelength, v);
        } else {
            PyErr_SetString(PyExc_TypeError, "key of tuple type must be a slice/int tuple of size 2(for 2D matrix only)");
            return -1;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "key must be an int, a slice, or a slice/int tuple of size 2(for 2D matrix only)");
        return -1;
    }
}

// assume offsets are valid
int check_set_scalar(matrix* mat, int row, int col, PyObject *v) {
    double value;
    if (PyFloat_Check(v)) {
        value = PyFloat_AsDouble(v);
    } else if (PyLong_Check(v)) {
        value = PyLong_AsLong(v);
    } else {
        PyErr_SetString(PyExc_ValueError, "Argument must be float or int");
        return -1;
    }

    set(mat, row, col, value);

    return 0;
}

// assume length is valid
int check_1D(int length, PyObject *v) {
    if (!PyList_Check(v)) {
        PyErr_SetString(PyExc_TypeError, "Expected a list");
        return -1;
    }

    // Get the length of the list
    Py_ssize_t len = PyList_Size(v);
    if (len != length) {
        PyErr_SetString(PyExc_ValueError, "Argument list length does not index length");
        return -1;
    }

    // Iterate through each item in the list
    for (Py_ssize_t i = 0; i < len; i++) {
        PyObject *item = PyList_GetItem(v, i);
        
        // Check if the item is an int or float
        if (!PyLong_Check(item) && !PyFloat_Check(item)) {
            PyErr_SetString(PyExc_ValueError, "All items in the list must be int or float");
            return -1;
        }
    }

    return 0;
}

// assume input arguments are all valid (pass check_1D)
void set_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v) {
    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *item = PyList_GetItem(v, i);
        
        double value;
        if (PyFloat_Check(item)) {
            value = PyFloat_AsDouble(item);
        } else {
            // must be a double or int as checked in the above iteration
            value = PyLong_AsLong(item);
        }
        set(mat, row_offset, col_offset + i, value);
    }
}

// assume offsets and lengths are valid
int check_set_row_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v) {
    if (check_1D(length, v) != 0) {
        return -1;
    }

    set_1D(mat, row_offset, col_offset, length, v);
    return 0;
}

// assume offsets and lengths are valid
int check_set_col_1D(matrix* mat, int row_offset, int col_offset, int length, PyObject *v) {
    if (check_1D(length, v) != 0) {
        return -1;
    }

    for (Py_ssize_t i = 0; i < length; i++) {
        PyObject *item = PyList_GetItem(v, i);
        
        double value;
        if (PyFloat_Check(item)) {
            value = PyFloat_AsDouble(item);
        } else {
            // must be a double or int as checked in the above iteration
            value = PyLong_AsLong(item);
        }
        set(mat, row_offset + i, col_offset, value);
    }
    return 0;
}

// assume offsets and lengths are valid
int check_set_2D(matrix* mat, int row_offset, int col_offset, int rows, int cols, PyObject* v) {
    if (!PyList_Check(v)) {
        PyErr_SetString(PyExc_TypeError, "Expected a 2D list at outer list check");
        return -1;
    }

    Py_ssize_t arg_rows = PyList_Size(v);
    if (arg_rows != rows) {
        PyErr_SetString(PyExc_ValueError, "Argument list rows does not match index slice rows");
        return -1;
    }

    for (Py_ssize_t i = 0; i < rows; i++) {
        PyObject *row_obj = PyList_GetItem(v, i);
        if (check_1D(cols, row_obj) != 0) {
            return -1;
        }
    }

    for (Py_ssize_t i = 0; i < rows; i++) {
        PyObject *row_obj = PyList_GetItem(v, i);
        set_1D(mat, row_offset + i, col_offset, cols, row_obj);
    }

    return 0;
}

PyMappingMethods Matrix61c_mapping = {
    NULL,
    (binaryfunc) Matrix61c_subscript,
    (objobjargproc) Matrix61c_set_subscript,
};

/* INSTANCE ATTRIBUTES*/
PyMemberDef Matrix61c_members[] = {
    {
        "shape", T_OBJECT_EX, offsetof(Matrix61c, shape), 0,
        "(rows, cols)"
    },
    {NULL}  /* Sentinel */
};

PyTypeObject Matrix61cType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "numc.Matrix",
    .tp_basicsize = sizeof(Matrix61c),
    .tp_dealloc = (destructor)Matrix61c_dealloc,
    .tp_repr = (reprfunc)Matrix61c_repr,
    .tp_as_number = &Matrix61c_as_number,
    .tp_flags = Py_TPFLAGS_DEFAULT |
    Py_TPFLAGS_BASETYPE,
    .tp_doc = "numc.Matrix objects",
    .tp_methods = Matrix61c_methods,
    .tp_members = Matrix61c_members,
    .tp_as_mapping = &Matrix61c_mapping,
    .tp_init = (initproc)Matrix61c_init,
    .tp_new = Matrix61c_new
};


struct PyModuleDef numcmodule = {
    PyModuleDef_HEAD_INIT,
    "numc",
    "Numc matrix operations",
    -1,
    Matrix61c_class_methods
};

/* Initialize the numc module */
PyMODINIT_FUNC PyInit_numc(void) {
    PyObject* m;

    if (PyType_Ready(&Matrix61cType) < 0)
        return NULL;

    m = PyModule_Create(&numcmodule);
    if (m == NULL)
        return NULL;

    Py_INCREF(&Matrix61cType);
    PyModule_AddObject(m, "Matrix", (PyObject *)&Matrix61cType);
    printf("CS61C Fall 2020 Project 4: numc imported!\n");
    fflush(stdout);
    return m;
}