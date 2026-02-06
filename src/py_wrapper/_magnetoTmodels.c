#include <Python.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
/* Use same array API symbol as numpy f2py for fortranobject.c compatibility */
#define PY_ARRAY_UNIQUE_SYMBOL _npy_f2py_ARRAY_API
#include <numpy/arrayobject.h>

#define F90WRAP_F_SYMBOL(name) name##_

void f90wrap_abort_(char *message, int len_message)
{
    /* Acquire GIL since we're calling Python C-API from Fortran */
    PyGILState_STATE gstate = PyGILState_Ensure();
    
    if (message == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "f90wrap_abort called");
        PyGILState_Release(gstate);
        return;
    }
    while (len_message > 0 && message[len_message - 1] == ' ') {
        --len_message;
    }
    if (len_message <= 0) {
        PyErr_SetString(PyExc_RuntimeError, "f90wrap_abort called");
        PyGILState_Release(gstate);
        return;
    }
    PyObject* unicode = PyUnicode_FromStringAndSize(message, len_message);
    if (unicode == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "f90wrap_abort called");
        PyGILState_Release(gstate);
        return;
    }
    PyErr_SetObject(PyExc_RuntimeError, unicode);
    Py_DECREF(unicode);
    PyGILState_Release(gstate);
}

void f90wrap_abort__(char *message, int len_message)
{
    f90wrap_abort_(message, len_message);
}

/* External f90wrap helper functions */
extern void F90WRAP_F_SYMBOL(f90wrap_compute__ta16)(int* f90wrap_n0, int* f90wrap_n1, int* f90wrap_n2, int* f90wrap_n3, \
    int* f90wrap_n4, int* f90wrap_n5, int* f90wrap_n6, int* f90wrap_n7, int* f90wrap_n8, int* f90wrap_n9, int* \
    f90wrap_n10, int* f90wrap_n11, int* date, double* data, double* x, double* y, double* z, double* bx, double* by, \
    double* bz, int* nx, int* ny, int* nz);

static PyObject* wrap_compute_ta16(PyObject* self, PyObject* args, PyObject* kwargs)
{
    int f90wrap_n0_val = 0;
    int f90wrap_n1_val = 0;
    int f90wrap_n2_val = 0;
    int f90wrap_n3_val = 0;
    int f90wrap_n4_val = 0;
    int f90wrap_n5_val = 0;
    int f90wrap_n6_val = 0;
    int f90wrap_n7_val = 0;
    int f90wrap_n8_val = 0;
    int f90wrap_n9_val = 0;
    int f90wrap_n10_val = 0;
    int f90wrap_n11_val = 0;
    PyObject* py_date = NULL;
    PyObject* py_data = NULL;
    PyObject* py_x = NULL;
    PyObject* py_y = NULL;
    PyObject* py_z = NULL;
    PyObject* py_bx = NULL;
    PyObject* py_by = NULL;
    PyObject* py_bz = NULL;
    PyObject* py_nx = NULL;
    int nx_val = 0;
    PyArrayObject* nx_scalar_arr = NULL;
    int nx_scalar_copyback = 0;
    int nx_scalar_is_array = 0;
    PyObject* py_ny = NULL;
    int ny_val = 0;
    PyArrayObject* ny_scalar_arr = NULL;
    int ny_scalar_copyback = 0;
    int ny_scalar_is_array = 0;
    PyObject* py_nz = NULL;
    int nz_val = 0;
    PyArrayObject* nz_scalar_arr = NULL;
    int nz_scalar_copyback = 0;
    int nz_scalar_is_array = 0;
    static char *kwlist[] = {"f90wrap_n0", "f90wrap_n1", "f90wrap_n2", "f90wrap_n3", "f90wrap_n4", "f90wrap_n5", \
        "f90wrap_n6", "f90wrap_n7", "f90wrap_n8", "f90wrap_n9", "f90wrap_n10", "f90wrap_n11", "date", "data", "x", "y", "z", \
        "bx", "by", "bz", "nx", "ny", "nz", NULL};
    
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "|iiiiiiiiiiiiOOOOOOOOOOO", kwlist, &f90wrap_n0_val, &f90wrap_n1_val, \
        &f90wrap_n2_val, &f90wrap_n3_val, &f90wrap_n4_val, &f90wrap_n5_val, &f90wrap_n6_val, &f90wrap_n7_val, \
        &f90wrap_n8_val, &f90wrap_n9_val, &f90wrap_n10_val, &f90wrap_n11_val, &py_date, &py_data, &py_x, &py_y, &py_z, \
        &py_bx, &py_by, &py_bz, &py_nx, &py_ny, &py_nz)) {
        return NULL;
    }
    
    PyArrayObject* date_arr = NULL;
    int* date = NULL;
    /* Extract date array data */
    if (!PyArray_Check(py_date)) {
        PyErr_SetString(PyExc_TypeError, "Argument date must be a NumPy array");
        return NULL;
    }
    date_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_date, NPY_INT, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (date_arr == NULL) {
        return NULL;
    }
    date = (int*)PyArray_DATA(date_arr);
    int n0_date = (int)PyArray_DIM(date_arr, 0);
    
    PyArrayObject* data_arr = NULL;
    double* data = NULL;
    /* Extract data array data */
    if (!PyArray_Check(py_data)) {
        PyErr_SetString(PyExc_TypeError, "Argument data must be a NumPy array");
        return NULL;
    }
    data_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_data, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (data_arr == NULL) {
        return NULL;
    }
    data = (double*)PyArray_DATA(data_arr);
    int n0_data = (int)PyArray_DIM(data_arr, 0);
    
    PyArrayObject* x_arr = NULL;
    double* x = NULL;
    /* Extract x array data */
    if (!PyArray_Check(py_x)) {
        PyErr_SetString(PyExc_TypeError, "Argument x must be a NumPy array");
        return NULL;
    }
    x_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_x, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (x_arr == NULL) {
        return NULL;
    }
    x = (double*)PyArray_DATA(x_arr);
    int n0_x = (int)PyArray_DIM(x_arr, 0);
    f90wrap_n0_val = n0_x;
    
    PyArrayObject* y_arr = NULL;
    double* y = NULL;
    /* Extract y array data */
    if (!PyArray_Check(py_y)) {
        PyErr_SetString(PyExc_TypeError, "Argument y must be a NumPy array");
        return NULL;
    }
    y_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_y, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (y_arr == NULL) {
        return NULL;
    }
    y = (double*)PyArray_DATA(y_arr);
    int n0_y = (int)PyArray_DIM(y_arr, 0);
    f90wrap_n1_val = n0_y;
    
    PyArrayObject* z_arr = NULL;
    double* z = NULL;
    /* Extract z array data */
    if (!PyArray_Check(py_z)) {
        PyErr_SetString(PyExc_TypeError, "Argument z must be a NumPy array");
        return NULL;
    }
    z_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_z, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (z_arr == NULL) {
        return NULL;
    }
    z = (double*)PyArray_DATA(z_arr);
    int n0_z = (int)PyArray_DIM(z_arr, 0);
    f90wrap_n2_val = n0_z;
    
    PyArrayObject* bx_arr = NULL;
    PyObject* py_bx_arr = NULL;
    int bx_needs_copyback = 0;
    double* bx = NULL;
    /* Extract bx array data */
    if (!PyArray_Check(py_bx)) {
        PyErr_SetString(PyExc_TypeError, "Argument bx must be a NumPy array");
        return NULL;
    }
    bx_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_bx, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (bx_arr == NULL) {
        return NULL;
    }
    bx = (double*)PyArray_DATA(bx_arr);
    int n0_bx = (int)PyArray_DIM(bx_arr, 0);
    int n1_bx = (int)PyArray_DIM(bx_arr, 1);
    int n2_bx = (int)PyArray_DIM(bx_arr, 2);
    f90wrap_n3_val = n0_bx;
    f90wrap_n4_val = n1_bx;
    f90wrap_n5_val = n2_bx;
    Py_INCREF(py_bx);
    py_bx_arr = py_bx;
    if (PyArray_DATA(bx_arr) != PyArray_DATA((PyArrayObject*)py_bx) || PyArray_TYPE(bx_arr) != \
        PyArray_TYPE((PyArrayObject*)py_bx)) {
        bx_needs_copyback = 1;
    }
    
    PyArrayObject* by_arr = NULL;
    PyObject* py_by_arr = NULL;
    int by_needs_copyback = 0;
    double* by = NULL;
    /* Extract by array data */
    if (!PyArray_Check(py_by)) {
        PyErr_SetString(PyExc_TypeError, "Argument by must be a NumPy array");
        return NULL;
    }
    by_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_by, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (by_arr == NULL) {
        return NULL;
    }
    by = (double*)PyArray_DATA(by_arr);
    int n0_by = (int)PyArray_DIM(by_arr, 0);
    int n1_by = (int)PyArray_DIM(by_arr, 1);
    int n2_by = (int)PyArray_DIM(by_arr, 2);
    f90wrap_n6_val = n0_by;
    f90wrap_n7_val = n1_by;
    f90wrap_n8_val = n2_by;
    Py_INCREF(py_by);
    py_by_arr = py_by;
    if (PyArray_DATA(by_arr) != PyArray_DATA((PyArrayObject*)py_by) || PyArray_TYPE(by_arr) != \
        PyArray_TYPE((PyArrayObject*)py_by)) {
        by_needs_copyback = 1;
    }
    
    PyArrayObject* bz_arr = NULL;
    PyObject* py_bz_arr = NULL;
    int bz_needs_copyback = 0;
    double* bz = NULL;
    /* Extract bz array data */
    if (!PyArray_Check(py_bz)) {
        PyErr_SetString(PyExc_TypeError, "Argument bz must be a NumPy array");
        return NULL;
    }
    bz_arr = (PyArrayObject*)PyArray_FROM_OTF(
        py_bz, NPY_FLOAT64, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
    if (bz_arr == NULL) {
        return NULL;
    }
    bz = (double*)PyArray_DATA(bz_arr);
    int n0_bz = (int)PyArray_DIM(bz_arr, 0);
    int n1_bz = (int)PyArray_DIM(bz_arr, 1);
    int n2_bz = (int)PyArray_DIM(bz_arr, 2);
    f90wrap_n9_val = n0_bz;
    f90wrap_n10_val = n1_bz;
    f90wrap_n11_val = n2_bz;
    Py_INCREF(py_bz);
    py_bz_arr = py_bz;
    if (PyArray_DATA(bz_arr) != PyArray_DATA((PyArrayObject*)py_bz) || PyArray_TYPE(bz_arr) != \
        PyArray_TYPE((PyArrayObject*)py_bz)) {
        bz_needs_copyback = 1;
    }
    
    int* nx = &nx_val;
    if (PyArray_Check(py_nx)) {
        nx_scalar_arr = (PyArrayObject*)PyArray_FROM_OTF(
            py_nx, NPY_INT, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
        if (nx_scalar_arr == NULL) {
            return NULL;
        }
        if (PyArray_SIZE(nx_scalar_arr) != 1) {
            PyErr_SetString(PyExc_ValueError, "Argument nx must have exactly one element");
            Py_DECREF(nx_scalar_arr);
            return NULL;
        }
        nx_scalar_is_array = 1;
        nx = (int*)PyArray_DATA(nx_scalar_arr);
        nx_val = nx[0];
        if (PyArray_DATA(nx_scalar_arr) != PyArray_DATA((PyArrayObject*)py_nx) || PyArray_TYPE(nx_scalar_arr) != \
            PyArray_TYPE((PyArrayObject*)py_nx)) {
            nx_scalar_copyback = 1;
        }
    } else if (PyNumber_Check(py_nx)) {
        nx_val = (int)PyLong_AsLong(py_nx);
        if (PyErr_Occurred()) {
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument nx must be a scalar number or NumPy array");
        return NULL;
    }
    int* ny = &ny_val;
    if (PyArray_Check(py_ny)) {
        ny_scalar_arr = (PyArrayObject*)PyArray_FROM_OTF(
            py_ny, NPY_INT, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
        if (ny_scalar_arr == NULL) {
            return NULL;
        }
        if (PyArray_SIZE(ny_scalar_arr) != 1) {
            PyErr_SetString(PyExc_ValueError, "Argument ny must have exactly one element");
            Py_DECREF(ny_scalar_arr);
            return NULL;
        }
        ny_scalar_is_array = 1;
        ny = (int*)PyArray_DATA(ny_scalar_arr);
        ny_val = ny[0];
        if (PyArray_DATA(ny_scalar_arr) != PyArray_DATA((PyArrayObject*)py_ny) || PyArray_TYPE(ny_scalar_arr) != \
            PyArray_TYPE((PyArrayObject*)py_ny)) {
            ny_scalar_copyback = 1;
        }
    } else if (PyNumber_Check(py_ny)) {
        ny_val = (int)PyLong_AsLong(py_ny);
        if (PyErr_Occurred()) {
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument ny must be a scalar number or NumPy array");
        return NULL;
    }
    int* nz = &nz_val;
    if (PyArray_Check(py_nz)) {
        nz_scalar_arr = (PyArrayObject*)PyArray_FROM_OTF(
            py_nz, NPY_INT, NPY_ARRAY_F_CONTIGUOUS | NPY_ARRAY_FORCECAST);
        if (nz_scalar_arr == NULL) {
            return NULL;
        }
        if (PyArray_SIZE(nz_scalar_arr) != 1) {
            PyErr_SetString(PyExc_ValueError, "Argument nz must have exactly one element");
            Py_DECREF(nz_scalar_arr);
            return NULL;
        }
        nz_scalar_is_array = 1;
        nz = (int*)PyArray_DATA(nz_scalar_arr);
        nz_val = nz[0];
        if (PyArray_DATA(nz_scalar_arr) != PyArray_DATA((PyArrayObject*)py_nz) || PyArray_TYPE(nz_scalar_arr) != \
            PyArray_TYPE((PyArrayObject*)py_nz)) {
            nz_scalar_copyback = 1;
        }
    } else if (PyNumber_Check(py_nz)) {
        nz_val = (int)PyLong_AsLong(py_nz);
        if (PyErr_Occurred()) {
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError, "Argument nz must be a scalar number or NumPy array");
        return NULL;
    }
    /* Call f90wrap helper */
    F90WRAP_F_SYMBOL(f90wrap_compute__ta16)(&f90wrap_n0_val, &f90wrap_n1_val, &f90wrap_n2_val, &f90wrap_n3_val, \
        &f90wrap_n4_val, &f90wrap_n5_val, &f90wrap_n6_val, &f90wrap_n7_val, &f90wrap_n8_val, &f90wrap_n9_val, \
        &f90wrap_n10_val, &f90wrap_n11_val, date, data, x, y, z, bx, by, bz, nx, ny, nz);
    if (PyErr_Occurred()) {
        Py_XDECREF(date_arr);
        Py_XDECREF(data_arr);
        Py_XDECREF(x_arr);
        Py_XDECREF(y_arr);
        Py_XDECREF(z_arr);
        Py_XDECREF(py_bx_arr);
        Py_XDECREF(py_by_arr);
        Py_XDECREF(py_bz_arr);
        return NULL;
    }
    
    if (nx_scalar_is_array) {
        if (nx_scalar_copyback) {
            if (PyArray_CopyInto((PyArrayObject*)py_nx, nx_scalar_arr) < 0) {
                Py_DECREF(nx_scalar_arr);
                return NULL;
            }
        }
        Py_DECREF(nx_scalar_arr);
    }
    if (ny_scalar_is_array) {
        if (ny_scalar_copyback) {
            if (PyArray_CopyInto((PyArrayObject*)py_ny, ny_scalar_arr) < 0) {
                Py_DECREF(ny_scalar_arr);
                return NULL;
            }
        }
        Py_DECREF(ny_scalar_arr);
    }
    if (nz_scalar_is_array) {
        if (nz_scalar_copyback) {
            if (PyArray_CopyInto((PyArrayObject*)py_nz, nz_scalar_arr) < 0) {
                Py_DECREF(nz_scalar_arr);
                return NULL;
            }
        }
        Py_DECREF(nz_scalar_arr);
    }
    Py_DECREF(date_arr);
    Py_DECREF(data_arr);
    Py_DECREF(x_arr);
    Py_DECREF(y_arr);
    Py_DECREF(z_arr);
    if (bx_needs_copyback) {
        if (PyArray_CopyInto((PyArrayObject*)py_bx, bx_arr) < 0) {
            Py_DECREF(bx_arr);
            Py_DECREF(py_bx_arr);
            return NULL;
        }
    }
    Py_DECREF(bx_arr);
    if (by_needs_copyback) {
        if (PyArray_CopyInto((PyArrayObject*)py_by, by_arr) < 0) {
            Py_DECREF(by_arr);
            Py_DECREF(py_by_arr);
            return NULL;
        }
    }
    Py_DECREF(by_arr);
    if (bz_needs_copyback) {
        if (PyArray_CopyInto((PyArrayObject*)py_bz, bz_arr) < 0) {
            Py_DECREF(bz_arr);
            Py_DECREF(py_bz_arr);
            return NULL;
        }
    }
    Py_DECREF(bz_arr);
    /* Build result tuple, filtering out NULL objects */
    int result_count = 0;
    if (py_bx_arr != NULL) result_count++;
    if (py_by_arr != NULL) result_count++;
    if (py_bz_arr != NULL) result_count++;
    if (result_count == 0) {
        Py_RETURN_NONE;
    }
    if (result_count == 1) {
        if (py_bx_arr != NULL) return py_bx_arr;
        if (py_by_arr != NULL) return py_by_arr;
        if (py_bz_arr != NULL) return py_bz_arr;
    }
    PyObject* result_tuple = PyTuple_New(result_count);
    if (result_tuple == NULL) {
        if (py_bx_arr != NULL) Py_DECREF(py_bx_arr);
        if (py_by_arr != NULL) Py_DECREF(py_by_arr);
        if (py_bz_arr != NULL) Py_DECREF(py_bz_arr);
        return NULL;
    }
    int tuple_index = 0;
    if (py_bx_arr != NULL) {
        PyTuple_SET_ITEM(result_tuple, tuple_index++, py_bx_arr);
    }
    if (py_by_arr != NULL) {
        PyTuple_SET_ITEM(result_tuple, tuple_index++, py_by_arr);
    }
    if (py_bz_arr != NULL) {
        PyTuple_SET_ITEM(result_tuple, tuple_index++, py_bz_arr);
    }
    return result_tuple;
}

/* Method table for _magnetoTmodels module */
static PyMethodDef _magnetoTmodels_methods[] = {
    {"f90wrap_compute__ta16", (PyCFunction)wrap_compute_ta16, METH_VARARGS | METH_KEYWORDS, "Wrapper for ta16"},
    {NULL, NULL, 0, NULL}  /* Sentinel */
};

static struct PyModuleDef _magnetoTmodelsmodule = {
    PyModuleDef_HEAD_INIT,
    "magnetoTmodels",
    "Direct-C wrapper for _magnetoTmodels module",
    -1,
    _magnetoTmodels_methods
};

/* Module initialization */
PyMODINIT_FUNC PyInit__magnetoTmodels(void)
{
    import_array();  /* Initialize NumPy */
    PyObject* module = PyModule_Create(&_magnetoTmodelsmodule);
    if (module == NULL) {
        return NULL;
    }
    return module;
}
