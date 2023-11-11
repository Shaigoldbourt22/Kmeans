

#define PY_SSIZE_T_CLEAN

#include <Python.h>
#include "spkmeans.h"
#include "spkmeans.c"
#include <unistd.h>


static PyObject *fit(PyObject *self, PyObject *args);
static PyObject *pass_to_choose_algo(PyObject *points_py,PyObject *matrix_py, int goal, int k_is_zero);



static PyMethodDef capiMethods[] = {
        {"fit",
                (PyCFunction) fit,
                     METH_VARARGS,
                PyDoc_STR("choosing algorithm")},
        {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "mykmeanssp",
        NULL,
        -1,
        capiMethods
};

PyMODINIT_FUNC
PyInit_spkmeansmodule(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}


static PyObject *fit(PyObject *self, PyObject *args)
{
    int  goal, k_is_zero;
    PyObject *points_py, *matrix_py;


    if(!PyArg_ParseTuple(args, "OOii",&points_py, &matrix_py, &goal, &k_is_zero)) {
        return NULL;
    }

    if (!PyList_Check(points_py) || !PyList_Check(matrix_py)) {
        return NULL;
    }
    return (PyObject*)Py_BuildValue("O", pass_to_choose_algo(points_py, matrix_py, goal, k_is_zero));
}




static PyObject *pass_to_choose_algo(PyObject *points_py ,PyObject *matrix_py, int goal, int k_is_zero){

    double **points, **mat, **T;
    double *array;
    PyObject *ret_list, *list, *point;
    Py_ssize_t i, j, n, k, d=0, new_k;
    n = PyList_Size(points_py);
    k=n;
    if(k_is_zero!=1)       /*   k != 0   */
        k = PyList_Size(matrix_py);
    points = (double**)malloc(((int)n)*sizeof(double*));
    T = (double**)malloc(((int)n)*sizeof(double*));
    array = (double*) calloc((int)n,sizeof(double));
    if(!points || !array || !T) {
        printf("An Error Has Occurred");
        exit(1);
    }
    mat = (double**)malloc(((int)k) * sizeof(double *));
    if(!mat){
        printf("An Error Has Occurred");
        exit(1);
    }


    for(i=0; i<n; i++){
        point = PyList_GetItem(points_py, i);
        if(((int)i)==0)
            d=PyList_Size(point);
        points[(int)i] = (double*) calloc((int)d, sizeof(double));
        T[(int)i] = (double*) calloc((int)k, sizeof(double));
        if(!points[(int)i] || !T[(int)i]) {
            printf("An Error Has Occurred");
            exit(1);
        }
        for(j=0; j<d; j++)
            points[(int)i][(int)j] = PyFloat_AsDouble(PyList_GetItem(point, j));
    }



    for (i = 0; i < k; i++) {
        point = PyList_GetItem(matrix_py, i);
        mat[(int)i] = (double*) calloc(d, sizeof(double));
        if (!mat[(int)i]) {
            printf("An Error Has Occurred");
            exit(1);
        }
        for (j = 0; j < d; j++)
            mat[(int)i][(int)j] = PyFloat_AsDouble(PyList_GetItem(point, j));
    }


    if(k_is_zero || goal!=1)
        new_k = (Py_ssize_t)choose_algo(goal, points, mat, T, n, d, 0);
    else
        new_k = (Py_ssize_t)choose_algo(goal, points, mat, T, n, d, k);


    ret_list = PyList_New(n);
    if(ret_list == NULL) {
        printf("An Error Has Occurred");
        exit(1);

    }
    for (i = 0; i < n; i++) {
        list = PyList_New(new_k);
        for (j = 0; j < new_k; j++)
            PyList_SetItem(list, j, PyFloat_FromDouble(T[(int)i][(int)j]));
        PyList_SetItem(ret_list, i, list);
    }


  for(i=0; i<n; i++) {
        free(points[(int)i]);
        free(T[(int)i]);
    }
    for(i=0; i<k; i++) {
        free(mat[(int)i]);
    }
    free(points);
    free(mat);
    free(T);

    return ret_list;
}

