/* supply some frequently used functions as C extensions. */

#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "_funcs.h"

static PyObject *
priors_theta(PyObject *self, PyObject *args)
{
  double x;

  if (!PyArg_ParseTuple(args, "d", &x))
    return NULL;

  return Py_BuildValue("d", _theta(x));
}

static PyObject *
priors_A_delta_if_cbar_f(PyObject *self, PyObject *args)
{
  double t, delta, cbar, Q;
  int h, k;

  if (!PyArg_ParseTuple(args, "ddddii",
                        &t, &delta, &cbar, &Q, &h, &k))
    return NULL;

  return Py_BuildValue("d", _A_delta_if_cbar_f(t, delta, cbar, Q, h, k));
}

static PyObject *
priors_pr_cbar_A(PyObject *self, PyObject *args)
{
  double cbar;

  if (!PyArg_ParseTuple(args, "d", &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cbar_A(cbar));
}

static PyObject*
priors_pr_cbar_B(PyObject *self, PyObject *args)
{
  double cbar;

  if (!PyArg_ParseTuple(args, "d", &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cbar_B(cbar));
}

static PyObject*
priors_pr_cbar_C(PyObject *self, PyObject *args)
{
  double cbar;

  if (!PyArg_ParseTuple(args, "d", &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cbar_C(cbar));
}

static PyObject *
priors_pr_cn_if_cbar_A(PyObject *self, PyObject *args)
{
  double cn, cbar;

  if (!PyArg_ParseTuple(args, "dd", &cn, &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cn_if_cbar_A(cn, cbar));
}

static PyObject *
priors_pr_cn_if_cbar_B(PyObject *self, PyObject *args)
{
  double cn, cbar;

  if (!PyArg_ParseTuple(args, "dd", &cn, &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cn_if_cbar_B(cn, cbar));
}

static PyObject *
priors_pr_cn_if_cbar_C(PyObject *self, PyObject *args)
{
  double cn, cbar;

  if (!PyArg_ParseTuple(args, "dd", &cn, &cbar))
    return NULL;

  return Py_BuildValue("d", _pr_cn_if_cbar_C(cn, cbar));
}

static PyObject *
priors_fct_sinc(PyObject *self, PyObject *args)
{
  double delta, cbar, Q;
  int h, k;

  if (!PyArg_ParseTuple(args, "dddii", &delta, &cbar, &Q, &h, &k))
    return NULL;

  return Py_BuildValue("d", _fct_sinc(delta, cbar, Q, h, k));
}


static PyMethodDef priorsMethods[] = {
  {"_theta", priors_theta,
   METH_VARARGS, NULL},
  
  {"_A_delta_if_cbar_f", priors_A_delta_if_cbar_f,
   METH_VARARGS, NULL},

  {"_pr_cbar_A", priors_pr_cbar_A,
   METH_VARARGS, NULL},

  {"_pr_cbar_B", priors_pr_cbar_B,
   METH_VARARGS, NULL},

  {"_pr_cbar_C", priors_pr_cbar_C,
   METH_VARARGS, NULL},

  {"_pr_cn_if_cbar_A", priors_pr_cn_if_cbar_A,
   METH_VARARGS, NULL},

  {"_pr_cn_if_cbar_B", priors_pr_cn_if_cbar_B,
   METH_VARARGS, NULL},

  {"_pr_cn_if_cbar_C", priors_pr_cn_if_cbar_C,
   METH_VARARGS, NULL},

  {"_fct_sinc", priors_fct_sinc,
   METH_VARARGS, NULL},

  {NULL, NULL, 0, NULL}
};


static struct PyModuleDef priorsmodule = {
  PyModuleDef_HEAD_INIT,
  "priors",                     /* module name */
  NULL,                         /* module doc strings */
  -1,                           /* size of pre-interpreter state */
  priorsMethods                 /* method table */
};


PyMODINIT_FUNC PyInit_priors(void)
{
  PyObject *module = PyModule_Create(&priorsmodule);

  return module;
}
