// Copyright 2012 Patrick Varilly, Stefano Angioletti-Uberti
//
//   This program is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program.  If not, see <http://www.gnu.org/licenses/>.

/// \file _generic.c
/// \brief Inner loops of DNACC self-consistent thoery
/// \author Patrick Varilly
/// \date Wed 1 Feb 2012

#include "Python.h"
#include "numpy/arrayobject.h"
#include <math.h>

static char py_do_calculate_doc[] =
    "do_calculate( p_free, boltz_bind_csr_indptr, boltz_bind_csr_indices,\n"
    "              boltz_bind_csr_data, prefactor,\n"
    "              maxDelta, maxSteps, weights )\n"
    "\n"
    "Inner loop to recalculate self-consistent values of p_free\n"
    "\n"
    "Roughly equivalent Python code\n"
    "\n"
    "boltz_bind_csr = scipy.sparse.csr_matrix( ... )\n"
    "boltz_bind_csr_indptr = boltz_bind_csr.indptr\n"
    "boltz_bind_csr_indices = boltz_bind_csr.indices\n"
    "boltz_bind_csr_data = boltz_bind_csr.data\n"
    "\n"
    "def do_calculate(p_free, boltz_bind_csr_indptr, boltz_bind_csr_indices,\n"
    "                 boltz_bind_csr_data, prefactor, maxDelta, maxSteps):\n"
    "\n"
    "       delta = 1.0\n"
    "       step = 0\n"
    "       while delta > maxDelta and step < maxSteps:\n"
    "           delta = 0.0\n"
    "           for i in xrange(0,len(p_free)):\n"
    "               oldP = p_free[i]\n"
    "               x = mtx[i,:] * p_free\n"
    "               p_free[i] = weights[i] / (1.0 + prefactor * x)\n"
    "               delta = max( delta, abs( p_free[i] - oldP ) / weights[i] )\n"
    "           step += 1\n"
    "\n"
    "       return delta,step";
static PyObject* 
py_do_calculate( PyObject *self, PyObject *args )
{
    PyArrayObject *array_p_free;
    PyArrayObject *array_boltz_bind_csr_indptr;
    PyArrayObject *array_boltz_bind_csr_indices;
    PyArrayObject *array_boltz_bind_csr_data;
    PyArrayObject *array_weights;
    double prefactor;
    npy_intp N;
    double *p_free, *data, *weights;
    int *indptr, *indices;
    int valid;
    int i, j, step, maxSteps, jidx;
    double delta, x, boltzBindIJ, oldP, thisDelta, maxDelta;

    if( !PyArg_ParseTuple( args, "OOOOddiO:do_calculate",
			   &array_p_free,
			   &array_boltz_bind_csr_indptr,
			   &array_boltz_bind_csr_indices,
			   &array_boltz_bind_csr_data,
			   &prefactor,
			   &maxDelta, &maxSteps,
			   &array_weights ) ) {
	return NULL;
    }

    // Check that all array dimensions make sense
    if( PyArray_NDIM( array_p_free ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_indptr ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_indices ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_data ) != 1 ) return NULL;

    N = PyArray_DIMS( array_p_free )[0];

    if( PyArray_DIMS( array_boltz_bind_csr_indptr )[0] != N+1 ) return NULL;

    // Get at array data
    // NB: There's probably some wild condition I'm not checking for here,
    //   but since this function gets used once in code that I wrote, it's
    //   safe enough.
    p_free = (double*)array_p_free->data;
    indptr = (int*) array_boltz_bind_csr_indptr->data;
    indices = (int*) array_boltz_bind_csr_indices->data;
    data = (double*) array_boltz_bind_csr_data->data;

    if( array_weights == Py_None ) {
	weights = NULL;
    }
    else {
	weights = (double*) array_weights->data;
	if( PyArray_NDIM( array_weights ) != 1 ) return NULL;
	if( PyArray_DIMS( array_weights )[0] != N ) return NULL;
    }

    // CSR arrays are laid out so that column indices for row i are
    // stored in indices[indptr[i]:indptr[i+1]] and their corresponding
    // values are stored in data[indptr[i]:indptr[i+1]].
	
    // Check that all p_free values are sensible
    valid = 1;
    for( i = 0; i < N; i++ )
	if( p_free[i] < 0 || p_free[i] > (weights ? weights[i] : 1.0) )
	    valid = 0;
    if( !valid )
	for( i = 0; i < N; i++ )
	    p_free[i] = (weights ? weights[i] : 1.0);

    // Main loop iteration
    delta = maxDelta + 1.0;
    for( step = 0; step < maxSteps && delta > maxDelta; step++ ) {

	delta = 0.0;
	for( i = 0; i < N; i++ ) {

	    // x = sum_j boltz_bind[i,j] * p_free[j]
	    // p_free[i] = weights[i] / (1.0 + prefactor * x)
	    
	    x = 0.0;
	    for( jidx = indptr[i]; jidx < indptr[i+1]; jidx++ ) {
		j = indices[jidx];
		boltzBindIJ = data[jidx];
		x += boltzBindIJ * p_free[j];
	    }
	    oldP = p_free[i];
	    p_free[i] = (weights ? weights[i] : 1.0)/(1.0 + prefactor * x);

	    thisDelta = fabs( oldP - p_free[i] );
	    if( weights && weights[i] != 0.0 ) thisDelta /= weights[i];
	    if( thisDelta > delta )
		delta = thisDelta;
	}
    }

    return Py_BuildValue( "di", delta, step );
}


static char py_do_add_up_doc[] =
    "do_add_up( p_free, boltz_bind_csr_indptr, boltz_bind_csr_indices,\n"
    "           boltz_bind_csr_data )\n"
    "\n"
    "Compute sum_(i <= j) p_free[i] * boltz_bind[i,j] * p_free[j]\n"
    "\n"
    "Roughly equivalent Python code\n"
    "\n"
    "  np.sum( p_free * (sps.triu(boltz_bind_csr) * p_free) )";
static PyObject* 
py_do_add_up( PyObject *self, PyObject *args )
{
    PyArrayObject *array_p_free;
    PyArrayObject *array_boltz_bind_csr_indptr;
    PyArrayObject *array_boltz_bind_csr_indices;
    PyArrayObject *array_boltz_bind_csr_data;
    npy_intp N;
    double *p_free, *data;
    int *indptr, *indices;
    int i, j, jidx;
    double result, cumSum;
    double boltzBindIJ;

    if( !PyArg_ParseTuple( args, "OOOO:do_add_up",
			   &array_p_free,
			   &array_boltz_bind_csr_indptr,
			   &array_boltz_bind_csr_indices,
			   &array_boltz_bind_csr_data ) ) {
	return NULL;
    }

    // Check that all array dimensions make sense
    if( PyArray_NDIM( array_p_free ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_indptr ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_indices ) != 1 ) return NULL;
    if( PyArray_NDIM( array_boltz_bind_csr_data ) != 1 ) return NULL;

    N = PyArray_DIMS( array_p_free )[0];

    if( PyArray_DIMS( array_boltz_bind_csr_indptr )[0] != N+1 ) return NULL;

    // Get at array data
    // NB: There's probably some wild condition I'm not checking for here,
    //   but since this function gets used once in code that I wrote, it's
    //   safe enough.
    p_free = (double*)array_p_free->data;
    indptr = (int*) array_boltz_bind_csr_indptr->data;
    indices = (int*) array_boltz_bind_csr_indices->data;
    data = (double*) array_boltz_bind_csr_data->data;

    // CSR arrays are laid out so that column indices for row i are
    // stored in indices[indptr[i]:indptr[i+1]] and their corresponding
    // values are stored in data[indptr[i]:indptr[i+1]].

    // Main loop iteration
    result = 0.0;
    for( i = 0; i < N; i++ ) {

	if( p_free[i] == 0.0 ) continue;

	// cumSum = sum_(j >= i) boltz_bind[i,j] * p_free[j]
	cumSum = 0.0;
	for( jidx = indptr[i]; jidx < indptr[i+1]; jidx++ ) {
	    j = indices[jidx];
	    boltzBindIJ = data[jidx];
	    if( i <= j )
		cumSum += boltzBindIJ * p_free[j];
	}

	result += p_free[i] * cumSum;
    }

    return Py_BuildValue( "d", result );
}

static PyMethodDef _innerloopmethods[] = {
    {"do_calculate", py_do_calculate, METH_VARARGS, py_do_calculate_doc},
    {"do_add_up", py_do_add_up, METH_VARARGS, py_do_add_up_doc},
    {NULL, NULL, 0, NULL}
};

void init_generic(void)
{
    Py_InitModule( "_generic", _innerloopmethods );
}
