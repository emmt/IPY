/*
 * linop.i --
 *
 * General linear operator class for Yorick.
 *
 *-----------------------------------------------------------------------------
 *
 * This file is part of IPY package, "Inverse Problems with Yorick", available
 * at <https://github.com/emmt/IPY>.
 *
 * Copyright (C) 2007-2012, Éric Thiébaut.
 * Copyright (C) 2013-2015, MiTiV project <http://mitiv.univ-lyon1.fr>
 * Copyright (C) 2016, Éric Thiébaut.
 *
 * This file is part of free software IPY: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * IPY is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 *-----------------------------------------------------------------------------
 */

/* Yeti is required. */
if (! is_func(h_new)) include, "yeti.i", 1;

local _linop_identity, _linop_diagonal, _linop_sparse, _linop_full;
local _linop_function, _linop_parametric_function;
local LINOP_DIRECT, LINOP_ADJOINT, LINOP_INVERSE;
local LINOP_INVERSE_ADJOINT, LINOP_ADJOINT_INVERSE;
local LINOP_AUTO, LINOP_IDENTITY, LINOP_SCALAR, LINOP_USERFUNC;
local LINOP_DIAGONAL, LINOP_SPARSE, LINOP_FULL;
local linop_new, linop_apply, _linop_type_table;
/* DOCUMENT op = linop_new();
         or op = linop_new(id);
         or op = linop_new(id, a);
         or op = linop_new(other);
         or op = linop_new(s);
         or op = linop_new(f);
         or op = linop_new(f, p);
         or op(x);
         or op(x, job);
         or linop_apply(op, x);
         or linop_apply(op, x, job);

     The function linop_new() creates a new linear operator object OP which can
     be used as a function (or with the linop_apply() function) to compute the
     dot product of the 'vector' X by the 'matrix' (or its adjoint if JOB is
     LINOP_ADJOINT) which corresponds to the linear operator.  The different
     possibilities to define a linear operator object are described in the
     following.  If keyword VIRTUAL is true, the function implementing the
     linear operator is not directly referenced and the actual function
     defintion is resolved whenever the linear operator is evaluated used (see
     symlink_to_name).

     The function linop_apply() applies the linear operator OP to a 'vector' X
     according to the value of JOB as follows (the symbolic names of job
     constants are indicated in parenthesis and are recommended to improve the
     readability of the code):

        JOB = 0 (LINOP_DIRECT) or unspecified to apply the direct operator,
              1 (LINOP_ADJOINT) to apply the adjoint operator,
              2 (LINOP_INVERSE) to apply the inverse operator,
              3 (LINOP_INVERSE_ADJOINT  or LINOP_ADJOINT_INVERSE)
                to apply the inverse adjoint operator.

     Thanks to hash-table evaluator, it is not necessary to use linop_apply()
     and the usage of the linear operator object is simplified as follows:

         OP(x)        is the same as: linop_apply(OP, x)
         OP(x, job)   is the same as: linop_apply(OP, x, job)


   RECURSION

     If linop_new() is called with a single argument which is already a linear
     operator object, then this object is returned (not a copy).  Hence:

         op = linop_new(other);

     where OTHER is already a linear operator object, simply creates a new
     reference to the object OTHER. (FIXME: we should have the possibility to
     'clone' such an object).


   IDENTITY MATRIX

     A linear operator object implementing the identity operation can be
     defined one of:

         op = linop_new();
         op = linop_new(LINOP_IDENTITY);
         op = linop_new("identity");


   MULTIPLICATION BY A SCALAR

     A linear operator object implementing the multiplication by a scalar can
     be defined by one of:

         op = linop_new(LINOP_SCALAR, a);
         op = linop_new("scalar", a);

     where A gives the scalar parameter of the 'matrix'.


   DIAGONAL MATRIX

     A linear operator object implementing the multiplication by a diagonal
     matrix can be defined by one of:

         op = linop_new(LINOP_DIAGONAL, a);
         op = linop_new("diagonal", a);

     where A gives the coefficients of the diagonal of the 'matrix'; A can be
     multi-dimensional, it must however be conformable with the 'vector' X when
     linop_apply() is used.


   FULL MATRIX

     A linear operator object can be implemented given a full matrix as:

         op = linop_new(LINOP_FULL, a);
         op = linop_new("full", a);

     where A is the array of matrix coefficients which are used as with
     mvmult() function (which to see).


   SPARSE MATRIX

     A linear operator object can be implemented as a sparse matrix:

         op = linop_new(s);

     where S is a sparse matrix object (see sparse_matrix).


   USER DEFINED FUNCTION

     The linear operator object functionalities (dot product with the
     corresponding 'matrix' or its adjoint) can be implemented by a user
     defined function F.  There are two different possibilities depending
     whether or not the function needs aditional data (for instance to store
     the coefficients of the 'matrix').

     If argument P is omitted, the pseudo-code for F must be:

         func f(x, job) {
           if (! job) return A.x;
           else if (job == 1) return A'.x;
           else if (job == 2) return (1/A).x;
           else if (job == 3) return (1/A)'.x;
           error, "unsupported value for JOB";
         }

     where A is the 'matrix' corresponding to the linear operator, where the
     dot and the prime indicate dot product and matrix transposition
     respectively and where (1/A) indicates matrix inverse.  Note that,
     depending on your needs, not all operations must be implemented in the
     function F.  For instance, if only direct and matrix adjoint products
     are implemented, the function can be something like:

         func f(x, job) {
           if (! job) return A.x;
           else if (job == 1) return A'.x;
           error, "unsupported value for JOB";
         }

     If argument P is specified, the pseudo-code for F must be:

         func f(p, x, job) {
           if (! job) return A(p).x;
           else if (job == 1) return A(p)'.x;
           else if (job == 2) return (1/A(p)).x;
           else if (job == 3) return (1/A(p))'.x;
           error, "unsupported value for JOB";
         }

     where A(p) is the 'matrix' which depends on the 'parameters' P.  Note that
     this is purely a notation: P can be anything needed by the user-defined
     operator.


   SEE ALSO:
     is_linop, sparse_matrix, mvmult, h_new, linop_new_fft, linop_new_fftw.
 */
func linop_new(a1, a2, virtual=)
{
  /*
   * Layout of a linear operator object as build by the linop_new function:
   *     op.f   = user-defined function
   *     op.p   = client data for f (or nil)
   *     op.a   = coefficients of the matrix
   *     op.s   = sparse matrix object
   *     op.w   = wrapper function
   */
  t1 = identof(a1);
  t2 = identof(a2);

  if (is_scalar(a1)) {
    if (t1 == Y_STRING) {
      id = _linop_type_table(a1);
    } else if (Y_CHAR <= t1 && t1 <= Y_LONG) {
      id = long(a1);
    }
    if (id == LINOP_IDENTITY) {
      if (t2 != Y_VOID) {
        error, "exceeding argument to define identity operator";
      }
      return _linop_init(h_new(class="linop", type=LINOP_IDENTITY),
                         "_linop_identity");
    } else if (id == LINOP_SCALAR) {
      if (! is_scalar(a2)) {
        error, "non-calar coefficient";
      }
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical scalar coefficient";
      }
      return _linop_init(h_new(class="linop", type=LINOP_SCALAR, a=a2),
                         "_linop_scalar");
    } else if (id == LINOP_DIAGONAL) {
      // FIXME: optimize if all coefficients are zero or one
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical diagonal coefficients";
      }
      return _linop_init(h_new(class="linop", type=LINOP_DIAGONAL, a=a2),
                         "_linop_diagonal");
    } else if (id == LINOP_FULL) {
      if (t2 < Y_CHAR || t2 > Y_COMPLEX) {
        error, "non-numerical matrix coefficients";
      }
      return _linop_init(h_new(class="linop", type=LINOP_FULL, a=a2),
                         "_linop_full");
    } else if (id == LINOP_AUTO) {
      t1 = t2;
      a1 = unref(a2);
      t2 = Y_VOID;
    } else {
      error, "invalid linear operator identifier";
    }
  }

  if (t1 == Y_FUNCTION || t1 == Y_BUILTIN) {
    if (t2 == Y_VOID) {
      return _linop_init(h_new(class="linop", type=LINOP_USERFUNC, f=a1),
                         "_linop_function");
    } else {
      return _linop_init(h_new(class="linop", type=LINOP_USERFUNC, f=a1, p=a2),
                             "_linop_parametric_function");
    }
  } else if (t1 == Y_VOID) {
    if (t2 == Y_VOID) {
      return _linop_init(h_new(class="linop", type=LINOP_IDENTITY),
                         "_linop_identity");
    }
  } else if (t1 == Y_OPAQUE) {
    if (is_sparse_matrix(a1)) {
      if (t2 != Y_VOID) {
        error, "exceeding argument to define sparse linear operator";
      }
      // FIXME: not needed for sparse matrix?
      return _linop_init(h_new(class="linop", type=LINOP_SPARSE, s=a1),
                         "_linop_sparse");
    } else if (is_linop(a1)) {
      if (t2 != Y_VOID) {
        error, "exceeding argument";
      }
      return a1;
    }
  }
  error, "bad argument(s) to define sparse linear operator";
}

/* Job values for linear operators: */
LINOP_DIRECT          = 0;
LINOP_ADJOINT         = 1;
LINOP_INVERSE         = 2;
LINOP_INVERSE_ADJOINT = (LINOP_ADJOINT|LINOP_INVERSE);
LINOP_ADJOINT_INVERSE = (LINOP_ADJOINT|LINOP_INVERSE);

/* Type of linear operators: */
LINOP_AUTO     = 0;
LINOP_IDENTITY = 1;
LINOP_DIAGONAL = 2;
LINOP_SPARSE   = 3;
LINOP_FULL     = 4;
LINOP_SCALAR   = 5;
LINOP_USERFUNC = 6;

/* For compatibility. */
LINOP_TRANSPOSE = LINOP_ADJOINT;
LINOP_INVERSE_TRANSPOSE = LINOP_INVERSE_ADJOINT;
LINOP_TRANSPOSE_INVERSE = LINOP_ADJOINT_INVERSE;

_linop_type_table = h_new(auto=LINOP_AUTO,
                          identity=LINOP_IDENTITY,
                          diagonal=LINOP_DIAGONAL,
                          sparse=LINOP_SPARSE,
                          full=LINOP_FULL,
                          userfunc=LINOP_USERFUNC,
                          scalar=LINOP_SCALAR);

func linop_apply(op, x, job) /* DOCUMENTED */
{
  /* Call the wrapper. */
  return op.w(op, x, job);
}

func _linop_init(op, evalname)
{
  extern virtual;
  h_evaluator, op, evalname;
  return h_set(op, w=(virtual
                        ? symlink_to_name(evalname)
                        : symbol_def(evalname)));
}

func _linop_identity(op, x, job)
{
  return x;
}

func _linop_scalar(op, x, job)
{
  if (! job || job == 1) {
    return op.a*x;
  } else if (job == 2 || job == 3) {
    return (1.0/op.a)*x;
  }
  error, "unexpected JOB";
}

func _linop_diagonal(op, x, job)
{
  if (! job || job == 1) {
    return op.a*x;
  } else if (job == 2 || job == 3) {
    /* Speed-up: compute/get fast matrix inverse. */
    local ainv; eq_nocopy, ainv, op.ainv;
    if (is_void(ainv)) {
      ainv = 1.0/op.a;
      h_set, op, ainv = ainv;
    }
    return ainv*x;
  }
  error, "unexpected JOB";
}

func _linop_sparse(op, x, job)
{
  if (! job || job == 1) {
    return op.s(x, job);
  }
  error, "unsupported value for JOB in sparse linear operator";
}

func _linop_function(op, x, job)
{
  return op.f(x, job);
}

func _linop_parametric_function(op, x, job)
{
  return op.f(op.p, x, job);
}

func _linop_full(op, x, job)
{
  if (! job || job == 1) {
    return mvmult(op.a, x, job);
  }
  error, "unsupported value for JOB for full matrix linear operator";
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

local linop_is_identity, linop_is_scalar, linop_is_diagonal, linop_is_sparse;
local linop_is_full, linop_is_userfunc, linop_is_virtual;
func is_linop(op)
/* DOCUMENT is_linop(op);
         or linop_is_identity(op);
         or linop_is_scalar(op);
         or linop_is_diagonal(op);
         or linop_is_sparse(op);
         or linop_is_full(op);
         or linop_is_userfunc(op);
         or linop_is_virtual(op);

     The function is_linop() checks whether object OP is a linear operator.

     The function linop_is_identity() checks whether object OP is a linear
     operator which implements the identity.

     The function linop_is_scalar() checks whether object OP is a linear
     operator which implements a multiplication by a scalar.

     The function linop_is_diagonal() checks whether object OP is a linear
     operator which implements a diagonal operator.

     The function linop_is_sparse() checks whether object OP is a linear
     operator which implements a sparse operator.

     The function linop_is_full() checks whether object OP is a linear operator
     which is a full matrix.

     The function linop_is_userfunc() checks whether object OP is a linear
     operator which is implemented by a user-defined function.

     The function linop_is_virtual() checks whether object OP is a linear
     operator implemented by a virtual evaluator.


   SEE ALSO: linop_new, linop_get_type, linop_get_coefs.
 */
{
  return (is_sparse_matrix(op) ||
          (is_hash(op) && is_string((class = op.class))
           && is_scalar(class) && class == "linop"));
}

func linop_is_virtual(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && is_symlink(op.w));
}

func linop_is_identity(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && op.type == LINOP_IDENTITY);
}

func linop_is_scalar(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && op.type == LINOP_SCALAR);
}

func linop_is_diagonal(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && op.type == LINOP_DIAGONAL);
}

func linop_is_sparse(op) /* DOCUMENTED */
{
  if (is_hash(op)) {
    return (op.class == "linop" && op.type == LINOP_SPARSE);
  } else {
    return is_sparse_matrix(op);
  }
}

func linop_is_full(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && op.type == LINOP_FULL);
}

func linop_is_userfunc(op) /* DOCUMENTED */
{
  return (is_hash(op) && op.class == "linop" && op.type == LINOP_USERFUNC);
}

func linop_get_type(op)
/* DOCUMENT linop_get_type(op);

     The function linop_get_type() returns the type of the linear operator OP
     (one of: LINOP_IDENTITY, LINOP_SCALAR, LINOP_DIAGONAL, LINOP_SPARSE,
     LINOP_FULL or LINOP_USERFUNC) or nil if OP is not a linear operator.

   SEE ALSO: linop_new, is_linop.
 */
{
  if (is_hash(op) && op.class == "linop") {
    return op.type;
  } else if (is_sparse_matrix(op)) {
    return LINOP_SPARSE;
  }
}

func linop_get_coefs(op)
/* DOCUMENT linop_get_coefs(op);

     The function linop_get_coefs() returns the coefficient(s) of the linear
     operator OP is it is of type LINOP_SCALAR, LINOP_DIAGONAL, LINOP_SPARSE,
     or LINOP_FULL; nothing is returned otherwise.

   SEE ALSO: linop_new, is_linop, linop_get_type, linop_make_matrix.
 */
{
  if (is_hash(op) && op.class == "linop" &&
      ((type = op.type) == LINOP_SCALAR ||
       type == LINOP_DIAGONAL ||
       type == LINOP_FULL)) {
    return op.a;
  } else if (is_sparse_matrix(op)) {
    return op.coefs;
  }
}

func linop_make_matrix(op, x, job, multi=)
/* DOCUMENT a = linop_make_matrix(op, x);
         or a = linop_make_matrix(op, x, job);

    Use linear operator OP (see linop_new) with input "vectors" of same data
    type (real or complex) and dimension list as X to build a "matrix" A with
    the same coefficients as the linear operator OP(x, JOB) -- see linop_new
    for the meaning of optional argument JOB.  The result A is always a regular
    array.  By default, A is a real M-by-N array where:

        M =   numberof(Y)  if Y is real,
            2*numberof(Y)  if Y is complex,

    where Y = OP(X, JOB), and

        N =   numberof(X)  if X is real,
            2*numberof(X)  if X is complex.

    If keyword MULTI is true, the dimension list of A is 2, if Y is complex,
    followed by dimsof(Y), followed by 2, if X is complex, followed by
    dimsof(X).


   SEE ALSO: mvmult.
 */
{
  xident = identof(x);
  if (xident == Y_COMPLEX) {
    xcast = linop_cast_real_as_complex;
    xdims = make_dimlist(2, dimsof(x));
    n = 2*numberof(x);
  } else if (xident <= Y_DOUBLE) {
    xcast = double;
    xdims = dimsof(x);
    n = numberof(x);
  } else {
    error, "invalid data type for X";
  }
  x = array(double, xdims);
  for (j = 1; j <= n; ++j) {
    x(j) = 1.0;
    if (j == 1) {
      y = op(xcast(x), job);
      yident = identof(y);
      if (yident == Y_COMPLEX) {
        ycast = linop_cast_complex_as_real;
      } else if (yident <= Y_DOUBLE) {
        ycast = double;
      } else {
        error, "invalid data type for Y";
      }
      y = ycast(y);
      ydims = dimsof(y);
      m = numberof(y);
      a = array(double, m, n);
      a(,j) = y(*);
    } else {
      a(,j) = ycast(op(xcast(x), job))(*);
    }
    x(j) = 0.0;
  }
  if (multi) {
    return linop_reshape(unref(a), ydims, xdims);
  }
  return a;
}

local linop_cast_real_as_complex, linop_cast_complex_as_real;
/* DOCUMENT z = linop_cast_real_as_complex(x);
       -or- x = linop_cast_complex_as_real(z);

     The first function converts a 2-by-any real array X into a complex array Z
     such that:

        Z.re = X(1,..)
        Z.im = X(2,..)

     The second function does the inverse operation.

   SEE ALSO: reshape, linop_reshape.
 */
func linop_cast_complex_as_real(z) /* DOCUMENTED */
{
  if (! is_complex(z)) error, "expecting complex argument";
  reshape, z, &z, double, 2, dimsof(z);
  return z;
}
func linop_cast_real_as_complex(r) /* DOCUMENTED */
{
  if (identof(r) > Y_DOUBLE || (n = numberof((dims = dimsof(r)))) < 2
      || dims(2) != 2) {
    error, "expecting integer or real array with leading dimension equals to 2";
  }
  dims = dims(2:);
  dims(1) = n - 2;
  reshape, r, &double(r), complex, dims;
  return r;
}

func linop_reshape(a, type_or_dims, ..)
/* DOCUMENT linop_reshape(a, type, dim1, dim2, ...);
       -or- linop_reshape(a, dim1, dim2, ...);

     Make input array A into an array with dimension list DIM1, DIM2, ...  If
     TYPE is specified, the result will be of that data type
     _w_i_t_h_o_u_t___c_o_n_v_e_r_s_i_o_n_.

   SEE ALSO:
     make_dimlist, reshape, linop_cast_real_as_complex,
     linop_cast_complex_as_real.
 */
{
  local ref;
  if (identof(type_or_dims) == Y_STRUCTDEF) {
    type = type_or_dims;
    dims = [0];
  } else {
    type = structof(a);
    dims = type_or_dims;
    make_dimlist, dims;
  }
  while (more_args()) {
    make_dimlist, dims, next_arg();
  }
  number = 1;
  for (k = dims(1) + 1; k >= 2; --k) {
    number *= dims(k);
  }
  if (sizeof(a) != sizeof(type)*number) {
    error, "size mismatch";
  }
  reshape, ref, &a, type, dims;
  return ref;
}

/*---------------------------------------------------------------------------*/
/* DEPRECATED STUFF */

func linop_new_fftw(nil, dims=, measure=, real=)
/* DOCUMENT op = linop_new_fftw(...)

     Return a new new linear operator object which can be used to compute FFT
     by means of FFTW, the "fastest FFT in the world".

     Keyword DIMS can be used to pre-specify the dimension list of the arrays
     to be transformed.  If left unspecified, the actual dimension list will be
     initialized the first time the linear operator is applied.  In any cases,
     a given operator can only be used onto arrays with same dimension lists.

     Keywords REAL and MEASURE have the same meaning as for fftw_plan (which to
     see).  Note that fftw_plan is only called as needed and cached into OP to
     save computation time.  Also note that if you use REAL=1, you must
     correctly initialize the dimension list of array to be transformed either
     when linop_new_fftw is called or by computing the first FFTW with JOB=0 or
     3 (*not* 1 or 2).

     Examples:

         OP = linop_new_fftw();
         linop_apply(OP, x)    // compute FFT of X
         linop_apply(OP, x, 0) // idem
         linop_apply(OP, x, 1) // apply adjoint of FFT to X
         linop_apply(OP, x, 2) // compute inverse FFT of X
         linop_apply(OP, x, 3) // apply adjoint of inverse FFT to X

         OP.nevals = number of FFT computed so far by OP

     If you have a recent version of Yeti (which implements hash-table
     evaluators), then it is not necessary to use linop_apply and the usage of
     the linear operator object is simplified as follows:

         OP(x)        is the same as: linop_apply(OP, x)
         OP(x, job)   is the same as: linop_apply(OP, x, job)


   SEE ALSO: linop_new, fftw, fftw_plan, linop_new_fft.
 */
{
  if (! is_func(fftw_plan)) {
    include, "yeti_fftw.i", 1;
  }
  if (! is_void(nil)) error, "no non-keyword argument allowed";
  if (is_void(dims)) {
    state = 0;
  } else {
    state = 1;
    for (k=numberof(dims), number=1; k >= 2; --k) {
      number *= dims(k);
    }
    scl = (1.0/number);
  }
  op = h_new(class="linop", type=LINOP_USERFUNC,
              w=_linop_fftw_wrapper, real=(real ? 1n : 0n),
              dims=dims, scl=scl, measure=measure, nevals=0, state=state);
  h_evaluator, op, "_linop_fftw_wrapper";
  return op;
}

func _linop_fftw_wrapper(op, x, job)
{
  if (! job || job == 3) {
    /* forward transform or backward conjugate transpose */
    if (! ((state = op.state) & 2)) {
      /* compute forward FFTW plan */
      if (! (state & 1)) {
        h_set, op, dims=dimsof(x), scl=(1.0/numberof(x));
      }
      h_set, op, state=(state |= 3),
        fwd=fftw_plan(op.dims, +1, real=op.real, measure=op.measure);
    }
    z = fftw(x, op.fwd);
    h_set, op, nevals = op.nevals + 1;
    return (job == 3 ? op.scl*z : z);
  } else if (job == 1 || job == 2) {
    /* forward conjugate transpose or backward transform */
    if (! ((state = op.state) & 4)) {
      /* compute backward FFTW plan */
      if (! (state & 1)) {
        if (op.real) {
          error, "you must initialize dimension list first (see doc)";
        }
        h_set, op, dims=dimsof(x), scl=(1.0/numberof(x));
      }
      h_set, op, state=(state |= 5),
        bck=fftw_plan(op.dims, -1, real=op.real, measure=op.measure);
    }
    z = fftw(x, op.bck);
    h_set, op, nevals = op.nevals + 1;
    return (job == 2 ? op.scl*z : z);
  }
  error, "unsupported value for JOB in FFTW linear operator";
}

func linop_new_fft(dims, ldir, rdir, real=)
/* DOCUMENT op = linop_new_fft(dims, ldir, rdir)

     Return a new new linear operator object which can be used to compute FFT
     by means of Swarztrauber's FFT.  DIMS is the dimension list of the arrays
     to be transformed and optional arguments LDIR and RDIR indicate the
     dimensions to transform and in which directions (see fft and fft_setup for
     more detailed explanations).  The returned operator can only be used onto
     arrays with same dimension lists.

     Keyword REAL can be set true to specify a real to complex transform.

     The FFT linear operator is more flexible than the FFTW one (can transform
     for only a subset of the dimensions and with different directions) but is
     slower.  Otherwise the two should behave the same and you can see the
     documentation of linop_new_fftw for examples.


   SEE ALSO: linop_new, fft, linop_new_fftw.
 */
{
  real = (real ? 1n : 0n);
  if (is_void(dims)) {
    dims = [0];
  } else if (! is_integer(dims) || ! is_vector(dims) ||
             numberof(dims) != dims(1) + 1 || min(dims) <= 0) {
    error, "invalid dimension list";
  }
  ndims = numberof(dims) - 1;
  dims = long(dims);
  ltyp = _linop_fft_get_dir(ldir);
  rtyp = _linop_fft_get_dir(rdir);
  if (! ltyp && ! rtyp) {
    ltyp = 1;
    ldir = 1;
  } else if (ltyp < 0 || rtyp < 0) {
    error, "bad FFT directions";
  }
  llen = numberof(ldir);
  rlen = numberof(rdir);
  if (llen + rlen > ndims) {
    error, "more FFT directions than number of dimensions";
  }
  if (noneof(ltyp) && noneof(rdir)) {
    wrapper = _linop_identity;
    ndirs = 0;
    scale = 1.0;
    setup = list = length = dirs = top = [];
  } else {
    /* compute stride and number of elements */
    stride = array(long, ndims);
    length = dims(2:0);
    stride(1) = 1;
    for (j = 1; j < ndims; ++j) {
      stride(j + 1) = stride(j)*length(j);
    }
    number = stride(ndims)*length(ndims);

    /* select directions of transform */
    dirs = array(long, ndims);
    if (ltyp == 1 && ! rtyp) {
      dirs(*) = ldir;
    } else {
      if (llen) dirs(1:llen) = ldir;
      if (rlen) dirs(1-rlen:0) = rdir;
    }
    list = where(dirs);
    ndirs = numberof(list);
    length = length(list);
    stride = stride(list);

    /* compute FFT workspaces */
    top = number/(stride*length);
    number = 1;
    setup = array(pointer, ndirs);
    for (j = 1 ; j <= ndirs; ++j) {
      len = length(j);
      number *= len;
      if (j > 1 && (k = where(len == length)(1)) < j) {
        setup(j) = setup(k);
      } else {
        ws = array(double, 6*len + 15);
        fft_init, len, ws;
        setup(j) = &ws;
      }
    }
    scale = 1.0/number; /* scale for inverse FFT */
    wrapper = _linop_fft_wrapper;
  }
  op = h_new(class="linop", type=LINOP_USERFUNC,
              w=_linop_fft_wrapper, dims=dims, nevals=0, real=real,
              scale=scale, list=list, ndirs=ndirs, dirs=dirs,
              stride=stride, length=length, top=top, setup=setup);
  h_evaluator, op, "_linop_fft_wrapper";
  return op;
}

func _linop_fft_get_dir(dir)
{
  if (is_void(dir)) return 0;
  if (is_integer(dir) && min(dir) >= -1 && max(dir) <= +1) {
    if (is_scalar(dir)) return 1;
    if (is_vector(dir)) return 2;
  }
  return -1;
}

func _linop_fft_wrapper(op, x, job)
{
  local dims, dirs, setup, length, stride, top;
  if (! job || job == 3) {
    real = 0n;
  } else if (job == 1 || job == 2) {
    real = op.real;
  } else {
    error, "unsupported value for JOB in FFT linear operator";
  }
  if ((type = identof(x)) > Y_COMPLEX) {
    error, "non-numerical argument";
  }
  eq_nocopy, dims, op.dims;
  if (numberof((xdims = dimsof(x))) != numberof(dims) ||
      anyof(xdims != dims)) {
    error, "incompatible dimensions of argument";
  }
  h_set, op, nevals = (op.nevals + 1);
  if (! (ndirs = op.ndirs)) {
    if (real) {
      if (type == Y_DOUBLE) {
        x = x; /* make a copy */
        return x;
      }
      return double(x);
    } else {
      if (type == Y_COMPLEX) {
        x = x; /* make a copy */
        return x;
      }
      return complex(x);
    }
  }
  if (type == Y_COMPLEX) {
    x = x; /* make a copy for in-place FFT */
  } else {
    x = complex(x);
  }

  /* do the requested transforms in-place */
  if (job == 1 || job == 2) {
    dirs = -op.dirs;
  } else {
    eq_nocopy, dirs, op.dirs;
  }
  eq_nocopy, setup, op.setup;
  eq_nocopy, length, op.length;
  eq_nocopy, stride, op.stride;
  eq_nocopy, top, op.top;
  for (j = 1; j <= ndirs; ++j) {
    fft_raw, dirs(j), x, stride(j), length(j), top(j), setup(j);
  }
  if (real) {
    x = double(x);
  }
  return ((job == 3 || job == 2) ? op.scale*x : x);
}
