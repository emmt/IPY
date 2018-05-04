/*
 * ipy.i --
 *
 * Implement fundamental tools for IPY package, "Inverse Problems with
 * Yorick".
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2013-2015: Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

if (is_func(plug_in)) plug_in, "ipy";

extern ipy_norm_1;
extern ipy_norm_2;
extern ipy_norm_inf;
/* DOCUMENT ipy_norm_1(x);
         or ipy_norm_2(x);
         or ipy_norm_inf(x);

     Compute the L1, L2 (Euclidean) or infinite norms of X.  These are defined
     as:

         ipy_norm_1(x) = sum(abs(x))
         ipy_norm_2(x) = sqrt(sum(abs(x)^2))
         ipy_norm_inf(x) = max(abs(x))

     the returned value is always a double.

   SEE ALSO: ipy_inner. */

extern ipy_inner;
/* DOCUMENT ipy_inner(x, y);
         or ipy_inner(w, x, y);
     Compute the (weighted) inner product of X  and Y.  This is also called the
     dot product.  The result  is the same as sum(X*Y) or  sum(W*X*Y) but it is
     more than twice as fast and  does not require additional memory (except if
     conversions are needed).  The returned value is always a double.

   SEE ALSO: ipy_norm_2. */

extern ipy_combine;
/* DOCUMENT ipy_combine(a1, x1, a2, x2, ...);
         or ipy_combine(dst, a1, x1, a2, x2, ...);
         or ipy_combine, dst, a1, x1, a2, x2, ...;

     Compute the  linear combination:  A1*X1 +  A2*X2 + ...   where the  An are
     scalar  weights  and  the  Xn  are "vectors"  which  must  have  the  same
     dimensions.

     When called with  an odd number of arguments, the  first one (denoted DST)
     must be a simple variable reference which is used to store the result.  As
     far as possible,  the value of DST will be  re-used to avoid re-allocating
     memory.  When called as a subroutine, the destination variable DST must be
     specified.   When  called  as  a   function,  the  result  of  the  linear
     combination is returned.

     Type promotion is  based on the type  of the operands (not on  that of the
     weights) and the result  is of type double if any operand  is of type long
     or double; the result is of type float otherwise.

     Currently, the routine  is implemented for float or double  arrays and for
     up  to 5  operands.  Compare  to equivalent  Yorick code,  the speedup  is
     better than twice the number of operands.

  SEE ALSO: ipy_scale.
*/

func ipy_scale(alpha, x) { return (alpha != 1.0 ? ipy_combine(alpha, x) : x); }
/* DOCUMENT ipy_scale(alpha, x)
     This function returns X multiplied by scalar ALPHA.

   SEE ALSO: ipy_identity, ipy_combine.
 */

/*---------------------------------------------------------------------------*/
/* LINEAR ALGEBRA AND OPERATORS */

func ipy_new_weighting_operator(w)
/* DOCUMENT ipy_new_weighting_operator(w);
     Create a weighting operator.
   SEE ALSO: ipy_new_diagonal_operator.
 */
{
  if (min(w) < 0.0) error, "weights must be nonnegative";
  return ipy_new_diagonal_operator(w);
}

local ipy_new_diagonal_operator;
local ipy_is_diagonal_operator;
local ipy_get_diagonal_of_diagonal_operator;
/* DOCUMENT op = ipy_new_diagonal_operator(w);
         or ipy_is_diagonal_operator(op);
         or w = ipy_get_diagonal_of_diagonal_operator(op);

     Create a diagonal operator, check whether a object is a diagonal operator
     or retrieve the diagonal elements of a diagonal operator.  Argument W
     gives the diagonal elements of the operator.

   SEE ALSO: ipy_new_weighting_operator.
 */

func ipy_new_diagonal_operator(w)
{
  if (! is_array(w)) error, "invalid argument";
  return closure("_ipy_apply_diagonal_operator", w);
}

func _ipy_apply_diagonal_operator(w, x, job) { return w*x; }

func ipy_is_diagonal_operator(op)
{
  return (is_func(op) == 5n &&
          op.function_name == "_ipy_apply_diagonal_operator");
}

func ipy_get_diagonal_of_diagonal_operator(op) { return op.data; }


func ipy_identity(x, job) { return x; }
/* DOCUMENT ipy_identity(x);
         or ipy_identity(x, job);
     This function returns X.  The prototype of this function makes it
     elligible as a linear operator (LinOp).

   SEE ALSO: ipy_scale, linop_new.
 */

func ipy_mirror(a, job)
/* DOCUMENT ipy_mirror(a);
         or ipy_mirror(a, job);
     Returns array A mirrored along all its dimensions.  If JOB is true,
     returns the transpose/inverse of the operation.
   SEE ALSO: linop.
 */
{
  type = structof(a);
  dims = dimsof(a);
  rank = dims(1);
  if (job) {
    /* transpose operation */
    if (rank >= 1) {
      if (anyof(dims(2:0)&1)) {
        error, "bad dimension(s)";
      }
      r1 = 1:dims(2)/2;
      if (rank >= 2) {
        r2 = 1:dims(3)/2;
        if (rank >= 3) {
          r3 = 1:dims(4)/2;
          if (rank >= 4) {
            r4 = 1:dims(5)/2;
            if (rank >= 5) {
              r5 = 1:dims(6)/2;
              if (rank >= 6) {
                r6 = 1:dims(7)/2;
                if (rank >= 7) {
                  r7 = 1:dims(8)/2;
                  if (rank >= 8) {
                    r8 = 1:dims(9)/2;
                    if (rank >= 9) {
                      r9 = 1:dims(10)/2;
                      if (rank >= 10) {
                        r10 = 1:dims(11)/2;
                        if (rank >= 11) {
                          error, "too many dimensions";
                        }
                        return a(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10);
                      }
                      return a(r1,r2,r3,r4,r5,r6,r7,r8,r9);
                    }
                    return a(r1,r2,r3,r4,r5,r6,r7,r8);
                  }
                  return a(r1,r2,r3,r4,r5,r6,r7);
                }
                return a(r1,r2,r3,r4,r5,r6);
              }
              return a(r1,r2,r3,r4,r5);
            }
            return a(r1,r2,r3,r4);
          }
          return a(r1,r2,r3);
        }
        return a(r1,r2);
      }
      return a(r1);
    }
  } else {
    /* direct operation */
    if (rank >= 1) {
      if (rank >= 2) {
        if (rank >= 3) {
          if (rank >= 4) {
            if (rank >= 5) {
              if (rank >= 6) {
                if (rank >= 7) {
                  if (rank >= 8) {
                    if (rank >= 9) {
                      if (rank >= 10) {
                        if (rank >= 11) {
                          error, "too many dimensions";
                        }
                        n = dims(11);
                        dims(11) *= 2;
                        b = array(type, dims);
                        b(,,,,,,,,,1:n, ..) = a;
                        b(,,,,,,,,,n+1:0, ..) = a(,,,,,,,,,::-1,..);
                        eq_nocopy, a, b;
                      }
                      n = dims(10);
                      dims(10) *= 2;
                      b = array(type, dims);
                      b(,,,,,,,,1:n, ..) = a;
                      b(,,,,,,,,n+1:0, ..) = a(,,,,,,,,::-1,..);
                      eq_nocopy, a, b;
                    }
                    n = dims(9);
                    dims(9) *= 2;
                    b = array(type, dims);
                    b(,,,,,,,1:n, ..) = a;
                    b(,,,,,,,n+1:0, ..) = a(,,,,,,,::-1,..);
                    eq_nocopy, a, b;
                  }
                  n = dims(8);
                  dims(8) *= 2;
                  b = array(type, dims);
                  b(,,,,,,1:n, ..) = a;
                  b(,,,,,,n+1:0, ..) = a(,,,,,,::-1,..);
                  eq_nocopy, a, b;
                }
                n = dims(7);
                dims(7) *= 2;
                b = array(type, dims);
                b(,,,,,1:n, ..) = a;
                b(,,,,,n+1:0, ..) = a(,,,,,::-1,..);
                eq_nocopy, a, b;
              }
              n = dims(6);
              dims(6) *= 2;
              b = array(type, dims);
              b(,,,,1:n, ..) = a;
              b(,,,,n+1:0, ..) = a(,,,,::-1,..);
              eq_nocopy, a, b;
            }
            n = dims(5);
            dims(5) *= 2;
            b = array(type, dims);
            b(,,,1:n, ..) = a;
            b(,,,n+1:0, ..) = a(,,,::-1,..);
            eq_nocopy, a, b;
          }
          n = dims(4);
          dims(4) *= 2;
          b = array(type, dims);
          b(,,1:n, ..) = a;
          b(,,n+1:0, ..) = a(,,::-1,..);
          eq_nocopy, a, b;
        }
        n = dims(3);
        dims(3) *= 2;
        b = array(type, dims);
        b(,1:n, ..) = a;
        b(,n+1:0, ..) = a(,::-1,..);
        eq_nocopy, a, b;
      }
      n = dims(2);
      dims(2) *= 2;
      b = array(type, dims);
      b(1:n, ..) = a;
      b(n+1:0, ..) = a(::-1,..);
      eq_nocopy, a, b;
    }
  }
  return a;
}

func ipy_split_index(idx, dims, zero_based)
/* DOCUMENT coord = ipy_split_index(idx, dims);
         or coord = ipy_split_index(idx, dims, zero_based);

     Split index  IDX into  coordinates according to  the dimension  list DIMS.
     The  result  is  RANK-by-dimsof(IDX).    By  default,  1-based  index  and
     coordinates  are  assumed  (as  in Yorick);  if  third  optional  argument
     ZERO_BASED is true, then 0-based index and coordinates are assumed.

     For example, to get the coordinates of the maxima of array A:

        coord = ipy_split_index(where(a == max(a)), dimsof(a));

   SEE ALSO:
 */
{
  rank = dims(1);
  offset = array(long, rank, dimsof(idx));
  if (! zero_based) {
    idx -= 1;
  }
  for (d = 1; d < rank; ++d) {
    len = dims(1 + d);
    offset(d,..) = idx%len;
    idx /= len;
  }
  offset(0,..) = idx;
  return (zero_based ? offset : offset + 1);
}

/*---------------------------------------------------------------------------*/
/* DECONVOLUTION TOOLS */

func ipy_recenter_psf(psf, pos)
/* DOCUMENT ipy_recenter_psf(psf, pos);

     This function returns PSF with its central element rolled at the origin of
     coordinates for the FFT.  The position  of the central element is given by
     POS.  If POS  is unspecified or if POS="fft", the  PSF is assumed centered
     at the origin of FFT coordinates (i.e. at the first element of the array);
     if POS="max",  the center  is at  the (first)  peak value  of the  PSF; if
     POS="min", the  center is  at the  (first) smallest value  of the  PSF; if
     POS="avg", the  center is at the  barycentre of the PSF;  otherwise POS is
     interpreted as  the index  of the  central element  in PSF  (POS can  be a
     scalar or a vector of indices to  give the coordinates of the center along
     all the dimensions of PSF).

   SEE ALSO roll, ipy_centroid.
*/
{
  if (is_void(pos)) return psf;
  scl = is_scalar(pos);
  if (scl && is_string(pos)) {
    if (pos == "fft") {
      return psf;
    } else if (pos == "max") {
      pos = psf(*)(mxx);
    } else if (pos == "avg") {
      pos = lround(1.0 + ipy_centroid(psf));
      scl = 0n;
    } else if (pos == "min") {
      pos = psf(*)(mnx);
    } else {
      error, "bad value for PSF central position";
    }
  } else if (! is_integer(pos)) {
    error, "PSF central position must be a string or an integer";
  }
  ntot = numberof(psf);
  dims = dimsof(psf);
  rank = dims(1);
  dims = (rank >= 1 ? dims(2:0) : []);
  if (scl) {
    if (pos <= 0) pos += ntot;
    if (pos < 1 || pos > ntot) {
      error, "out of range index for PSF central position";
    }
    off = array(long, rank);
    pos = pos - 1; /* Yorick has 1-based indexes, also make sure POS is a long
                      integer */
    for (k = 1; k <= rank; ++k) {
      len = dims(k);
      off(k) = pos % len;
      pos /= len;
    }
  } else {
    if (numberof(pos) != rank) {
      error, ("index for PSF central position must be a scalar " +
              "or have the same size as PSF rank");
    }
    if (min(pos) < 1 || anyof(pos > dims)) {
      error, "out of range index for PSF central position";
    }
    off = pos - 1;
  }
  return (anyof(off) ? roll(psf, - off) : psf);
}

func ipy_centroid(a)
/* DOCUMENT c = ipy_centroid(a);

     This function returns the barycentric coordinates of array A.  Coordinates
     are  0-based and  may  be fractional.   Array  A may  have  any number  of
     dimensions and its  sum should be non-zero (otherwise  all coordinates are
     zero).

   SEE ALSO: ipy_recenter_psf.
 */
{
  dims = dimsof(a);
  rank = dims(1);
  if (rank < 1) {
    return 0.0;
  }
  c = array(double, rank);
  q = sum(a);
  if (q != 0.0) {
    // FIXME: would be much faster if we compute partial sums
    q = 1.0/q;
    for (k = 1; k <= rank; ++k) {
      lenght = dims(k+1);
      if (lenght > 1) {
        d = array(1, k+1);
        d(1) = k;
        d(0) = lenght;
        (x = array(double, d))(*) = indgen(0:lenght-1);
        c(k) = q*sum(a*x);
      }
    }
  }
  return c;
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

func ipy_include(feature, file)
/* DOCUMENT ipy_include, feature, file;
     If FEATURE is true, does nothing; else, include FILE immediately.
     For instance:
        ipy_include, is_func(h_new), "yeti.i";

   SEE ALSO: include.
 */
{
  if (! feature) include, file, 1;
}

func ipy_rms(a, b)
/* DOCUMENT ipy_rms(a,b);
     Returns the root mean squared (RMS) value of the difference between A and
     B; that is: sqrt(avg((A - B)^2));
   SEE ALSO: sqrt.
 */
{
  r = a - b;
  return sqrt(avg(r*r));
}

local ipy_cast_real_as_complex, ipy_cast_complex_as_real;
/* DOCUMENT z = ipy_cast_real_as_complex(x);
       -or- x = ipy_cast_complex_as_real(z);

     The first function converts a 2-by-any real array X into a complex
     array Z such that:

        Z.re = X(1,..)
        Z.im = X(2,..)

     The second function does the inverse operation.

   SEE ALSO: reshape, ipy_reshape.
 */
func ipy_cast_complex_as_real(z) /* DOCUMENTED */
{
  if (! is_complex(z)) error, "expecting complex argument";
  reshape, z, &z, double, 2, dimsof(z);
  return z;
}
func ipy_cast_real_as_complex(r) /* DOCUMENTED */
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

func ipy_reshape(a, type_or_dims, ..)
/* DOCUMENT ipy_reshape(a, type, dim1, dim2, ...);
       -or- ipy_reshape(a, dim1, dim2, ...);

     Make input array A into an array with dimension list DIM1, DIM2, ...
     If TYPE is specified, the result will be of that data type
     _w_i_t_h_o_u_t___c_o_n_v_e_r_s_i_o_n_.

   SEE ALSO:
     make_dimlist, reshape, ipy_shift_dims, ipy_cast_real_as_complex,
     ipy_cast_complex_as_real.
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

func ipy_shift_dims(a, n)
/* DOCUMENT ipy_shift_dims(a, n)
     Returns array A with N new leading dimensions (each new dimension has a
     lenght of 1).

   SEE ALSO: reshape, ipy_reshape.
 */
{
  if (n > 0) {
    old_dims = dimsof(a);
    old_rank = old_dims(1);
    new_rank = n + old_rank;
    new_dims = array(1, new_rank + 1);
    new_dims(1) = new_rank;
    if (old_rank >= 1) {
      last = 1 - old_rank : 0;
      new_dims(last) = old_dims(last);
    }
    reshape, a, structof(a), new_dims;
  }
  return a;
}

local ipy_format, ipy_warn, ipy_error;
/* DOCUMENT ipy_format, fmt, ...;
         or ipy_format(fmt, ...);
         or ipy_warn, fmt, ...;
         or ipy_error, fmt, ...;

     When called as a subroutine, ipy_format write a formatted message to
     standard output; when called as a function, it returns the formatted
     string.  FMT must be a string with "%" directives for each other
     arguments.

     The subroutines ipy_warn and ipy_error are simple wrappers to print
     warnings or raise errors with a formatted message.

   SEE ALSO: write, swrite, error.
 */
func ipy_format(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) /* DOCUMENTED */
{
  if (am_subroutine()) {
    if (is_void(a1)) {
      write, format="%s", a0;
    } else if (is_void(a2)) {
      write, format=a0, a1;
    } else if (is_void(a3)) {
      write, format=a0, a1, a2;
    } else if (is_void(a4)) {
      write, format=a0, a1, a2, a3;
    } else if (is_void(a5)) {
      write, format=a0, a1, a2, a3, a4;
    } else if (is_void(a6)) {
      write, format=a0, a1, a2, a3, a4, a5;
    } else if (is_void(a7)) {
      write, format=a0, a1, a2, a3, a4, a5, a6;
    } else if (is_void(a8)) {
      write, format=a0, a1, a2, a3, a4, a5, a6, a7;
    } else if (is_void(a9)) {
      write, format=a0, a1, a2, a3, a4, a5, a6, a7, a8;
    } else {
      write, format=a0, a1, a2, a3, a4, a5, a6, a7, a8, a9;
    }
  } else if (is_void(a1)) {
    return swrite(format="%s", a0);
  } else if (is_void(a2)) {
    return swrite(format=a0, a1);
  } else if (is_void(a3)) {
    return swrite(format=a0, a1, a2);
  } else if (is_void(a4)) {
    return swrite(format=a0, a1, a2, a3);
  } else if (is_void(a5)) {
    return swrite(format=a0, a1, a2, a3, a4);
  } else if (is_void(a6)) {
    return swrite(format=a0, a1, a2, a3, a4, a5);
  } else if (is_void(a7)) {
    return swrite(format=a0, a1, a2, a3, a4, a5, a6);
  } else if (is_void(a8)) {
    return swrite(format=a0, a1, a2, a3, a4, a5, a6, a7);
  } else if (is_void(a9)) {
    return swrite(format=a0, a1, a2, a3, a4, a5, a6, a7, a8);
  } else {
    return swrite(format=a0, a1, a2, a3, a4, a5, a6, a7, a8, a9);
  }
}
func ipy_warn(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) /* DOCUMENTED */
{
  write, format="WARNING - %s\n", ipy_format(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
}
func ipy_error(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9) /* DOCUMENTED */
{
  error, ipy_format(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9);
}
errs2caller, ipy_error;

/*---------------------------------------------------------------------------*/
/* LINEAR CONJUGATE GRADIENTS */

local ipy_conjgrad_iterations;
func ipy_conjgrad(A, b, x, precond=, maxiter=, gtol=, quiet=)
/* DOCUMENT x = ipy_conjgrad(A, b);
         or x = ipy_conjgrad(A, b, x0);

     The function ipy_conjgrad()  solves the linear system A.x =  b by means of
     conjugate   gradient   method.    The   argument  A   is   the   symmetric
     (semi-)definite left-hand side "matrix" which  is used as a function, i.e.
     A(x) returns A.x, and the argument b is the right-hand side "vector".  The
     optional  argument X0  is  the  initial solution;  if  not specified,  the
     initial solution is an array of zeros.

     Keyword PRECOND  can be  used to  specify a preconditioner  M, that  is an
     operator which approximates  the inverse of A.  If specified,  M is called
     as function (like A).

     Keyword GTOL can be used to specify a convergence threshold.  The value of
     GTOL is  either a  scalar GATOL  or a vector  of two  values [GATOL,GRTOL]
     where GATOL and GRTOL are an  absolute and a relative gradient tolerances.
     If GTOL is a scalar, the GRTOL=0  is assumed; if GTOL is unspecified, then
     GATOL=0  and  grtol=1E-3`  are  the  default  tolerances.   The  conjugate
     gradient  iterations  are  performed  until  the  Euclidean  norm  of  the
     (preconditioned)  residuals M.(A.x  -  b)  is less  or  equal the  largest
     between GATOL and GRTOL times the Euclidean norm of the initial residuals.

     Keyword MAXITER  can be used to  specify the maximum number  of iterations
     (50 by default) to perform until  returning the result.  A warning message
     is printed if  the number of iterations is exceeded  (unless keyword QUIET
     is true).

     Variable ipy_conjgrad_iterations can  be used to figure out  the number of
     conjugate gradient iterations.  For instance:

         local ipy_conjgrad_iterations;
         x = ipy_conjgrad(A, b, quiet=1);
         if (ipy_conjgrad_iterations > 100) error, "too many iterations";


   SEE ALSO: ipy_inner, ipy_combine. */
{
  extern ipy_conjgrad_iterations;
  local p, q, r, z;

  /* Check arguments. */
  use_precond = ! is_void(precond);
  if (is_void(maxiter)) maxiter = 50;
  if (is_void(gtol)) {
    gatol = 0.0;
    grtol = 1e-3;
  } else {
    if (identof(gtol) > Y_DOUBLE){
      error, "bad type for GTOL";
    }
    if (numberof(gtol) > 2) {
      error, "GTOL must have 1 or 2 elements";
    }
    if (min(gtol) < 0.0) {
      error, "GTOL must be nonnegative";
    }
    if (is_scalar(gtol)) {
      gatol = double(gtol(1));
      grtol = 0.0;
    } else {
      gatol = double(gtol(1));
      grtol = double(gtol(2));
    }
  }

  /* Compute initial residuals. */
  if (is_void(x)) {
    eq_nocopy, r, unref(b);
    x = array(double, dimsof(b));
  } else {
    r = unref(b) - A(x);
  }
  if (use_precond) {
    z = precond(r);
  } else {
    eq_nocopy, z, r;
  }
  rho = ipy_inner(z, z);
  epsilon = max(gatol, grtol*sqrt(rho));

  /* Conjugate gradient iterations. */
  beta = 0.0;
  ipy_conjgrad_iterations = 0;
  while (sqrt(rho) > epsilon) {
    /* Check number of iterations. */
    if (ipy_conjgrad_iterations >= maxiter) {
      if (! quiet) warn, swrite(format="too many iterations (%d)", maxiter);
      break;
    }

    /* Next search direction. */
    if (beta == 0.0) {
      eq_nocopy, p, z;
    } else {
      ipy_combine, p, 1.0, z, beta, p;
    }

    /* Make optimal step along search direction. */
    q = A(p);
    gamma = ipy_inner(p, q);
    if (gamma <= 0.0) error, "operator A is not positive definite";
    alpha = rho/gamma;
    ipy_combine, x, 1.0, x, +alpha, p;
    ipy_combine, r, 1.0, r, -alpha, q;
    rho_prev = rho;
    if (use_precond) {
      z = precond(r);
    } else {
      eq_nocopy, z, r;
    }
    rho = ipy_inner(z, z);
    beta = rho/rho_prev;
    ++ipy_conjgrad_iterations;
  }
  return x;
}

/*---------------------------------------------------------------------------*/
/* BENCHMARKING */

local __ipy_benchmark_proc;
func ipy_benchmark(__ipy_benchmark_script,
                   __ipy_benchmark_repeat,
                   __ipy_benchmark_warmup)
/* DOCUMENT t = ipy_benchmark(script);
         or t = ipy_benchmark(script, repeat);
         or t = ipy_benchmark(script, repeat, warmup);
         or ipy_benchmark, script;
         or ipy_benchmark, script, repeat;
         or ipy_benchmark, script, repeat, warmup;

     Measure the  time spent  by executing  Yorick code  in the  SCRIPT string.
     Argument SCRIPT can be a scalar  string or an array of strings (typically,
     one per line of code).  Optional argument REPEAT gives the number of times
     the script is executed; if omitted, its default value is 100.  When called
     as a function,  the returned value is the elapsed  times [CPU,SYS,WALL] in
     seconds  (see `timer`).   Optional argument  `warmup` can  be set  with an
     integer to  warmup your CPU  before timing the  script.  When called  as a
     subroutine, the results are printed to standard output.

     The function `__ipy_benchmark_proc`  is created on the fly  to perform the
     benchmark.   A   consequence  is   that  it  is   not  possible   to  call
     `ipy_benchmark`   from    the   script.    Also   symbols    prefixed   by
     `__ipy_benchmark_` should not be used by  the script (they are reserved to
     implement the `ipy_benchmark` function).


   SEE ALSO: timer, include.
 */
{
  extern __ipy_benchmark_proc;
  if (is_void(__ipy_benchmark_repeat)) __ipy_benchmark_repeat = 100;
  if (! is_string(__ipy_benchmark_script)) error, "SCRIPT must be a string";
  include, grow("func __ipy_benchmark_proc(__ipy_benchmark__counter) {",
                "  while (--__ipy_benchmark__counter >= 0) { ",
                __ipy_benchmark_script,
                "  }",
                "}"), 1;
  if (__ipy_benchmark_warmup) {
    __ipy_benchmark_arr = random(10000); // not too big to fit in cache
    while (__ipy_benchmark_warmup-- > 0) {
      __ipy_benchmark_tmp = sum(__ipy_benchmark_arr);
    }
  }
  __ipy_benchmark_t0 = __ipy_benchmark_t1 = array(double, 3);
  timer, __ipy_benchmark_t0;
  __ipy_benchmark_proc, __ipy_benchmark_repeat;
  timer, __ipy_benchmark_t1;
  __ipy_benchmark_proc = [];
  __ipy_benchmark_t = (__ipy_benchmark_t1 -
                       __ipy_benchmark_t0)/__ipy_benchmark_repeat;
  if (am_subroutine()) {
    write, format="cpu=%gms, system=%gms, wall=%gms (measured for %d iteration%s)\n",
      __ipy_benchmark_t(1)*1e3,
      __ipy_benchmark_t(2)*1e3,
      __ipy_benchmark_t(3)*1e3,
      __ipy_benchmark_repeat,
      (__ipy_benchmark_repeat > 1 ? "s" : "");
  }
  return __ipy_benchmark_t;
}

/*---------------------------------------------------------------------------*/
/* CONSTANTS AND INITIALIZATION */

local IPY_HUGE, IPY_TINY, IPY_EPSILON, IPY_PI;
local IPY_TRUE, IPY_FALSE;
extern ipy_init;
/* DOCUMENT Constants defined by IPY package

     IPY_TRUE    - true value;
     IPY_FALSE   - false value;
     IPY_HUGE    - largest positive real;
     IPY_TINY    - smallest strictly positive real;
     IPY_EPSILON - relative machine precision (smallest strictly positive real
                   such that 1 + EPSILON is numerically different from 1);
     IPY_PI      - value of pi;
 */
ipy_init;
