/*
 * ipy.c --
 *
 * Implements basic "vector" operations in IPY.
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

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include <pstdlib.h>
#include <play.h>
#include <yapi.h>

/* Define some macros to get rid of some GNU extensions when not compiling
   with GCC. */
#if ! (defined(__GNUC__) && __GNUC__ > 1)
#   define __attribute__(x)
#   define __inline__
#   define __FUNCTION__        ""
#   define __PRETTY_FUNCTION__ ""
#endif

#define TRUE  1
#define FALSE 0

PLUG_API void y_error(const char *) __attribute__ ((noreturn));

static int same_dims(const long xdims[], const long ydims[]);

/* Define a Yorick global symbol with an scalar value. */
static void define_int_const(const char *name, int value);
static void define_double_const(const char *name, double value);

/*---------------------------------------------------------------------------*/
/* BUILT-IN FUNCTIONS */

void
Y_ipy_norm_1(int argc)
{
  double sum = 0.0;
  long ntot, i;
  int type;

  if (argc != 1) y_error("expecting exactly one argument");
  type = yarg_typeid(0);
  if (type > Y_DOUBLE) y_error("bad argument type");
  if (type == Y_DOUBLE || type == Y_LONG) {
    const double *v = ygeta_d(0, &ntot, NULL);
    for (i = 0; i < ntot; ++i) {
      sum += fabs(v[i]);
    }
  } else {
    const float *v = ygeta_f(0, &ntot, NULL);
    for (i = 0; i < ntot; ++i) {
      sum += fabsf(v[i]);
    }
  }
  ypush_double(sum);
}

void
Y_ipy_norm_2(int argc)
{
  double sum = 0.0;
  long ntot, i;
  int type;

  if (argc != 1) y_error("expecting exactly one argument");
  type = yarg_typeid(0);
  if (type > Y_DOUBLE) y_error("bad argument type");
  if (type == Y_DOUBLE || type == Y_LONG) {
    const double *v = ygeta_d(0, &ntot, NULL);
    for (i = 0; i < ntot; ++i) {
      sum += v[i]*v[i];
    }
  } else {
    const float *v = ygeta_f(0, &ntot, NULL);
    for (i = 0; i < ntot; ++i) {
      sum += v[i]*v[i];
    }
  }
  ypush_double(sqrt(sum));
}

void
Y_ipy_norm_inf(int argc)
{
  double result;
  long ntot, i;
  int type;

  if (argc != 1) y_error("expecting exactly one argument");
  type = yarg_typeid(0);
  if (type > Y_DOUBLE) y_error("bad argument type");
  if (type == Y_DOUBLE || type == Y_LONG) {
    const double *v = ygeta_d(0, &ntot, NULL);
    double val, maxval = 0;
    for (i = 0; i < ntot; ++i) {
      val = fabs(v[i]);
      if (val > maxval) maxval = val;
    }
    result = maxval;
  } else {
    const float *v = ygeta_f(0, &ntot, NULL);
    float val, maxval = 0;
    for (i = 0; i < ntot; ++i) {
      val = fabsf(v[i]);
      if (val > maxval) maxval = val;
    }
    result = maxval;
  }
  ypush_double(result);
}

void
Y_ipy_inner(int argc)
{
  double sum = 0.0;
  long dims[Y_DIMSIZE];
  long tdims[Y_DIMSIZE];
  void *wptr;
  void *xptr;
  void *yptr;
  long ntot, i;
  int wtype, xtype, ytype, type;

  if (argc != 2 && argc != 3) y_error("expecting 2 or 3 arguments");
  xptr = ygeta_any(1, &ntot, dims, &xtype);
  if (xtype > Y_DOUBLE) y_error("bad data type for X");
  type = ((xtype == Y_DOUBLE || xtype == Y_LONG) ? Y_DOUBLE : Y_FLOAT);

  yptr = ygeta_any(0, NULL, tdims, &ytype);
  if (ytype > Y_DOUBLE) y_error("bad data type for Y");
  if (! same_dims(dims, tdims)) y_error("X and Y have not same dimensions");
  if (ytype == Y_DOUBLE || ytype == Y_LONG) type = Y_DOUBLE;

  if (argc >= 3) {
    wptr = ygeta_any(2, NULL, tdims, &wtype);
    if (wtype > Y_DOUBLE) y_error("bad data type for W");
    if (! same_dims(dims, tdims)) y_error("W and X have not same dimensions");
    if (wtype == Y_DOUBLE || wtype == Y_LONG) type = Y_DOUBLE;
    if (wtype != type) wptr = ygeta_coerce(2, wptr, ntot, dims, wtype, type);
  } else {
    wptr = NULL;
  }
  if (xtype != type) xptr = ygeta_coerce(1, xptr, ntot, dims, xtype, type);
  if (ytype != type) yptr = ygeta_coerce(0, yptr, ntot, dims, ytype, type);

  if (argc == 3) {
    if (type == Y_FLOAT) {
      const float *w = (const float *)wptr;
      const float *x = (const float *)xptr;
      const float *y = (const float *)yptr;
      for (i = 0; i < ntot; ++i) {
        sum += w[i]*x[i]*y[i];
      }
    } else {
      const double *w = (const double *)wptr;
      const double *x = (const double *)xptr;
      const double *y = (const double *)yptr;
      for (i = 0; i < ntot; ++i) {
        sum += w[i]*x[i]*y[i];
      }
    }
  } else {
    if (type == Y_FLOAT) {
      const float *x = (const float *)xptr;
      const float *y = (const float *)yptr;
      for (i = 0; i < ntot; ++i) {
        sum += x[i]*y[i];
      }
    } else {
      const double *x = (const double *)xptr;
      const double *y = (const double *)yptr;
      for (i = 0; i < ntot; ++i) {
        sum += x[i]*y[i];
      }
    }
  }

  ypush_double(sum);
}

void
Y_ipy_combine(int argc)
{
  double a[5];
  long dims[Y_DIMSIZE];
  long tdims[Y_DIMSIZE];
  void *rptr;
  void *v[5];
  long ntot, i, index;
  int rtype, vtype[5], type, iarg[5], j, jp, n;
  int fresh;

  if (argc <  2) y_error("too few arguments");
  if (argc > 11) y_error("too many arguments");
  if ((argc&1) == 0) {
    if (yarg_subroutine()) y_error("missing destination variable");
    index = -1;
    iarg[0] = argc - 2;
  } else {
    index = yget_ref(argc - 1);
    if (index < 0) y_error("destination must be a simple variable");
    iarg[0] = argc - 3;
  }
  a[0] = ygets_d(iarg[0] + 1);
  v[0] = ygeta_any(iarg[0], &ntot, dims, &vtype[0]);
  if (vtype[0] > Y_DOUBLE) y_error("bad operand data type");
  type = ((vtype[0] == Y_DOUBLE || vtype[0] == Y_LONG) ? Y_DOUBLE : Y_FLOAT);
  n = argc/2;
  for (j = 1; j < n; ++j) {
    iarg[j] = iarg[j-1] - 2;
    a[j] = ygets_d(iarg[j] + 1);
    v[j] = ygeta_any(iarg[j], NULL, tdims, &vtype[j]);
    if (vtype[j] > Y_DOUBLE) y_error("bad operand data type");
    if (vtype[j] == Y_DOUBLE || vtype[j] == Y_LONG) type = Y_DOUBLE;
    if (! same_dims(dims, tdims)) {
      y_error("all operands must have same dimensions");
    }
  }
  for (j = 0; j < n; ++j) {
    if (vtype[j] != type) {
      v[j] = ygeta_coerce(iarg[j], v[j], ntot, dims, vtype[j], type);
    }
  }
  rptr = NULL;
  if ((argc&1) != 0) {
    /* Get current value of destination. */
    if (yarg_typeid(argc-1) == type) {
      rptr = ygeta_any(argc-1, NULL, tdims, &rtype);
      if (rtype == type && same_dims(tdims, dims)) {
        index = -1; /* no needs to redefine */
      } else {
        rptr = NULL;
      }
    }
    if (rptr == NULL) {
      /* Discard destination value because it cannot be reused. */
      ypush_nil();
      yarg_swap(argc, 0);
      yarg_drop(1);
    }
  }

  /* Create the destination if needed. */
  fresh = (rptr == NULL);
  if (fresh) {
    if (type == Y_FLOAT) {
      rptr = ypush_f(dims);
    } else {
      rptr = ypush_d(dims);
    }
  }

  /* Only keep operands with a non-zero weight. */
  for (jp = 0, j = 0; j < n; ++j) {
    if (a[j] != 0.0) {
      if (jp != j) {
        a[jp] = a[j];
        v[jp] = v[j];
      }
      ++jp;
    }
  }
  n = jp;

  /* Apply the operation. */
  if (n == 0) {
    if (! fresh) {
      memset(rptr, 0, ntot*(type == Y_FLOAT ? sizeof(float)
                            : sizeof(double)));
    }
  } else if (n == 1) {
    if (type == Y_FLOAT) {
      float *r = (float *)rptr;
      const float *v0 = (const float *)v[0];
      float a0 = a[0];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i];
      }
    } else {
      double *r = (double *)rptr;
      const double *v0 = (const double *)v[0];
      double a0 = a[0];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i];
      }
    }
  } else if (n == 2) {
    if (type == Y_FLOAT) {
      float *r = (float *)rptr;
      const float *v0 = (const float *)v[0];
      const float *v1 = (const float *)v[1];
      float a0 = a[0];
      float a1 = a[1];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i];
      }
    } else {
      double *r = (double *)rptr;
      const double *v0 = (const double *)v[0];
      const double *v1 = (const double *)v[1];
      double a0 = a[0];
      double a1 = a[1];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i];
      }
    }
  } else if (n == 3) {
    if (type == Y_FLOAT) {
      float *r = (float *)rptr;
      const float *v0 = (const float *)v[0];
      const float *v1 = (const float *)v[1];
      const float *v2 = (const float *)v[2];
      float a0 = a[0];
      float a1 = a[1];
      float a2 = a[2];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i];
      }
    } else {
      double *r = (double *)rptr;
      const double *v0 = (const double *)v[0];
      const double *v1 = (const double *)v[1];
      const double *v2 = (const double *)v[2];
      double a0 = a[0];
      double a1 = a[1];
      double a2 = a[2];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i];
      }
    }
  } else if (n == 4) {
    if (type == Y_FLOAT) {
      float *r = (float *)rptr;
      const float *v0 = (const float *)v[0];
      const float *v1 = (const float *)v[1];
      const float *v2 = (const float *)v[2];
      const float *v3 = (const float *)v[3];
      float a0 = a[0];
      float a1 = a[1];
      float a2 = a[2];
      float a3 = a[3];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i] + a3*v3[i];
      }
    } else {
      double *r = (double *)rptr;
      const double *v0 = (const double *)v[0];
      const double *v1 = (const double *)v[1];
      const double *v2 = (const double *)v[2];
      const double *v3 = (const double *)v[3];
      double a0 = a[0];
      double a1 = a[1];
      double a2 = a[2];
      double a3 = a[3];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i] + a3*v3[i];
      }
    }
  } else if (n == 5) {
    if (type == Y_FLOAT) {
      float *r = (float *)rptr;
      const float *v0 = (const float *)v[0];
      const float *v1 = (const float *)v[1];
      const float *v2 = (const float *)v[2];
      const float *v3 = (const float *)v[3];
      const float *v4 = (const float *)v[4];
      float a0 = a[0];
      float a1 = a[1];
      float a2 = a[2];
      float a3 = a[3];
      float a4 = a[4];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i] + a3*v3[i] + a4*v4[i];
      }
    } else {
      double *r = (double *)rptr;
      const double *v0 = (const double *)v[0];
      const double *v1 = (const double *)v[1];
      const double *v2 = (const double *)v[2];
      const double *v3 = (const double *)v[3];
      const double *v4 = (const double *)v[4];
      double a0 = a[0];
      double a1 = a[1];
      double a2 = a[2];
      double a3 = a[3];
      double a4 = a[4];
      for (i = 0; i < ntot; ++i) {
        r[i] = a0*v0[i] + a1*v1[i] + a2*v2[i] + a3*v3[i] + a4*v4[i];
      }
    }
  } else {
    y_error("there must be a bug!");
  }

  if (index >= 0) {
    /* Define destination variable (result is already on top of the stack). */
    yput_global(index, 0);
  } else if ((argc&1) != 0) {
    /* Drop stack elements to left result on top of the stack. */
    yarg_drop(argc - 1);
  }
}

#ifndef M_PI
# define M_PI 3.14159265358979323846 /* pi */
#endif

void
Y_ipy_init(int argc)
{
  /* Define constants. */
  define_double_const("IPY_HUGE", DBL_MAX);
  define_double_const("IPY_TINY", DBL_MIN);
  define_double_const("IPY_EPSILON", DBL_EPSILON);
  define_double_const("IPY_PI", M_PI);
  define_int_const("IPY_TRUE", 1);
  define_int_const("IPY_FALSE", 0);

  ypush_nil();
}

/*---------------------------------------------------------------------------*/
/* UTILITIES */

static int
same_dims(const long xdims[], const long ydims[])
{
  int i, xrank, yrank;
  if (xdims == ydims) return TRUE;
  xrank = (xdims == NULL ? 0 : xdims[0]);
  yrank = (ydims == NULL ? 0 : ydims[0]);
  if (xrank != yrank) return FALSE;
  for (i = 1; i <= xrank; ++i) {
    if (xdims[i] != ydims[i]) return FALSE;
  }
  return TRUE;
}

static void
define_int_const(const char *name, int value)
{
  ypush_int(value);
  yput_global(yget_global(name, 0), 0);
  yarg_drop(1);
}

static void
define_double_const(const char *name, double value)
{
  ypush_double(value);
  yput_global(yget_global(name, 0), 0);
  yarg_drop(1);
}
