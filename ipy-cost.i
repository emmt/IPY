/*
 * ipy-cost.i -
 *
 * Implement cost functions for iterative optimization and inverse problems.
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

func ipy_compute_cost(alpha, obj, x)
/* DOCUMENT cost = ipy_compute_cost(alpha, f, x);

      Compute cost function.  Return alpha*f(x) where, alpha is a nonnegative
      weight, f is a cost function object and x are the variables.

    SEE ALSO: ipy_compute_cost_and_gradient, ipy_apply_prox, ipy_new_cost.
*/
{
  if (alpha <= 0.0) {
    if (alpha < 0.0) error, "invalid multiplier";
    return 0.0;
  }
  return obj.f(alpha, obj.ctx, x);
}

func ipy_compute_cost_and_gradient(alpha, obj, x, &grad)
/* DOCUMENT cost = ipy_compute_cost_and_gradient(alpha, f, x, grad);

     Compute cost function  and its gradient.  Return  alpha*f(x) where, alpha
     is a nonnegative weight, f is a cost function object, x are the variables
     and grad is a variable to  store the gradient of alpha*f(x) (with respect
     to x).

    SEE ALSO: ipy_compute_cost, ipy_apply_prox, ipy_new_cost.
*/
{
  if (alpha <= 0.0) {
    if (alpha < 0.0) error, "invalid multiplier";
    grad = array(structof(x), dimsof(x));
    return 0.0;
  }
  return obj.fg(alpha, obj.ctx, x, grad);
}

func ipy_apply_prox(alpha, obj, x, xpinit, tol=, maxiter=)
/* DOCUMENT xp = ipy_apply_prox(alpha, f, x);
         or xp = ipy_apply_prox(alpha, f, x, xpinit);

      Apply proximal operator.  The result is the solution of:

         min_xp { alpha*f(xp) + (1/2) ||xp - x||^2 }

      where alpha is a  nonnegative weight, f is a cost  function object and x
      are, possible  non-feasible, variables.

      The optimization problem may be too difficult to be exactly solved, then
      some tolerance for solving the  problem must be provided.  The tolerance
      can be specified by the keyword TOL.

      For proximal operators which iteratively seek for the solution, optional
      argument  XPINIT provides  the initial  solution to  start with  and the
      following keywords can be used:

          maxiter = Maximum number of iterations.
          tol = Tolerance for convergence.

    SEE ALSO: ipy_compute_cost, ipy_compute_cost_and_gradient, ipy_new_cost.
*/
{
  if (alpha <= 0.0) {
    if (alpha < 0.0) error, "invalid multiplier";
    return x;
  }
  return obj.prox(alpha, obj.ctx, x, xpinit, maxiter=maxiter, tol=tol);
}

func ipy_new_cost(name, ctx=, f=, fg=, prox=)
/* DOCUMENT obj = ipy_new_cost(name, ...);

     Create a  cost function object.   Argument NAME is  the name of  the cost
     function, it may be used  for informative purposes.  The other components
     of the object are specified via the following keywords:

         ctx  = Any contextual data needed by the methods of the objects.
         f    = Method to compute the cost.
         fg   = Method to compute the cost and its gradient.
         prox = Method to apply the proximal operator.

     Note that a "method" is specified as either an object callable as a
     function or as the name of a function.  Not all methods need to be
     specified (the corresponding capability will simply be missing).  The
     definitions and behaviors of the different methods are:

         Definition                         Description
         --------------------------------------------------------------------
         func f(alpha, ctx, x)              return alpha*f(x)
         func fg(alpha, ctx, x, &grad)      return alpha*f(x) and store
                                            alpha*g(x) in grad
         func prox(alpha, ctx, x, xpinit,   return the solution of:
                   tol=, maxiter=)          min_xp { alpha*f(xp)
                                                     + (1/2)||xp -x||^2 }
         --------------------------------------------------------------------

     where  f(x) is  the  cost function,  g(x)  is its  gradient,  alpha is  a
     nonnegative  multiplier, ctx  is the  value of  the eponymous  keyword, x
     represents the  variables, grad is a  variable to store the  gradient and
     tol   is  some   tolerance.   See   documentation  of   ipy_compute_cost,
     ipy_compute_cost_and_gradient, and ipy_apply_prox for more details.

     In  principle,  only  cost  function   builders  call  this  function  to
     instanciate a cost fucntion object of a specific class.

     Note  that all  the methods  involving a  multiplier ALPHA  will only  be
     called with ALPHA > 0 (strictly positive).


   SEE ALSO: ipy_compute_cost, ipy_compute_cost_and_gradient, ipy_apply_prox.
 */
{
  flags = ((is_void(f)    ? 0 : 1) |
           (is_void(fg)   ? 0 : 2) |
           (is_void(prox) ? 0 : 4));
  obj = h_new(name = name, ctx = ctx, flags = flags,
              f = _ipy_get_method(f, "_ipy_bogus_f"),
              fg = _ipy_get_method(fg, "_ipy_bogus_fg"),
              prox = _ipy_get_method(prox, "_ipy_bogus_prox"));
}

func _ipy_get_method(f, def)
{
  if (is_void(f)) {
    eq_nocopy, f, def;
  }
  id = identof(f);
  if (id == Y_STRING && is_scalar(f)) {
    return symlink_to_name(f);
  }
  if (id == Y_OPAQUE || id == Y_FUNCTION || id == Y_BUILTIN) {
    return f;
  }
  error, "illegal method argument";
}
errs2caller, _ipy_get_method;

func _ipy_bogus_f(alpha, ctx, x)
{
  extern obj;
  error, ("`" + obj.name + "` cost does not implement `compute_cost`");
}

func _ipy_bogus_fg(alpha, c, x, &g)
{
  extern obj;
  error, ("`" + obj.name +
          "` cost does not implement `compute_cost_and_gradient`");
}

func _ipy_bogus_prox(alpha, c, x)
{
  extern obj;
  error, ("`" + obj.name + "` cost does not implement `apply_prox`");
}

errs2caller, _ipy_bogus_f, _ipy_bogus_fg, _ipy_bogus_prox;

func ipy_implements_f(obj)    { return ((obj.flags & 1) == 1); }
func ipy_implements_fg(obj)   { return ((obj.flags & 2) == 2); }
func ipy_implements_prox(obj) { return ((obj.flags & 4) == 4); }
/* DOCUMENT ipy_implements_f(obj);
         or ipy_implements_fg(obj);
         or ipy_implements_prox(obj);

     Check whether cost function object OBJ implements a specific method.

   SEE ALSO: ipy_new_cost.
 */

/*---------------------------------------------------------------------------*/
/* QUADRATIC PENALTY */

func ipy_new_quadratic_cost(nil, A=, b=, W=)
/* DOCUMENT obj = ipy_new_quadratic_cost(...);

     Create a cost function object implementing quadratic penalty.
     Using matrix notation, the quadratic penaly writes:

         f(x) = (1/2) (A.x - b)'.W.(A.x - b)

     where x are the variables, A is a linear operator, b is an array and W is
     a  linear  weighting  operator  (W must  be  definite  or  seimi-definite
     positive).  All components  of the quadratic penalty are  be specified by
     keywords:

         A = A function (or an object callable as a function) which implements
             the linear  operator A: A(x)  yields A.x and A(x,1)  yields A'.x.
             If not specified, A is assumed to be the identity.

         b = An array.  If not specified, b = 0 is assumed.

         W = A function (or an object callable as a function) which implements
             the  weighting  operator   W:  W(x)  yields  W.x;   or  an  array
             (conformable with the variables x) of nonegative weights.  If not
             specified, A is assumed to be the identity.

     Depending on which components are specified, the methods implementing the
     behavior of the cost function are optimized.


   SEE ALSO: ipy_new_cost.
 */
{
  if (! is_void(nil)) error, "unexpected non-nil argument";
  if (is_array(W)) W = ipy_new_weighting_operator(W);
  case = ((is_void(A) ? 0 : 1) |
          (noneof(b)  ? 0 : 2) |
          (is_void(W) ? 0 : 4));
  prefix = swrite(format="_ipy_quadratic%d_", case);
  return ipy_new_cost("quadratic",
                      ctx = h_new(A=A, b=b, W=W),
                      f = prefix+"f",
                      fg = prefix+"fg",
                      prox = prefix+"prox");
}

/* Quadratic0: f(x) = (1/2) ||x||^2 */

func _ipy_quadratic0_f(alpha, c, x)
{
  return (alpha/2.0)*ipy_inner(x,x);
}

func _ipy_quadratic0_fg(alpha, c, x, &g)
{
  g = ipy_scale(alpha, x);
  return (alpha/2.0)*ipy_inner(x,x);
}

func _ipy_quadratic0_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  return ipy_scale(1.0/(1.0 + alpha), x);
}

/* Quadratic1: f(x) = (1/2) ||A.x||^2 */

func _ipy_quadratic1_f(alpha, c, x)
{
  Ax = c.A(x);
  return (alpha/2.0)*ipy_inner(Ax,Ax);
}

func _ipy_quadratic1_fg(alpha, c, x, &g)
{
  Ax = c.A(x);
  g = ipy_scale(alpha, c.A(Ax,1));
  return (alpha/2.0)*ipy_inner(Ax,Ax);
}

func _ipy_quadratic1_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.A)) {
    local a;
    eq_nocopy, a, c.A.data;
    d = 1.0 + alpha*a*a;
    return (numberof(d) < numberof(x) ? (1.0/d)*x : x/d);
  }
  return ipy_conjgrad(h_functor("_ipy_quadratic_I_alphaAtA",
                                alpha=alpha, A=c.A),
                      x,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}

func _ipy_quadratic_I_alphaAtA(c, x)
{
  return ipy_combine(1.0, x, c.alpha, c.A(c.A(x),1));
}

func _ipy_quadratic_I_alphaW(c, x)
{
  return ipy_combine(1.0, x, c.alpha, c.W(x));
}

func _ipy_quadratic_I_alphaAtWA(c, x)
{
  return ipy_combine(1.0, x, c.alpha, c.A(c.W(c.A(x)),1));
}

/* Quadratic2: f(x) = (1/2) ||x - b||^2 */

func _ipy_quadratic2_f(alpha, c, x)
{
  r = x - c.b;
  return (alpha/2.0)*ipy_inner(r,r);
}

func _ipy_quadratic2_fg(alpha, c, x, &g)
{
  r = x - c.b;
  g = ipy_scale(alpha, r);
  return (alpha/2.0)*ipy_inner(r,r);
}

func _ipy_quadratic2_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  return ipy_combine(1.0/(1.0 + alpha), x, alpha/(1.0 + alpha), c.b);
}

/* Quadratic3: f(x) = (1/2) ||A.x - b||^2 */

func _ipy_quadratic3_f(alpha, c, x)
{
  r = c.A(x) - c.b;
  return (alpha/2.0)*ipy_inner(r,r);
}

func _ipy_quadratic3_fg(alpha, c, x, &g)
{
  r = c.A(x) - c.b;
  g = ipy_scale(alpha, c.A(r,1));
  return (alpha/2.0)*ipy_inner(r,r);
}

func _ipy_quadratic3_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.A)) {
    local a;
    eq_nocopy, a, c.A.data;
    d = 1.0 + alpha*a*a;
    n = x + alpha*a*c.b;
    return (numberof(d) < numberof(n) ? (1.0/d)*n : n/d);
  }
  if (is_void(c.Atb)) h_set, c, Atb=c.A(c.b,1);
  return ipy_conjgrad(h_functor("_ipy_quadratic_I_alphaAtA",
                                alpha=alpha, A=c.A),
                      x + alpha*c.Atb,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}

/* Quadratic4: f(x) = (1/2) x'.W.x */

func _ipy_quadratic4_f(alpha, c, x)
{
  return (alpha/2.0)*ipy_inner(x, c.W(x));
}

func _ipy_quadratic4_fg(alpha, c, x, &g)
{
  Wx = c.W(x);
  f = (alpha/2.0)*ipy_inner(x,Wx);
  g = ipy_scale(alpha, unref(Wx));
  return f;
}

func _ipy_quadratic4_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.W)) {
    d = 1.0 + alpha*c.W.data;
    return (numberof(d) < numberof(x) ? (1.0/d)*x : x/d);
  }
  return ipy_conjgrad(h_functor("_ipy_quadratic_I_alphaW",
                                alpha=alpha, W=c.W),
                      x,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}

/* Quadratic5: f(x) = (1/2) x'.A'.W.A.x) */

func _ipy_quadratic5_f(alpha, c, x)
{
  Ax = c.A(x);
  return (alpha/2.0)*ipy_inner(Ax, c.W(Ax));
}

func _ipy_quadratic5_fg(alpha, c, x, &g)
{
  Ax = c.A(x);
  WAx = c.W(Ax);
  g = ipy_scale(alpha, c.A(WAx,1));
  return (alpha/2.0)*ipy_inner(Ax, WAx);
}

func _ipy_quadratic5_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.W) && ipy_is_diagonal_operator(c.A)) {
    local a;
    eq_nocopy, a, c.A.data;
    d = 1.0 + alpha*a*c.W.data*a;
    return (numberof(d) < numberof(x) ? (1.0/d)*x : x/d);
  }
  return ipy_conjgrad(h_functor("_ipy_quadratic_I_alphaAtWA",
                                alpha=alpha, A=c.A, W=c.W),
                      x,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}

/* Quadratic6: f(x) = (1/2) (x - b)'.W.(x - b) */

func _ipy_quadratic6_f(alpha, c, x)
{
  r = unref(x) - c.b;
  return (alpha/2.0)*ipy_inner(r, c.W(r));
}

func _ipy_quadratic6_fg(alpha, c, x, &g)
{
  r = unref(x) - c.b;
  Wr = c.W(r);
  f = (alpha/2.0)*ipy_inner(r,Wr);
  g = ipy_scale(alpha, unref(Wr));
  return f;
}

func _ipy_quadratic6_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.W)) {
    local w;
    eq_nocopy, w, c.W.data;
    d = 1.0 + alpha*w;
    n = x + alpha*w*c.b;
    return (numberof(d) < numberof(n) ? (1.0/d)*n : n/d);
  }
  if (is_void(c.Wb)) h_set, c, Wb=c.W(c.b);
  return ipy_conjgrad(h_functor("_ipy_quadratic2_prox_LHS",
                                alpha=alpha, W=c.W),
                      x + alpha*c.Wb,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}

/* Quadratic7: f(x) = (1/2) (A.x -b)'.W.(A.x -b) */

func _ipy_quadratic7_f(alpha, c, x)
{
  r = c.A(x) - c.b;
  return (alpha/2.0)*ipy_inner(r, c.W(r));
}

func _ipy_quadratic7_fg(alpha, c, x, &g)
{
  r = c.A(x) - c.b;
  Wr = c.W(r);
  g = ipy_scale(alpha, c.A(Wr,1));
  return (alpha/2.0)*ipy_inner(r,Wr);
}

func _ipy_quadratic7_prox(alpha, c, x, xpinit, tol=, maxiter=)
{
  if (ipy_is_diagonal_operator(c.W) && ipy_is_diagonal_operator(c.A)) {
    local a, w;
    eq_nocopy, a, c.A.data;
    eq_nocopy, w, c.W.data;
    d = 1.0 + alpha*a*w*a;
    n = x + alpha*a*w*c.b;
    return (numberof(d) < numberof(n) ? (1.0/d)*n : n/d);
  }
  if (is_void(c.AtWb)) h_set, c, AtWb=c.A(c.W(c.b),1);
  return ipy_conjgrad(h_functor("_ipy_quadratic_I_alphaAtWA",
                                alpha=alpha, A=c.A, W=c.W),
                      x + alpha*c.AtWb,
                      (is_void(xpinit) ? x : xpinit),
                      maxiter=maxiter, gtol=[0.0,tol]);
}
