# Cost Functions

Solving an inverse problem usually amounts to minimizing some cost
function.

For problems of small size, some optimizers only require the
function(differentiable free solvers).  For smooth functions, the (trust
region) Newton method is the most efficient and require the function, its
gradient and the Hessian.

For smooth functions and large size problems, the optimizers require the
function and its gradient.

For non-smooth functions, the proximal operator may be used.

The total cost function is, in genral, a combination of several cost
functions.  For instance:
```
f(x) = fdata(x) + mu*fprior(x)
```
with `fdata(x)` the distance of the model to the data, `fprior(x)` a
regularization term and `mu` a nonnegative weight.


## Interface

Computing the cost function is done by:
```
cost = ipy_compute_cost(alpha, obj, x);
```
with `alpha` a nonnegative weight, `obj` an *object* representing the
objective function (and its parameters) and `x` the variables.  The
returned value is `cost = alpha*f(x)`.

Computing the cost function and its gradient is done by:
```
local grad;
cost = ipy_compute_cost_and_gradient(alpha, obj, x, grad);
```
where `grad` is a variable to store the gradient of `f(x)` times `alpha`.

Applying the proximal operator is done by:
```
xp = ipy_apply_prox(alpha, obj, x, tol);
```
which returns the solution of:
```
min_z { alpha*f(z) + (1/2) ||z - x||^2 }
```
with some tolerance `tol`.


## Implementation

To implement the above behavior, each *class* of cost function should
provide the following *methods* (some methods can be omitted to indicate
that this capabilty is missing):
```
func _prefix_class_f(alpha, ctx, x) {
    return alpha*f(x);
}
func _prefix_class_fg(alpha, ctx, x, &grad) {
    grad = alpha*g(x);
    return alpha*f(x);
}
func _prefix_class_prox(alpha, ctx, x, tol) {
    return arg_min(alpha*f(z) + 0.5*norm2(z - x));
}
```
with `ctx` some contextual data and where the body of the functions are
indicative.


To build a cost function:
```
fn = ipy_new_cost();
```
