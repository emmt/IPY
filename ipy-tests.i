/*
 * ipy-tests.i --
 *
 * Suite of tests and benchmarks for IPY package.
 *
 *-----------------------------------------------------------------------------
 *
 * Copyright (C) 2015: Éric Thiébaut <eric.thiebaut@univ-lyon1.fr>
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

func ipy_test_compare(fmt, val, ref)
{
  c = avg(abs(val - ref));
  d = avg(abs(ref));
  r = (c > 0.0 && d > 0.0 ? c/d : c);
  write, format=fmt, (is_scalar(val) ? val : avg(abs(val))), r;
}

func ipy_bench1
{
  dims = [6, 3,4,5,6,7,8];
  w = random_n(dims);
  x = random_n(dims);
  y = random_n(dims);
  n = numberof(x);

  testing = "\n*** Testing %s ***\n";

  write, format=testing, "L1 norm";
  ipy_test_compare, " ipy_norm_1(x) = %g (relative error: %e)\n",
    ipy_norm_1(x), sum(abs(x));
  expr = "sum(abs(x))";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_norm_1(x)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "L2 (Euclidean) norm";
  ipy_test_compare, " ipy_norm_2(x) = %g (relative error: %e)\n",
    ipy_norm_2(x), sqrt(sum(x*x));
  expr = "sqrt(sum(x*x))";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_norm_2(x)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "infinite norm";
  ipy_test_compare, " ipy_norm_inf(x) = %g (relative error: %e)\n",
    ipy_norm_inf(x), max(abs(x));
  expr = "max(abs(x))";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_norm_inf(x)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "inner product";
  ipy_test_compare, " <x,y> = %g (relative error: %e)\n",
    ipy_inner(x,y), sum(x*y);
  expr = "sum(x*y)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "(x(*))(+)*(y(*))(+)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_inner(x,y)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "weighted inner product";
  ipy_test_compare, " <w,x,y> = %g (relative error: %e)\n",
    ipy_inner(w,x,y), sum(w*x*y);
  expr = "sum(w*x*y)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_inner(w,x,y)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;


  v1 = random_n(dims);
  v2 = random_n(dims);
  v3 = random_n(dims);
  a1 = 1;
  a2 = 2;
  a3 = 3;

  write, format=testing, "linear combinations (1 operand)";
  ipy_test_compare, " a1*v1 = %g (relative error: %e)\n",
    ipy_combine(a1,v1), a1*v1;
  expr = "a1*v1";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(a1,v1)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "linear combinations (2 operands)";
  ipy_test_compare, " a1*v1+a2*v2 = %g (relative error: %e)\n",
    ipy_combine(a1,v1,a2,v2), a1*v1+a2*v2;
  expr = "a1*v1+a2*v2";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(a1,v1,a2,v2)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "linear combinations (3 operands)";
  ipy_test_compare, " a1*v1+a2*v2+a3*v3 = %g (relative error: %e)\n",
    ipy_combine(a1,v1,a2,v2,a3,v3), a1*v1+a2*v2+a3*v3;
  expr = "a1*v1+a2*v2+a3*v3";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(a1,v1,a2,v2,a3,v3)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  dst = array(double, dims);
  dst_addr = strtrim(swrite(&dst));
  write, format=testing, "linear combinations (1 operand and destination pre-allocated)";
  ipy_test_compare, " a1*v1 = %g (relative error: %e)\n",
    ipy_combine(dst, a1,v1), a1*v1;
  expr = "a1*v1";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(dst,a1,v1)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "linear combinations (2 operands and destination pre-allocated)";
  ipy_test_compare, " a1*v1+a2*v2 = %g (relative error: %e)\n",
    ipy_combine(dst,a1,v1,a2,v2), a1*v1+a2*v2;
  expr = "a1*v1+a2*v2";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(dst,a1,v1,a2,v2)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  write, format=testing, "linear combinations (3 operands and destination pre-allocated)";
  ipy_test_compare, " a1*v1+a2*v2+a3*v3 = %g (relative error: %e)\n",
    ipy_combine(dst,a1,v1,a2,v2,a3,v3), a1*v1+a2*v2+a3*v3;
  expr = "a1*v1+a2*v2+a3*v3";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;
  expr = "ipy_combine(dst,a1,v1,a2,v2,a3,v3)";
  write, format="Benchmark %s with %d elements: ", expr, n;
  ipy_benchmark, "t = " + expr, 10000, 1000;

  if (dst_addr != strtrim(swrite(&dst))) error, "address of destination has changed";
}
