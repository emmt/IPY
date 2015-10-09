/* Routines provided by ipy.i */
autoload, "ipy.i", ipy_norm_1, ipy_norm_2, ipy_norm_inf;
autoload, "ipy.i", ipy_inner, ipy_combine, ipy_scale;
autoload, "ipy.i", ipy_new_weighting_operator;
autoload, "ipy.i", ipy_new_diagonal_operator, ipy_is_diagonal_operator;
autoload, "ipy.i", ipy_get_diagonal_of_diagonal_operator;
autoload, "ipy.i", ipy_identity;
autoload, "ipy.i", ipy_mirror;
autoload, "ipy.i", ipy_split_index;
autoload, "ipy.i", ipy_recenter_psf;
autoload, "ipy.i", ipy_centroid;
autoload, "ipy.i", ipy_include;
autoload, "ipy.i", ipy_rms;
autoload, "ipy.i", ipy_reshape, ipy_cast_complex_as_real, ipy_cast_real_as_complex;
autoload, "ipy.i", ipy_shift_dims;
autoload, "ipy.i", ipy_format, ipy_warn, ipy_error;
autoload, "ipy.i", ipy_conjgrad;
autoload, "ipy.i", ipy_benchmark;

/* Routines provided by ipy-cost.i */
autoload, "ipy-cost.i", ipy_compute_cost;
autoload, "ipy-cost.i", ipy_compute_cost_and_gradient;
autoload, "ipy-cost.i", ipy_apply_prox;
autoload, "ipy-cost.i", ipy_new_cost;
autoload, "ipy-cost.i", ipy_implements_f;
autoload, "ipy-cost.i", ipy_implements_fg;
autoload, "ipy-cost.i", ipy_implements_prox;
autoload, "ipy-cost.i", ipy_new_quadratic_cost;
