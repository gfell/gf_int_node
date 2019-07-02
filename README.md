#  Interpolation for ordered arrays in **Fortran** 

This repository contains subroutines to find the index of the right
interpolating node for each data point in an ordered array. It speeds
up computation by exploiting the ordering of both the array of
interpolating nodes and of the array of data points.


## Contains

- *mod_int_node* for standard interpolation 
- *mod_integm_node* for interpolation in economic models with
  occasionally binding borrowing constraints when the first element of
  the array of interpolating nodes coincides with the lower bound on
  assets (as, e.g., in Carroll [2006.  "The method of endogenous
  gridpoints for solving dynamic stochastic optimization problems,"
  _Economics Letters_ 91]). The code returns index 0 to flag data
  points for which the borrowing constraint is binding and the
  solution can be computed analytically.
