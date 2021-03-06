% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simevol1.R
\name{simevol}
\alias{simevol}
\title{function for evolutionary simulation}
\usage{
simevol(z, en, fitness = NULL, mutate = NULL, pop_dynamics = NULL,
  set_parms = NULL, edge_die = 1e-06, edge_fit = 1e-12,
  tmax = 1e+08, times = 10^seq(0, 20, by = 1) - 1, out_interval = 10,
  show_interval = 10, amp_invf = 0.1, level_invf = 0.01,
  file_data = "test.dat", fitness_contour = TRUE, yorick_plot = NULL,
  plot_mask = NULL, reset_win = TRUE, trait_names = NULL,
  env_names = NULL, param = param0)
}
\arguments{
\item{z}{: traits}
}
\description{
function for evolutionary simulation.
}
\author{
Hiroshi C. Ito
}
