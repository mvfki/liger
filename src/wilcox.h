/*
  Mathlib : A C Library of Special Functions
  Copyright (C) 1999-2014  The R Core Team

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or (at
  your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, a copy is available at
  https://www.R-project.org/Licenses/

  SYNOPSIS

    #include <Rmath.h>
    double dwilcox(double x, double m, double n, int give_log)
    double pwilcox(double x, double m, double n, int lower_tail, int log_p)
    double qwilcox(double x, double m, double n, int lower_tail, int log_p);
    double rwilcox(double m, double n)

  DESCRIPTION

    dwilcox	The density of the Wilcoxon distribution.
    pwilcox	The distribution function of the Wilcoxon distribution.
    qwilcox	The quantile function of the Wilcoxon distribution.
    rwilcox	Random variates from the Wilcoxon distribution.

 */

/*
   Note: the checks here for R_CheckInterrupt also do stack checking.

   calloc/free are remapped for use in R, so allocation checks are done there.
   freeing is completed by an on.exit action in the R wrappers.

   The Wilcoxon distribution is calculated using work from Andreas Loeffler
   https://upload.wikimedia.org/wikipedia/commons/f/f5/LoefflerWilcoxonMannWhitneyTest.pdf
   https://upload.wikimedia.org/wikipedia/de/1/19/MannWhitney_151102.pdf
*/
#ifndef NEW_WILCOX_LIB
#define NEW_WILCOX_LIB


/* This counts the number of choices with statistic = k */

extern double Rf_dwilcox(double, double, double, int);

/* args have the same meaning as R function pwilcox */
extern double Rf_pwilcox(double q, double m, double n, int lower_tail, int log_p);


/* x is 'p' in R function qwilcox */

extern double Rf_qwilcox(double x, double m, double n, int lower_tail, int log_p);

extern double Rf_rwilcox(double m, double n);
#endif