/*
 * post-refinement.c
 *
 * Post refinement
 *
 * (c) 2006-2010 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdlib.h>
#include <gsl/gsl_poly.h>
#include <assert.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_linalg.h>

#include "image.h"
#include "post-refinement.h"
#include "peaks.h"
#include "symmetry.h"
#include "geometry.h"


/* Return the gradient of parameter 'k' given the current status of 'image'. */
double gradient(struct image *image, int k,
                       struct cpeak spot, double I_partial)
{
	double ds;
	double nom, den;

	ds = 2.0 * resolution(image->indexed_cell, spot.h, spot.k, spot.l);

	switch ( k ) {

	case REF_SCALE :
		return I_partial;

	case REF_DIV :
		nom = sqrt(2.0) * ds * sin(image->div);
		den = sqrt(1.0 - cos(image->div));
		return nom/den;

	}

	ERROR("No gradient defined for parameter %i\n", k);
	abort();
}


/* Apply the given shift to the 'k'th parameter of 'image'. */
void apply_shift(struct image *image, int k, double shift)
{
	switch ( k ) {

	case REF_SCALE :
		image->osf += shift;
		break;

	case REF_DIV :
		STATUS("Shifting div by %e\n", shift);
		image->div += shift;
		break;

	default :
		ERROR("No shift defined for parameter %i\n", k);
		abort();

	}
}


double mean_partial_dev(struct image *image, struct cpeak *spots, int n,
                        const char *sym, double *i_full, FILE *graph)
{
	int h;
	double delta_I = 0.0;

	for ( h=0; h<n; h++ ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full;
		float I_partial;
		float xc, yc;

		hind = spots[h].h;
		kind = spots[h].k;
		lind = spots[h].l;

		/* Don't attempt to use spots with very small
		 * partialities, since it won't be accurate. */
		if ( spots[h].p < 0.1 ) continue;

		/* Actual measurement of this reflection from this
		 * pattern? */
		/* FIXME: Coordinates aren't whole numbers */
		if ( integrate_peak(image, spots[h].x, spots[h].y,
		                    &xc, &yc, &I_partial, NULL, NULL, 1, 1) ) {
			continue;
		}
		I_partial *= image->osf;

		get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
		I_full = lookup_intensity(i_full, ha, ka, la);
		delta_I += fabs(I_partial - spots[h].p * I_full);

		if ( graph != NULL ) {
			fprintf(graph, "%3i %3i %3i %5.2f (at %5.2f,%5.2f)"
			               " %5.2f %5.2f\n",
			       hind, kind, lind, I_partial/spots[h].p, xc, yc,
			       spots[h].p, I_partial / I_full);
		}

	}

	return delta_I / (double)n;
}


/* Perform one cycle of post refinement on 'image' against 'i_full' */
double pr_iterate(struct image *image, double *i_full, const char *sym,
                  struct cpeak **pspots, int *n)
{
	gsl_matrix *M;
	gsl_vector *v;
	gsl_vector *shifts;
	int h, param;
	struct cpeak *spots = *pspots;

	M = gsl_matrix_calloc(NUM_PARAMS, NUM_PARAMS);
	v = gsl_vector_calloc(NUM_PARAMS);

	/* Construct the equations, one per reflection in this image */
	for ( h=0; h<*n; h++ ) {

		signed int hind, kind, lind;
		signed int ha, ka, la;
		double I_full, delta_I;
		float I_partial;
		float xc, yc;
		int k;

		hind = spots[h].h;
		kind = spots[h].k;
		lind = spots[h].l;

		/* Don't attempt to use spots with very small
		 * partialities, since it won't be accurate. */
		if ( spots[h].p < 0.1 ) continue;

		/* Actual measurement of this reflection from this
		 * pattern? */
		/* FIXME: Coordinates aren't whole numbers */
		if ( integrate_peak(image, spots[h].x, spots[h].y,
		                    &xc, &yc, &I_partial, NULL, NULL, 1, 1) ) {
			continue;
		}
		I_partial *= image->osf;

		get_asymm(hind, kind, lind, &ha, &ka, &la, sym);
		I_full = lookup_intensity(i_full, ha, ka, la);
		delta_I = I_partial - spots[h].p * I_full;

		for ( k=0; k<NUM_PARAMS; k++ ) {

			int g;
			double v_c;

			for ( g=0; g<NUM_PARAMS; g++ ) {

				double M_curr, M_c;

				M_curr = gsl_matrix_get(M, g, k);

				M_c = gradient(image, g, spots[h], I_partial)
				    * gradient(image, k, spots[h], I_partial);
				M_c *= pow(I_full, 2.0);

				gsl_matrix_set(M, g, k, M_curr + M_c);

			}

			v_c = delta_I * I_full * gradient(image, k, spots[h],
			                                  I_partial);
			gsl_vector_set(v, k, v_c);

		}

	}

	shifts = gsl_vector_alloc(NUM_PARAMS);
	gsl_linalg_HH_solve(M, v, shifts);
	for ( param=0; param<NUM_PARAMS; param++ ) {
		double shift = gsl_vector_get(shifts, param);
		apply_shift(image, param, shift);
	}

	gsl_matrix_free(M);
	gsl_vector_free(v);
	gsl_vector_free(shifts);

	free(spots);
	spots = find_intersections(image, image->indexed_cell, n, 0);
	*pspots = spots;
	return mean_partial_dev(image, spots, *n, sym, i_full, NULL);
}
