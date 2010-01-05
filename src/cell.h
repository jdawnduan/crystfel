/*
 * cell.h
 *
 * Unit Cell Calculations
 *
 * (c) 2007-2009 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */

#ifndef CELL_H
#define CELL_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {

	/* Crystallographic representation */
	double a;	/* m */
	double b;	/* m */
	double c;	/* m */
	double alpha;	/* Radians */
	double beta;	/* Radians */
	double gamma;	/* Radians */

	/* Cartesian representation */
	double ax;	double bx;	double cx;
	double ay;	double by;	double cy;
	double az;	double bz;	double cz;

	/* Cartesian representation of reciprocal axes */
	double axs;	double bxs;	double cxs;
	double ays;	double bys;	double cys;
	double azs;	double bzs;	double czs;

} UnitCell;

extern UnitCell *cell_new(void);

/* Lengths in m, angles in radians */
extern UnitCell *cell_new_from_parameters(double a, double b, double c,
				double alpha, double beta, double gamma);

extern void cell_set_cartesian(UnitCell *cell,
			double ax, double ay, double az,
			double bx, double by, double bz,
			double cx, double cy, double cz);

extern void cell_set_parameters(UnitCell *cell, double a, double b, double c,
				double alpha, double beta, double gamma);

extern void cell_get_parameters(UnitCell *cell, double *a, double *b, double *c,
                         double *alpha, double *beta, double *gamma);

extern void cell_get_cartesian(UnitCell *cell,
                               double *ax, double *ay, double *az,
                               double *bx, double *by, double *bz,
                               double *cx, double *cy, double *cz);

extern void cell_set_cartesian_x(UnitCell *cell,
                                double ax, double bx, double cx);
extern void cell_set_cartesian_y(UnitCell *cell,
                                 double ay, double by, double cy);
extern void cell_set_cartesian_z(UnitCell *cell,
                                 double az, double bz, double cz);

extern void cell_get_reciprocal(UnitCell *cell,
                               double *asx, double *asy, double *asz,
                               double *bsx, double *bsy, double *bsz,
                               double *csx, double *csy, double *csz);

extern double resolution(UnitCell *cell,
                         signed int h, signed int k, signed int l);

#endif	/* CELL_H */
