/*
 * predict-refine.h
 *
 * Prediction refinement
 *
 * Copyright © 2012-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2010-2019 Thomas White <taw@physics.org>
 *
 * This file is part of CrystFEL.
 *
 * CrystFEL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CrystFEL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CrystFEL.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef PREDICT_REFINE_H
#define PREDICT_REFINE_H

struct reflpeak;

/** Enumeration of parameters which may want to be refined */
enum gparam {
	GPARAM_ASX,
	GPARAM_ASY,
	GPARAM_ASZ,
	GPARAM_BSX,
	GPARAM_BSY,
	GPARAM_BSZ,
	GPARAM_CSX,
	GPARAM_CSY,
	GPARAM_CSZ,
	GPARAM_DET_TX,
	GPARAM_DET_TY,
	GPARAM_DET_TZ,
	GPARAM_DET_RX,  /* Detector panel (group) rotation about +x */
	GPARAM_DET_RY,  /* Detector panel (group) rotation about +y */
	GPARAM_DET_RZ,  /* Detector panel (group) rotation about +z */
};


#include "crystal.h"
#include "crystfel-mille.h"

struct reflpeak {
	Reflection *refl;
	struct imagefeature *peak;
	double Ih;   /* normalised */
	struct detgeom_panel *panel;  /* panel the reflection appears on
                                       * (we assume this never changes) */
};

/**
 * \file predict-refine.h
 * Prediction refinement: refinement of indexing solutions before integration.
 */

extern int refine_prediction(struct image *image, Crystal *cr, Mille *mille);

extern int refine_radius(Crystal *cr, struct image *image);

extern double r_dev(struct reflpeak *rp);

extern double fs_dev(struct reflpeak *rp, struct detgeom *det);

extern double ss_dev(struct reflpeak *rp, struct detgeom *det);

extern double r_gradient(UnitCell *cell, int k, Reflection *refl,
                         struct image *image);

extern double fs_gradient(int param, Reflection *refl, UnitCell *cell,
                          struct detgeom_panel *p);

extern double ss_gradient(int param, Reflection *refl, UnitCell *cell,
                          struct detgeom_panel *p);

#endif	/* PREDICT_REFINE_H */
