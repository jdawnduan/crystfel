/*
 * crystfel-mille.c
 *
 * Interface to Millepede geometry refinement
 *
 * Copyright © 2023 Deutsches Elektronen-Synchrotron DESY,
 *                  a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Thomas White <taw@physics.org>
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

#include <libcrystfel-config.h>

#include <stdlib.h>
#include <assert.h>

#include "image.h"
#include "geometry.h"
#include "cell-utils.h"
#include "predict-refine.h"
#include "profile.h"

#include <mille_c_wrap.h>


static const enum gparam rv[] =
{
	GPARAM_ASX,
	GPARAM_ASY,
	GPARAM_ASZ,
	GPARAM_BSX,
	GPARAM_BSY,
	GPARAM_BSZ,
	GPARAM_CSX,
	GPARAM_CSY,
	GPARAM_CSZ,
};


int mille_label(int hierarchy_level, int member_index, enum gparam param)
{
	int label;

	assert(member_index < 1000);

	label = 100000*hierarchy_level + 100*member_index;
	switch ( param ) {
		case GPARAM_DET_TX : return label+1;  /* x-shift */
		case GPARAM_DET_TY : return label+2;  /* y-shift */
		case GPARAM_DET_TZ : return label+3;  /* z-shift */
		case GPARAM_DET_RX : return label+4;  /* Rotation around x */
		case GPARAM_DET_RY : return label+5;  /* Rotation around y */
		case GPARAM_DET_RZ : return label+6;  /* Rotation around z */
		default : abort();
	}
}


void write_mille(Mille *mille, int n, UnitCell *cell,
                 struct reflpeak *rps, struct image *image)
{
	float local_gradients[9];
	float global_gradients[64];
	int labels[6];
	int i;
	double asx, asy, asz, bsx, bsy, bsz, csx, csy, csz;

	cell_get_reciprocal(cell, &asx, &asy, &asz,
	                          &bsx, &bsy, &bsz,
	                          &csx, &csy, &csz);

	for ( i=0; i<n; i++ ) {

		signed int h, k, l;
		double xl, yl, zl, kpred;
		int j;
		const struct detgeom_panel_group *group;

		get_indices(rps[i].refl, &h, &k, &l);
		kpred = get_kpred(rps[i].refl);
		xl = h*asx + k*bsx + l*csx;
		yl = h*asy + k*bsy + l*csy;
		zl = h*asz + k*bsz + l*csz;

		/* fs terms: local */
		for ( j=0; j<9; j++ ) {
			local_gradients[j] = fs_gradient(rv[j], rps[i].refl,
			                                 cell, rps[i].panel);
		}

		/* fs terms: global */
		j = 0;
		group = rps[i].panel->group;
		while ( group != NULL ) {

			global_gradients[j] = -1.0;
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_TX);
			j++;

			global_gradients[j] = -xl / (kpred+zl);
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_TZ);
			j++;

			global_gradients[j] = 0.0;
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_RZ);
			j++;

			group = group->parent;
		}

		mille_add_measurement(mille,
		                      9, local_gradients,
		                      j, global_gradients, labels,
		                      fs_dev(&rps[i], image->detgeom),
		                      0.65*rps[i].panel->pixel_pitch);

		/* ss terms: local */
		for ( j=0; j<9; j++ ) {
			local_gradients[j] = ss_gradient(rv[j], rps[i].refl,
			                                 cell, rps[i].panel);
		}

		/* ss terms: global */
		j = 0;
		group = rps[i].panel->group;
		while ( group != NULL ) {

			global_gradients[j] = -1.0;
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_TY);
			j++;

			global_gradients[j] = -yl / (kpred+zl);
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_TZ);
			j++;

			global_gradients[j] = 0.0;
			labels[j] = mille_label(group->hierarchy_level, group->member_index, GPARAM_DET_RZ);
			j++;

			group = group->parent;
		}

		mille_add_measurement(mille,
		                      9, local_gradients,
		                      j, global_gradients, labels,
		                      ss_dev(&rps[i], image->detgeom),
		                      0.65*rps[i].panel->pixel_pitch);
	}
}


Mille *crystfel_mille_new(const char *outFileName,
                          int asBinary,
                          int writeZero)
{
	return mille_new(outFileName, asBinary, writeZero);
}


void crystfel_mille_free(Mille *m)
{
	mille_free(m);
}


void crystfel_mille_delete_last_record(Mille *m)
{
	mille_delete_last_record(m);
}


void crystfel_mille_write_record(Mille *m)
{
	mille_write_record(m);
}
