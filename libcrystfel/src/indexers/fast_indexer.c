/*
 * fast_indexer.c
 *
 * Invoke the Fast Feedback Indexer library
 *
 * Copyright © 2023 Paul Scherrer Institute
 * Copyright © 2017-2021 Deutsches Elektronen-Synchrotron DESY,
 *                       a research centre of the Helmholtz Association.
 *
 * Authors:
 *   2023 Filip Leonarski <filip.leonarski@psi.ch>
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
#include <stdio.h>
#include <ffbidx/c-wrapper.h>

#include "cell-utils.h"
#include "fast_indexer.h"

#ifdef HAVE_FAST_INDEXER

struct fast_indexer_private_data {
    UnitCell *cellTemplate;
};

static void makeRightHanded(UnitCell *cell)
{
    // From xgandalf.c
    double ax, ay, az, bx, by, bz, cx, cy, cz;
    cell_get_cartesian(cell, &ax, &ay, &az, &bx, &by, &bz, &cx, &cy, &cz);

    if ( !right_handed(cell) ) {
        cell_set_cartesian(cell, -ax, -ay, -az, -bx, -by, -bz, -cx, -cy, -cz);
    }
}

int run_fast_indexer(struct image *image, void *ipriv) {
    int npk;
    int i;

    npk = image_feature_count(image->features);
    if ( npk < 5 ) return 0;

    float *x = (float *) calloc(npk, sizeof(float));
    float *y = (float *) calloc(npk, sizeof(float));
    float *z = (float *) calloc(npk, sizeof(float));

    for ( i=0; i<npk; i++ ) {

        struct imagefeature *f;
        double r[3];

        f = image_get_feature(image->features, i);
        if ( f == NULL ) {
            ERROR("Empty feature ???");
            continue;
        }

        detgeom_transform_coords(&image->detgeom->panels[f->pn],
                                 f->fs, f->ss, image->lambda,
                                 0.0, 0.0, r);
        x[i] = r[0] * 1e-10;
        y[i] = r[1] * 1e-10;
        z[i] = r[2] * 1e-10;
    }

    struct fast_indexer_private_data *prv_data = (struct fast_indexer_private_data *) ipriv;

    float cell[9];

    double cell_internal_double[9];

    cell_get_cartesian(prv_data->cellTemplate,
                       &cell_internal_double[0],&cell_internal_double[1],&cell_internal_double[2],
                       &cell_internal_double[3],&cell_internal_double[4],&cell_internal_double[5],
                       &cell_internal_double[6],&cell_internal_double[7],&cell_internal_double[8]);

    for (int i = 0; i < 9; i++)
        cell[i] = cell_internal_double[i] * 1e10;

    struct ffbidx_settings settings;
    settings.max_spots = 100;

    fast_feedback_crystfel(&settings, cell, x, y, z, npk);

    free(x);
    free(y);
    free(z);

    UnitCell *uc;

    uc = cell_new();

    cell_set_cartesian(uc,
                        cell[0] * 1e-10, cell[3] * 1e-10, cell[6] * 1e-10,
                        cell[1] * 1e-10, cell[4] * 1e-10, cell[7] * 1e-10,
                        cell[2] * 1e-10, cell[5] * 1e-10, cell[8] * 1e-10);

    makeRightHanded(uc);

    cell_set_lattice_type(uc, cell_get_lattice_type(prv_data->cellTemplate));
    cell_set_centering(uc, cell_get_centering(prv_data->cellTemplate));
    cell_set_unique_axis(uc, cell_get_unique_axis(prv_data->cellTemplate));

    if ( validate_cell(uc) ) {
        ERROR("fast_indexer: problem with returned cell!\n");
        cell_free(uc);
        return 0;
    }

    //cell_free(uc);
    //return 0;


    Crystal *cr = crystal_new();
    if ( cr == NULL ) {
        ERROR("Failed to allocate crystal.\n");
        return 0;
    }
    crystal_set_cell(cr, uc);
    crystal_set_det_shift(cr, 0, 0);
    image_add_crystal(image, cr);
    return 1;
}

void *fast_indexer_prepare(IndexingMethod *indm, UnitCell *cell) {
    if ( cell == NULL ) {
        ERROR("Unit cell information is required for fast indexer.\n");
        return NULL;
    }

    struct fast_indexer_private_data *prv_data = (struct fast_indexer_private_data *) malloc(sizeof(struct fast_indexer_private_data));

    prv_data->cellTemplate = cell;

    *indm &= INDEXING_METHOD_MASK | INDEXING_USE_CELL_PARAMETERS;
    return prv_data;
}

void fast_indexer_cleanup(void *pp) {
    free(pp);
}

const char *fast_indexer_probe(UnitCell *cell) {
    return "fast_indexer";
}

#endif