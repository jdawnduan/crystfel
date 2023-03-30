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
#include <stdlib.h>

#include "cell-utils.h"
#include "fast_indexer.h"

#ifdef HAVE_FAST_INDEXER

struct fast_indexer_private_data {
    UnitCell *cellTemplate;
    struct fast_feedback_options opts;
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

    struct fast_indexer_private_data *prv_data = (struct fast_indexer_private_data *) ipriv;

    npk = image_feature_count(image->features);
    if ( npk < prv_data->opts.min_peaks )
        return 0;

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


    float cell[9];

    double cell_internal_double[9];

    cell_get_cartesian(prv_data->cellTemplate,
                       &cell_internal_double[0],&cell_internal_double[1],&cell_internal_double[2],
                       &cell_internal_double[3],&cell_internal_double[4],&cell_internal_double[5],
                       &cell_internal_double[6],&cell_internal_double[7],&cell_internal_double[8]);

    for (int i = 0; i < 9; i++)
        cell[i] = cell_internal_double[i] * 1e10;

    struct ffbidx_settings settings;
    settings.max_spots_to_index = prv_data->opts.max_peaks;
    settings.min_spots_for_solution = prv_data->opts.min_peaks;
    settings.threshold = prv_data->opts.threshold_for_solution;
    settings.output_cells = prv_data->opts.output_cells;

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

void *fast_indexer_prepare(IndexingMethod *indm, UnitCell *cell, struct fast_feedback_options *opts) {
    if ( cell == NULL ) {
        ERROR("Unit cell information is required for fast indexer.\n");
        return NULL;
    }

    struct fast_indexer_private_data *prv_data = (struct fast_indexer_private_data *) malloc(sizeof(struct fast_indexer_private_data));

    prv_data->cellTemplate = cell;
    prv_data->opts = *opts;

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

static void fast_indexer_show_help()
{
    printf("Parameters for the fast feedback indexing algorithm:\n"
           "     --fast-feedback-indexer-max-peaks\n"
           "                            Maximum number of peaks used for indexing.\n"
           "                            All peaks are used for refinement.\n"
           "                            Default: 250\n"
           "     --fast-feedback-indexer-min-peaks\n"
           "                            Maximum number of indexed peaks to accept solution.\n"
           "                            Default: 9\n"
           "     --fast-feedback-indexer-threshold\n"
           "                            Threshold to accept solution as indexed.\n"
           "                            Default: 0.02\n"
           "     --fast-feedback-indexer-output-cells\n"
           "                            Number of output cells.\n"
           "                            Default: 1\n"
    );
}


int fast_indexer_default_options(struct fast_feedback_options **opts_ptr)
{
    struct fast_feedback_options *opts;

    opts = malloc(sizeof(struct fast_feedback_options));
    if ( opts == NULL ) return ENOMEM;

    opts->max_peaks = 100;
    opts->min_peaks = 9;
    opts->threshold_for_solution = 0.02f;
    opts->output_cells = 1;
    *opts_ptr = opts;
    return 0;
}


static error_t fast_indexer_parse_arg(int key, char *arg, struct argp_state *state)
{
    struct fast_feedback_options **opts_ptr = state->input;
    int r;

    switch ( key ) {
        case ARGP_KEY_INIT :
            r = fast_indexer_default_options(opts_ptr);
            if ( r ) return r;
            break;

        case 1 :
            fast_indexer_show_help();
            return EINVAL;

        case 2 :
            if (sscanf(arg, "%u", &(*opts_ptr)->max_peaks) != 1) {
                ERROR("Invalid value for --fast-feedback-indexer-max-peaks\n");
                return EINVAL;
            }
            break;

        case 3 :
            if (sscanf(arg, "%u", &(*opts_ptr)->min_peaks) != 1) {
                ERROR("Invalid value for --fast-feedback-indexer-min-peaks\n");
                return EINVAL;
            }
            break;

        case 4 :
            if (sscanf(arg, "%f", &(*opts_ptr)->threshold_for_solution) != 1) {
                ERROR("Invalid value for --fast-feedback-indexer-threshold\n");
                return EINVAL;
            }
            if (((*opts_ptr)->threshold_for_solution <= 0.0f) || ((*opts_ptr)->threshold_for_solution > 1.0f)) {
                ERROR("Invalid value for --fast-feedback-indexer-threshold; must be in range 0.0-1.0\n");
                return EINVAL;
            }
            break;
        case 5 :
            if (sscanf(arg, "%u", &(*opts_ptr)->output_cells) != 1) {
                ERROR("Invalid value for --fast-feedback-indexer-output-cells\n");
                return EINVAL;
            }
            if (((*opts_ptr)->output_cells == 0) || ((*opts_ptr)->output_cells > 128)) {
                ERROR("Invalid value for --fast-feedback-indexer-output-cells; must be in range 1-128\n");
                return EINVAL;
            }
            break;
    }

    return 0;
}


static struct argp_option fast_indexer_options[] = {
        {"help-fast-feedback-indexer", 1, NULL, OPTION_NO_USAGE, "Show options for fast feedback indexing algorithm", 99},
        {"fast-feedback-indexer-max-peaks", 2, "ffbidx_maxn", OPTION_HIDDEN, NULL},
        {"fast-feedback-indexer-min-peaks", 3, "ffbidx_minn", OPTION_HIDDEN, NULL},
        {"fast-feedback-indexer-threshold", 4, "ffbidx_threshold", OPTION_HIDDEN, NULL},
        {"fast-feedback-indexer-output-cells", 5, "ffbidx_out_cells", OPTION_HIDDEN, NULL},
        {0}
};


struct argp fast_feedback_argp = { fast_indexer_options, fast_indexer_parse_arg,
                              NULL, NULL, NULL, NULL, NULL };
