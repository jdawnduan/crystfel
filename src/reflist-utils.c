/*
 * reflist-utils.c
 *
 * Utilities to complement the core reflist.c
 *
 * (c) 2006-2011 Thomas White <taw@physics.org>
 *
 * Part of CrystFEL - crystallography with a FEL
 *
 */


#include <stdio.h>


#include "reflist.h"
#include "cell.h"
#include "utils.h"
#include "reflist-utils.h"
#include "symmetry.h"


double *intensities_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	double *out = new_list_intensity();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double intensity = get_intensity(refl);

		get_indices(refl, &h, &k, &l);

		set_intensity(out, h, k, l, intensity);

	}

	return out;
}


double *phases_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	double *out = new_list_phase();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double phase = get_phase(refl);

		get_indices(refl, &h, &k, &l);

		set_phase(out, h, k, l, phase);

	}

	return out;

}


unsigned char *flags_from_list(RefList *list)
{
	Reflection *refl;
	RefListIterator *iter;
	unsigned char *out = new_list_flag();

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		set_flag(out, h, k, l, 1);

	}

	return out;

}


int check_list_symmetry(RefList *list, const char *sym)
{
	unsigned char *flags;
	Reflection *refl;
	RefListIterator *iter;

	flags = flags_from_list(list);

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		int j;
		int found = 0;
		signed int h, k, l;

		get_indices(refl, &h, &k, &l);

		for ( j=0; j<num_equivs(h, k, l, sym); j++ ) {

			signed int he, ke, le;
			get_equiv(h, k, l, &he, &ke, &le, sym, j);

			if ( abs(he) > INDMAX ) continue;
			if ( abs(le) > INDMAX ) continue;
			if ( abs(ke) > INDMAX ) continue;

			found += lookup_flag(flags, he, ke, le);

		}

		if ( found > 1 ) {
			free(flags);
			return 1;  /* Symmetry is wrong! */
		}

	}

	free(flags);

	return 0;
}


int find_equiv_in_list(RefList *list, signed int h, signed int k,
                       signed int l, const char *sym, signed int *hu,
                       signed int *ku, signed int *lu)
{
	int i;
	int found = 0;

	for ( i=0; i<num_equivs(h, k, l, sym); i++ ) {

		signed int he, ke, le;
		Reflection *f;
		get_equiv(h, k, l, &he, &ke, &le, sym, i);
		f = find_refl(list, he, ke, le);

		/* There must only be one equivalent.  If there are more, it
		 * indicates that the user lied about the input symmetry.
		 * This situation should have been checked for earlier by
		 * calling check_symmetry() with 'items' and 'mero'. */

		if ( (f != NULL) && !found ) {
			*hu = he;  *ku = ke;  *lu = le;
			return 1;
		}

	}

	return 0;
}


void write_reflections_to_file(FILE *fh, RefList *list, UnitCell *cell)
{
	Reflection *refl;
	RefListIterator *iter;

	fprintf(fh, "  h   k   l          I    phase   sigma(I) "
		     " 1/d(nm^-1)  counts  fs/px  ss/px\n");

	for ( refl = first_refl(list, &iter);
	      refl != NULL;
	      refl = next_refl(refl, iter) ) {

		signed int h, k, l;
		double intensity, esd_i, s;
		int red;
		double fs, ss;
		char res[16];

		get_indices(refl, &h, &k, &l);
		get_detector_pos(refl, &fs, &ss);
		intensity = get_intensity(refl);
		esd_i = get_esd_intensity(refl);
		red = get_redundancy(refl);

		if ( cell != NULL ) {
			s = resolution(cell, h, k, l);
			snprintf(res, 16, "%10.2f", s/1e9);
		} else {
			strcpy(res, "         -");
		}

		fprintf(fh,
		       "%3i %3i %3i %10.2f %s %10.2f  %s %7i %6.1f %6.1f\n",
		       h, k, l, intensity, "       -", esd_i, res, red,
		       fs, ss);

	}
}


int write_reflist(const char *filename, RefList *list, UnitCell *cell)
{
	FILE *fh;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "w");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open output file '%s'.\n", filename);
		return 1;
	}

	write_reflections_to_file(fh, list, cell);

	fclose(fh);

	return 0;
}


RefList *read_reflections_from_file(FILE *fh)
{
	char *rval = NULL;
	int first = 1;
	RefList *out;

	out = reflist_new();

	do {

		char line[1024];
		signed int h, k, l;
		float intensity, sigma, fs, ss;
		char phs[1024];
		char ress[1024];
		int cts;
		int r;
		Reflection *refl;

		rval = fgets(line, 1023, fh);
		if ( rval == NULL ) continue;
		chomp(line);

		if ( strcmp(line, REFLECTION_END_MARKER) == 0 ) return out;

		r = sscanf(line, "%i %i %i %f %s %f %s %i %f %f",
		           &h, &k, &l, &intensity, phs, &sigma, ress, &cts,
		           &fs, &ss);
		if ( (r != 10) && (!first) ) {
			reflist_free(out);
			return NULL;
		}

		first = 0;
		if ( r == 10 ) {

			double ph;
			char *v;

			refl = add_refl(out, h, k, l);
			set_int(refl, intensity);
			set_detector_pos(refl, fs, ss, 0.0);
			set_esd_intensity(refl, sigma);

			ph = strtod(phs, &v);
			if ( v != NULL ) set_ph(refl, ph);

			/* The 1/d value is actually ignored. */

		}

	} while ( rval != NULL );

	/* Got read error of some kind before finding PEAK_LIST_END_MARKER */
	return NULL;
}


RefList *read_reflections(const char *filename)
{
	FILE *fh;
	RefList *out;

	if ( filename == NULL ) {
		fh = stdout;
	} else {
		fh = fopen(filename, "r");
	}

	if ( fh == NULL ) {
		ERROR("Couldn't open input file '%s'.\n", filename);
		return NULL;
	}

	out = read_reflections_from_file(fh);

	fclose(fh);

	return out;
}
