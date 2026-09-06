//
//  fp.c
//  WDSat
//
//  Prime field (F_p) engine - Phase 1: reader, term lattice, and
//  linear solving by Gaussian elimination. Nonlinear rows are rejected
//  (Phases 2-3 add substitution propagation + search).
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "config.h"
#include "fp.h"

static fp_t fp_p;
static int_t fp_n_vars;
static int_t fp_n_eqs;
static int_t fp_n_terms;
static int_t fp_rank;

/// row e, slot 0 = constant (rhs), slot t = coefficient of term t
static fp_t fp_matrix[__FP_MAX_EQ__][__ID_SIZE__ + 1];
static bool fp_alive[__FP_MAX_EQ__];
static int_t fp_pivot_col[__FP_MAX_EQ__];
static fp_t fp_sol[__ID_SIZE__ + 1];

/// term lattice: term id 1..fp_n_terms; 1..fp_n_vars are the original
/// variables (degree 1), higher ids are monomials
static int_t fp_term_deg[__ID_SIZE__ + 1];
static int_t fp_term_fac[__ID_SIZE__ + 1][__MAX_DEGREE__];

static char fp_err[256];

const char *fp_errno_msg(void) {
	return(fp_err);
}

static void fp_seterr(const char *m) {
	strncpy(fp_err, m, sizeof(fp_err) - 1);
	fp_err[sizeof(fp_err) - 1] = 0;
}

/// @brief field ops (p < 2^24, products < 2^48)
static inline fp_t fp_sub(fp_t a, fp_t b) { return(a >= b ? a - b : a + fp_p - b); }
static inline fp_t fp_mul(fp_t a, fp_t b) { return((fp_t) ((uint64_t) a * b % fp_p)); }

static fp_t fp_pow(fp_t a, uint64_t e) {
	fp_t r = 1;
	while(e) {
		if(e & 1) r = fp_mul(r, a);
		a = fp_mul(a, a);
		e >>= 1;
	}
	return(r);
}

static inline fp_t fp_inv(fp_t a) {
	assert(a != 0);
	return(fp_pow(a, (uint64_t) fp_p - 2));
}

/// @brief find (or create) the term id for the sorted factor list
static int_t fp_term_id(const int_t *fac, int_t deg) {
	int_t i, j;
	if(deg == 1) return(fac[0]);
	for(i = 1; i <= fp_n_terms; i++) {
		if(fp_term_deg[i] != deg) continue;
		for(j = 0; j < deg; j++)
			if(fp_term_fac[i][j] != fac[j]) break;
		if(j == deg) return(i);
	}
	assert(fp_n_terms < __ID_SIZE__);
	fp_n_terms++;
	fp_term_deg[fp_n_terms] = deg;
	for(j = 0; j < deg; j++) fp_term_fac[fp_n_terms][j] = fac[j];
	return(fp_n_terms);
}

static int cmp_int(const void *a, const void *b) {
	return(*(const int_t *) a - *(const int_t *) b);
}

bool fp_read(const char *fn) {
	FILE *f;
	char line[4096];
	char *tok, *save;
	int_t e, k, d;

	f = fopen(fn, "r");
	if(f == NULL) {
		fp_seterr("cannot open input file");
		return(false);
	}
	if(fgets(line, sizeof(line), f) == NULL) goto fail;
	tok = strtok_r(line, " \t\r\n", &save);
	if(tok == NULL || strcmp(tok, "p") != 0) goto bad_header;
	tok = strtok_r(NULL, " \t\r\n", &save);
	if(tok == NULL || strcmp(tok, "fp") != 0) goto bad_header;
	fp_n_vars = atoi(strtok_r(NULL, " \t\r\n", &save));
	fp_n_eqs = atoi(strtok_r(NULL, " \t\r\n", &save));
	fp_p = (fp_t) atoi(strtok_r(NULL, " \t\r\n", &save));
	if(fp_n_vars < 1 || fp_n_vars > __MAX_ID__) { fp_seterr("bad n_vars"); goto fail; }
	if(fp_n_eqs < 1 || fp_n_eqs > __FP_MAX_EQ__) { fp_seterr("bad n_eqs"); goto fail; }
	if(fp_p < 2 || fp_p > __FP_MAX_P__) { fp_seterr("bad prime (need 2 <= p < 2^24)"); goto fail; }

	// variables are terms 1..n_vars
	fp_n_terms = fp_n_vars;
	for(int_t v = 1; v <= fp_n_vars; v++) fp_term_deg[v] = 1;
	memset(fp_matrix, 0, sizeof(fp_matrix));

	e = 0;
	while(fgets(line, sizeof(line), f) != NULL) {
		tok = strtok_r(line, " \t\r\n", &save);
		if(tok == NULL || tok[0] == '#') continue;
		if(strcmp(tok, "e") != 0) continue;
		if(e >= fp_n_eqs) { fp_seterr("too many equations"); goto fail; }
		int_t c0 = 0;
		int_t c = 0;
		bool have_c = false;
		// optional constant marker: K <const>; otherwise the first token
		// is the first coefficient (strtok cannot push back)
		tok = strtok_r(NULL, " \t\r\n", &save);
		if(tok != NULL && strcmp(tok, "K") == 0) {
			tok = strtok_r(NULL, " \t\r\n", &save);
			if(tok == NULL) { fp_seterr("dangling K"); goto fail; }
			c0 = (int_t) (atoi(tok) % (int) fp_p);
		}
		else if(tok != NULL) {
			c = (int_t) (atoi(tok) % (int) fp_p);
			have_c = true;
		}
		// coefficient/term pairs until the terminating 0
		while(have_c || (tok = strtok_r(NULL, " \t\r\n", &save)) != NULL) {
			if(!have_c) {
				if(strcmp(tok, "0") == 0) break;
				c = (int_t) (atoi(tok) % (int) fp_p);
			}
			have_c = false;
			tok = strtok_r(NULL, " \t\r\n", &save);
			if(tok == NULL || strcmp(tok, "0") == 0) { fp_seterr("dangling coefficient"); goto fail; }
			int_t tid;
			if(tok[0] == '.') {
				// monomial: .2 a b | .3 a b c | .4 a b c d
				d = atoi(tok + 1);
				if(d < 2 || d > __MAX_DEGREE__) { fp_seterr("bad monomial degree"); goto fail; }
				int_t fac[__MAX_DEGREE__];
				for(k = 0; k < d; k++) {
					tok = strtok_r(NULL, " \t\r\n", &save);
					if(tok == NULL || strcmp(tok, "0") == 0) { fp_seterr("truncated monomial"); goto fail; }
					fac[k] = atoi(tok);
					if(fac[k] < 1 || fac[k] > fp_n_vars) { fp_seterr("monomial factor out of range"); goto fail; }
				}
				qsort(fac, d, sizeof(int_t), cmp_int);
				tid = fp_term_id(fac, d);
			}
			else {
				tid = atoi(tok);
				if(tid < 1 || tid > fp_n_vars) { fp_seterr("variable out of range"); goto fail; }
			}
			if(c != 0) fp_matrix[e][tid] = (fp_t) c;
		}
		fp_matrix[e][0] = (fp_t) ((fp_p - (fp_t) c0) % fp_p); // move constant to rhs
		fp_alive[e] = true;
		e++;
	}
	if(e != fp_n_eqs) { fp_seterr("equation count mismatch"); goto fail; }
	fclose(f);
	return(true);

bad_header:
	fp_seterr("bad header (expected: p fp <n_vars> <n_eqs> <p>)");
fail:
	if(f != NULL) fclose(f);
	return(false);
}

/// @brief in-place Gaussian elimination (reduced row echelon) over the
/// variable columns 1..fp_n_terms; rows with monomial terms (deg >= 2)
/// are not linear and are rejected before this is called.
/// @return false on contradiction (0 = nonzero)
static bool fp_linear_solve(void) {
	int_t col, row, e, t;
	fp_t piv;

	fp_rank = 0;
	for(col = 1; col <= fp_n_terms; col++) {
		// pivot search
		row = -1;
		for(e = fp_rank; e < fp_n_eqs; e++)
			if(fp_alive[e] && fp_matrix[e][col] != 0) { row = e; break; }
		if(row < 0) continue;
		if(row != fp_rank) {
			for(t = 0; t <= fp_n_terms; t++) {
				const fp_t tmp = fp_matrix[row][t];
				fp_matrix[row][t] = fp_matrix[fp_rank][t];
				fp_matrix[fp_rank][t] = tmp;
			}
			const bool tmp_a = fp_alive[row];
			fp_alive[row] = fp_alive[fp_rank];
			fp_alive[fp_rank] = tmp_a;
		}
		piv = fp_matrix[fp_rank][col];
		if(piv != 1) {
			const fp_t ipiv = fp_inv(piv);
			for(t = col; t <= fp_n_terms; t++) fp_matrix[fp_rank][t] = fp_mul(fp_matrix[fp_rank][t], ipiv);
			fp_matrix[fp_rank][0] = fp_mul(fp_matrix[fp_rank][0], ipiv);
		}
		// eliminate column from all other rows
		for(e = 0; e < fp_n_eqs; e++) {
			if(e == fp_rank || !fp_alive[e]) continue;
			const fp_t c = fp_matrix[e][col];
			if(c == 0) continue;
			for(t = col; t <= fp_n_terms; t++)
				fp_matrix[e][t] = fp_sub(fp_matrix[e][t], fp_mul(fp_matrix[fp_rank][t], c));
			fp_matrix[e][0] = fp_sub(fp_matrix[e][0], fp_mul(fp_matrix[fp_rank][0], c));
		}
		fp_pivot_col[fp_rank] = col;
		fp_rank++;
	}
	// zero rows with nonzero constant => contradiction
	for(e = fp_rank; e < fp_n_eqs; e++) {
		if(!fp_alive[e]) continue;
		bool all_zero = (fp_matrix[e][0] == 0);
		for(t = 1; t <= fp_n_terms; t++)
			if(fp_matrix[e][t] != 0) { all_zero = false; break; }
		if(all_zero) {
			fp_alive[e] = false;
			continue;
		}
		if(fp_matrix[e][0] != 0) return(false);
	}
	return(true);
}

bool fp_solve(void) {
	int_t e, t, v;

	// Phase 1: reject nonlinear rows
	for(e = 0; e < fp_n_eqs; e++)
		for(t = 1; t <= fp_n_terms; t++)
			if(fp_matrix[e][t] != 0 && fp_term_deg[t] >= 2) {
				fprintf(stderr, "wdsat_fp: nonlinear system not supported in this build (phase 2+)\n");
				return(false);
			}

	if(!fp_linear_solve()) {
		printf("UNSAT\n");
		printf("conf:1\n");
		return(false);
	}

	// RREF: pivot row e is  v + sum(c*u over free u) = rhs, so with the
	// free variables set to 0 the pivot takes the rhs value
	for(t = 1; t <= fp_n_terms; t++) fp_sol[t] = 0;
	for(e = 0; e < fp_rank; e++) fp_sol[fp_pivot_col[e]] = fp_matrix[e][0];
	for(v = 1; v <= fp_n_vars; v++)
		printf("%u\n", (unsigned) fp_sol[v]);
	printf("conf:0\n");
	return(true);
}
