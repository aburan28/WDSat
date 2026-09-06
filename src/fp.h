//
//  fp.h
//  WDSat
//
//  Prime field (F_p) engine: polynomial systems over F_p for
//  Weil-descent ECDLP. See docs/FPRIME_ENGINE.md.
//

#ifndef __fp_h__
#define __fp_h__

#include <stdint.h>
#include <stdbool.h>

#include "wdsat_utils.h"

/// max prime supported (p < 2^24 keeps products in 48 bits)
#define __FP_MAX_P__ (1U << 24)
/// max number of equations
#define __FP_MAX_EQ__ 4096
/// journal arena cap in bytes
#define __FP_JOURNAL_CAP__ (256U * 1024U * 1024U)

/// field element
typedef uint32_t fp_t;

/// @fn bool fp_read(const char *fn)
/// @brief read an .fp file (header "p fp <n_vars> <n_eqs> <p>"), build the
/// term lattice and the row matrix. Returns false with a readable reason
/// from fp_errno_msg.
bool fp_read(const char *fn);

/// @fn bool fp_solve(void)
/// @brief solve the loaded system: unit propagation + DFS over free
/// variables. Prints the model (one value per line) or UNSAT, plus a
/// "conf:<n>" line, mirroring the GF(2) engine's output style.
bool fp_solve(void);

/// @brief last error message (fp_read only)
const char *fp_errno_msg(void);

#endif
