#ifndef REFINE_H
#define REFINE_H

#include "proj.h"

int refine (REAL eps_alpha2, REAL eps_th2, int sheet_cnt, sheet sheets[MAX_SHEET_CNT], work **scratch);

#endif /* !defined(REFINE_H) */
