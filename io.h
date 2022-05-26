int readin(FILE *infile, REAL *del2, REAL *t0, REAL *dt0, REAL *vel0, REAL *vel,
           int *sheet_cnt, sheet sheets[MAX_SHEET_CNT], work **scratch);
int dump(REAL del2, REAL t, REAL dt0, REAL vel0, REAL vel, int sheet_cnt,
         sheet sheets[MAX_SHEET_CNT], FILE *outfile);
