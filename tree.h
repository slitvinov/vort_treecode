/* various maximums & minimums */
#define MAX_P 8
#define MAX_CHILD_CNT 8

/* what type of error bound is to be used */
#define POTENTIAL 0
#define VELOCITY 1
#define BOUND POTENTIAL

struct str_tree_cell {
  struct str_tree_cell *child[MAX_CHILD_CNT];
  vect exp_pt, width_2, mmnts[MAX_P][MAX_P][MAX_P];
  REAL r, wght_nrm[MAX_P + 1];
  work *nodes;
  long node_cnt;
  int level, mmnts_stat, child_cnt;
};

typedef struct str_tree_cell tree_cell;

tree_cell *make_tree(int level, work *nodes, long node_cnt, int node_min);
int depth(tree_cell *cell_ptr);
int tree_cell_cnt(tree_cell *cell_ptr);
REAL child_parent_ratio(tree_cell *cell_ptr);
void tree_dump(FILE *f_ptr, tree_cell *cell_ptr);
void free_tree(void);
