#ifndef GEN_POINTS_H
#define GEN_POINTS_H

double **get_points(int argc, char *argv[], int *n_dims, long *np, long n_consumes, int index_dim);
void print_point(double *point, int n_dims);

#endif
