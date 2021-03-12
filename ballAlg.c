#include <omp.h>
#include <stdio.h>
#include "gen_points.h"

int n_dims;
long n_points;

/* TODO: Project goes here */

int main(int argc, char *argv[]) {
    double exec_time;

    exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);

    /* root = build_tree(); */

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.11f\n", exec_time);

    /* dump_tree(root); */
}
