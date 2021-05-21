#include <stdio.h>
#include <stdlib.h>

#define RANGE 10

void print_point(double *point, int n_dims) {
    int i;

    for (i = 0; i < n_dims - 1; i++) {
        printf("%lf,", point[i]);
    }
    printf("%lf\n", point[i]);
}

double **create_array_pts(int n_dims, long np)
{
    double *_p_arr;
    double **p_arr;

    _p_arr = (double *) malloc(n_dims * np * sizeof(double));
    p_arr = (double **) malloc(np * sizeof(double *));
    if((_p_arr == NULL) || (p_arr == NULL)){
        printf("Error allocating array of points, exiting.\n");
        exit(4);
    }

    for(long i = 0; i < np; i++)
        p_arr[i] = &_p_arr[i * n_dims];

    return p_arr;
}


void consume_rand(int n) {
    for (int i = 0; i < n; i++) {
        random();
    }
}

double **get_points(int argc, char *argv[], int *n_dims, long *np, long n_consumes)
{
    double **pt_arr;
    long i;
    int j;

    pt_arr = (double **) create_array_pts(*n_dims, *np);

    int seed = atoi(argv[3]);
    srandom(seed);

    consume_rand(n_consumes);

    for(i = 0; i < *np; i++)
        for(j = 0; j < *n_dims; j++)
            pt_arr[i][j] = RANGE * ((double) random()) / RAND_MAX;

#ifdef DEBUG
    for (i = 0; i < *np; i++)
        print_point(pt_arr[i], *n_dims);
#endif

    return pt_arr;
}
