#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "gen_points.h"

int n_dims;
long n_points;

typedef struct _node {
    long center;
    long radius;
    struct _node *L;
    struct _node *R;
} node_t;

double distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for(int d = 0; d < n_dims; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return sqrt(dist);
}

void get_furthest_points(double **pts, long *current_set, long set_size, long *a, long *b) {
    long i;
    double dist, max_distance = 0.0;

    /* Lock b as first point in set and find a */
    *b = *current_set;
    for (i = 1; i < set_size; i++) {
        if ((dist = distance(pts[*b], pts[current_set[i]])) > max_distance) {
            *a = current_set[i];
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (i = 0; i < set_size; i++) {
        if (i != *a && (dist = distance(pts[*a], pts[current_set[i]])) > max_distance) {
            *b = current_set[i];
            max_distance = dist;
        }
    }
}

/* Subtracts p2 from p1 and saves in result */
void sub_points(double *p1, double *p2, double *result) {
    long i;

    for (i = 0; i < n_dims; i++) {
        result[i] = p1[i] - p2[i]; 
    }
}

/* Adds p2 to p1 and saves in result */
void add_points(double *p1, double *p2, double *result) {
    long i;

    for (i = 0; i < n_dims; i++) {
        result[i] = p1[i] + p2[i]; 
    }
}

/* Computes inner product of p1 and p2 */
double inner_product(double *p1, double *p2) {
    long i;
    double result = 0;

    for (i = 0; i < n_dims; i++) {
        result += p1[i] * p2[i];
    }

    return result;
}

/* Multiplies p1 with constant */
void mul_point(double *p1, double constant) {
    long i;

    for (i = 0; i < n_dims; i++) {
        p1[i] *= constant;
    }
}

double *project(double *p, double *a, double *b) {
    double *auxiliary1 = (double*) malloc(n_dims * sizeof(double));
    double *auxiliary2 = (double*) malloc(n_dims * sizeof(double));
    assert(auxiliary1);
    assert(auxiliary2);
    double numerator, denominator, fraction;

    /* Numerator of formula */
    sub_points(p, a, auxiliary1);
    sub_points(b, a, auxiliary2);
    numerator = inner_product(auxiliary1, auxiliary2);

    /* Denominator of formula */
    sub_points(b, a, auxiliary1);
    sub_points(b, a, auxiliary2);
    denominator = inner_product(auxiliary1, auxiliary2);

    fraction = numerator / denominator;
    sub_points(b, a, auxiliary1);
    mul_point(auxiliary1, fraction);
    add_points(auxiliary1, a, auxiliary2);

    free(auxiliary1); 
    return auxiliary2;
}

void build_tree(double **pts, long *current_set, long set_size) {
    long a, b, i;
    double **projected = (double**) malloc(set_size * sizeof(double*));
    assert(projected);

    get_furthest_points(pts, current_set, set_size, &a, &b);

    /* Project points onto ab */
    for (i = 0; i < set_size; i++) {
        projected[i] = project(pts[current_set[i]], pts[a], pts[b]);
        /* print_point(projected[i], n_dims); */
    }

    /* Ask professor if it's fine to have memory leaks
     * If it is, remove this part
     * */
    for (i = 0; i < set_size; i++) {
        free(projected[i]);
    }
    free(projected);

    printf("a = ");
    print_point(pts[a], n_dims);
    printf("b = ");
    print_point(pts[b], n_dims);
}

int main(int argc, char *argv[]) {
    double exec_time;
    long *current_set;
    long i;

    exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);

    current_set = (long*) malloc(n_points * sizeof(long));
    assert(current_set);

    for (i = 0; i < n_points; i++) current_set[i] = i;

    /* node_t *root = build_tree(pts, current_set); */
    build_tree(pts, current_set, n_points);
    free(current_set);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.11f\n", exec_time);

    /* dump_tree(root); */

    free(*pts);
    free(pts);
}
