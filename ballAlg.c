#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "gen_points.h"

int n_dims;
long n_points;
long current_id = 0;

FILE *ptsOutputFile;
FILE *projOutputFile;

typedef struct _node
{
    long id;
    double *center;
    double radius;
    struct _node *L;
    struct _node *R;
} node_t;

double distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for (int d = 0; d < n_dims; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return sqrt(dist);
}

void mean(double *pt1, double *pt2, double *mean)
{
    for (long i = 0; i < n_dims; i++)
    {
        mean[i] = (pt1[i] + pt2[i]) / 2;
    }
}

int cmpfunc(const void *a, const void *b)
{
    return (double) ** (double**) a > (double) ** (double**) b;
}

/* Computes the median point of a set of points in a line */
void median(double **pts, long pts_size, double *median_pt)
{
    double *to_sort[pts_size];

    /* Copies points to array for sorting */
    for (long i = 0; i < pts_size; i++)
    {
        to_sort[i] = pts[i];
    }

    /* Sorts array */
    qsort(to_sort, pts_size, sizeof(double*), cmpfunc);

    if (pts_size % 2 != 0)
    {
        memcpy(median_pt, to_sort[pts_size / 2], sizeof(double) * n_dims);
    }
    else
    {
        mean(to_sort[pts_size / 2 - 1], to_sort[pts_size / 2], median_pt);
    }
}

void get_furthest_points(double **pts, long *current_set, long set_size, long *a, long *b)
{
    long i;
    double dist, max_distance = 0.0;

    /* Lock b as first point in set and find a */
    *b = 0;
    for (i = 1; i < set_size; i++)
    {
        if ((dist = distance(pts[current_set[*b]], pts[current_set[i]])) > max_distance)
        {
            *a = i;
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (i = 0; i < set_size; i++)
    {
        if (i != *a && (dist = distance(pts[current_set[*a]], pts[current_set[i]])) > max_distance)
        {
            *b = i;
            max_distance = dist;
        }
    }
}

/* Subtracts p2 from p1 and saves in result */
void sub_points(double *p1, double *p2, double *result)
{
    long i;

    for (i = 0; i < n_dims; i++)
    {
        result[i] = p1[i] - p2[i];
    }
}

/* Adds p2 to p1 and saves in result */
void add_points(double *p1, double *p2, double *result)
{
    long i;

    for (i = 0; i < n_dims; i++)
    {
        result[i] = p1[i] + p2[i];
    }
}

/* Computes inner product of p1 and p2 */
double inner_product(double *p1, double *p2)
{
    long i;
    double result = 0.0;

    for (i = 0; i < n_dims; i++)
    {
        result += p1[i] * p2[i];
    }

    return result;
}

/* Multiplies p1 with constant */
void mul_point(double *p1, double constant)
{
    long i;

    for (i = 0; i < n_dims; i++)
    {
        p1[i] *= constant;
    }
}

void project(double *p, double *a, double *b, double *result) {
    double *auxiliary = (double*) malloc(n_dims * sizeof(double));
    assert(auxiliary);
    double numerator, denominator, fraction;

    /* Numerator of formula */
    sub_points(p, a, result);
    sub_points(b, a, auxiliary);
    numerator = inner_product(result, auxiliary);

    /* Denominator of formula */
    denominator = inner_product(auxiliary, auxiliary);

    fraction = numerator / denominator;
    mul_point(auxiliary, fraction);
    add_points(auxiliary, a, result);

    free(auxiliary); 
}

node_t *build_tree(double **pts, long *current_set, long set_size)
{
    long a, b, i, L_count = 0, R_count = 0;

    double *projected_coords = (double*) malloc(n_dims * set_size * sizeof(double));
    double **projected = (double**) malloc(set_size * sizeof(double*));

    double *center = (double*) malloc(n_dims * sizeof(double));
    node_t *node = (node_t*) malloc(sizeof(node_t));

    assert(projected_coords);
    assert(projected);
    assert(center);
    assert(node);

    if (set_size == 1) {
        for (long d = 0; d < n_dims; d++) {
            center[d] = pts[current_set[0]][d];
        }
        node->id = current_id++;
        node->center = center;
        node->radius = 0.0;
        node->L = NULL;
        node->R = NULL;

        free(projected);
        free(projected_coords);

        return node;
    }

    get_furthest_points(pts, current_set, set_size, &a, &b);

    /* Project points onto ab */
    for (i = 0; i < set_size; i++) {
        project(pts[current_set[i]], pts[current_set[a]], pts[current_set[b]], &projected_coords[i * n_dims]);
        projected[i] = &projected_coords[i * n_dims];
    }

#ifdef DEBUG
    printf("a = ");
    print_point(pts[current_set[a]], n_dims);
    printf("b = ");
    print_point(pts[current_set[b]], n_dims);
    if (n_dims != 2) {
        printf("Visualization only works in 2d, skipping...");
    }
    else
    {
        projOutputFile = fopen("projections.csv", "w");
        fprintf(projOutputFile, "x,y,projx,projy\n");
        for (i = 0; i < set_size; i++)
        {
            fprintf(projOutputFile, "%f,%f,%f,%f\n", pts[current_set[i]][0], pts[current_set[i]][1], projected[i][0], projected[i][1]);
        }
        fclose(projOutputFile);
    }
#endif

    /* Find median point */
    median(projected, set_size, center);
    node->id = current_id++;
    node->center = center;
    node->radius = 0.0;

#ifdef DEBUG
    printf("median = ");
    print_point(center, n_dims);
#endif

    /* Split */
    long *L_set = (long*) malloc(sizeof(long));
    long *R_set = (long*) malloc(sizeof(long));
    assert(L_set);
    assert(R_set);

    for (i = 0; i < set_size; i++) {
        if (projected[i][0] < center[0]) {
            L_count++;
            L_set = (long*) realloc(L_set, L_count * sizeof(long));
            L_set[L_count - 1] = current_set[i];
        } else {
            R_count++;
            R_set = (long*) realloc(R_set, R_count * sizeof(long));
            R_set[R_count - 1] = current_set[i];
        }
    }

    free(projected_coords);
    free(projected);

    node->L = build_tree(pts, L_set, L_count);
    free(L_set);
    node->R = build_tree(pts, R_set, R_count);
    free(R_set);

    return node;
}

void print_node(node_t *node) {
    if (node->L) print_node(node->L);
    if (node->R) print_node(node->R);

    printf("%ld %ld %ld %lf",
            node->id,
            node->L ? node->L->id : -1,
            node->R ? node->R->id : -1,
            node->radius);

    for (long i = 0; i < n_dims; i++) {
        printf(" %lf", node->center[i]);
    }
    printf("\n");
}

void dump_tree(node_t *root) {
    printf("%d %ld\n", n_dims, current_id);
    print_node(root);
}

void free_tree(node_t *root) {
    if (root->L) free_tree(root->L);
    if (root->R) free_tree(root->R);

    free(root->center);
    free(root);
}

int main(int argc, char *argv[])
{
    double exec_time;
    long *current_set;
    long i;

    exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);

#ifdef DEBUG
    if (n_dims != 2)
    {
        printf("Visualization only works in 2d, skipping...");
    }
    else
    {
        ptsOutputFile = fopen("pts.csv", "w");
        fprintf(ptsOutputFile, "x,y\n");
        long i;
        for (i = 0; i < n_points; i++)
        {
            fprintf(ptsOutputFile, "%f,%f\n", pts[i][0], pts[i][1]);
        }
        fclose(ptsOutputFile);
    }
#endif

    current_set = (long *)malloc(n_points * sizeof(long));
    assert(current_set);

    for (i = 0; i < n_points; i++)
        current_set[i] = i;

    node_t *root = build_tree(pts, current_set, n_points);
    free(current_set);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.11f\n", exec_time);

    /* fastest way would be to print the tree as we build it */
    dump_tree(root);
    free_tree(root);

    free(*pts);
    free(pts);
}
