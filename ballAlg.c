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

typedef struct _point {
    double *original;
    double *projection;
} point_t;

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

void get_furthest_points(point_t *pts, long set_size, double **a, double **b)
{
    double dist, max_distance = 0.0;

    /* Lock b as first point in set and find a */
    *b = pts[0].original;
    for (long i = 1; i < set_size; i++)
    {
        if ((dist = distance(*b, pts[i].original)) > max_distance)
        {
            *a = pts[i].original;
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (long i = 0; i < set_size; i++)
    {
        if ((dist = distance(*a, pts[i].original)) > max_distance)
        {
            *b = pts[i].original;
            max_distance = dist;
        }
    }
}

/* Subtracts p2 from p1 and saves in result */
void sub_points(double *p1, double *p2, double *result)
{
    for (long i = 0; i < n_dims; i++)
    {
        result[i] = p1[i] - p2[i];
    }
}

/* Adds p2 to p1 and saves in result */
void add_points(double *p1, double *p2, double *result)
{
    for (long i = 0; i < n_dims; i++)
    {
        result[i] = p1[i] + p2[i];
    }
}

/* Computes inner product of p1 and p2 */
double inner_product(double *p1, double *p2)
{
    double result = 0.0;

    for (long i = 0; i < n_dims; i++)
    {
        result += p1[i] * p2[i];
    }

    return result;
}

/* Multiplies p1 with constant */
void mul_point(double *p1, double constant)
{
    for (long i = 0; i < n_dims; i++)
    {
        p1[i] *= constant;
    }
}

/* Project p onto line ab and save in result */
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

/* Sort compare function */
int cmpfunc(const void *a, const void *b)
{
    point_t *p1 = (point_t*) a;
    point_t *p2 = (point_t*) b;

    for (long i = 0; i < n_dims; i++) {
        if (p1->projection[i] > p2->projection[i]) {
            return 1;
        } else if (p1->projection[i] < p2->projection[i]) {
            return -1;
        }
    }

    return 0;
}

/* Computes the median point of a set of points in a line */
long median(point_t *pts, long pts_size, double *median_pt)
{
    /* Sorts array */
    qsort(pts, pts_size, sizeof(point_t), cmpfunc);

    if (pts_size % 2 != 0)
    {
        memcpy(median_pt, pts[pts_size / 2].projection, sizeof(double) * n_dims);
    }
    else
    {
        mean(pts[pts_size / 2 - 1].projection, pts[pts_size / 2].projection, median_pt);
    }
        
    /* Return pivot */
    return pts_size / 2;
}

node_t *build_tree(point_t *pts, long set_size)
{
    double *a, *b;


    node_t *node = (node_t*) malloc(sizeof(node_t));
    assert(node);

    node->center = (double*) malloc(n_dims * sizeof(double));
    assert(node->center);

    node->id = current_id++;

    /* Return leaf */
    if (set_size == 1) {
        memcpy(node->center, pts[0].original, sizeof(double) * n_dims);
        node->radius = 0.0;
        node->L = NULL;
        node->R = NULL;


        return node;
    }

    /* Find a and b */
    get_furthest_points(pts, set_size, &a, &b);

    /* Project points onto ab */
    double *projections = (double*) malloc(set_size * n_dims * sizeof(double));
    for (long i = 0; i < set_size; i++) {
        pts[i].projection = &projections[i * n_dims];
        project(pts[i].original, a, b, pts[i].projection);
    }

#ifdef DEBUG
    printf("a = ");
    print_point(a, n_dims);
    printf("b = ");
    print_point(b, n_dims);
    if (n_dims != 2) {
        printf("Visualization only works in 2d, skipping...");
    }
    else
    {
        projOutputFile = fopen("projections.csv", "w");
        fprintf(projOutputFile, "x,y,projx,projy\n");
        for (long i = 0; i < set_size; i++)
        {
            fprintf(projOutputFile, "%f,%f,%f,%f\n", pts[i].original[0], pts[i].original[1], pts[i].projection[0], pts[i].projection[1]);
        }
        fclose(projOutputFile);
    }
#endif

    /* Find median point */
    long pivot = median(pts, set_size, node->center);
    free(projections);

#ifdef DEBUG
    printf("median = ");
    print_point(node->center, n_dims);
#endif

    /* Calculate radius */
    double distA = distance(node->center, pts[0].original);
    double distB = distance(node->center, pts[set_size - 1].original);
    node->radius = distA > distB ? distA : distB;

    /* Split */
    node->L = build_tree(pts, pivot);
    node->R = build_tree(&pts[pivot], set_size - pivot);

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

    exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);
    point_t *points = (point_t*) malloc(n_points * sizeof(point_t));

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

    for (long i = 0; i < n_points; i++)
        points[i].original = pts[i];

    node_t *root = build_tree(points, n_points);


    exec_time += omp_get_wtime();
    fprintf(stderr, "%.11f\n", exec_time);

    free(points);

    dump_tree(root);
    free_tree(root);

    free(*pts);
    free(pts);
}
