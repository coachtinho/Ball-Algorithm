#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "gen_points.h"

int n_dims;
long n_points;

FILE *ptsOutputFile;
FILE *projOutputFile;

typedef struct _node
{
    long center;
    long radius;
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

/* used in qsort */
typedef struct _sort
{
    double *pt;
    double distance;
} sort_t;

void mean(double *pt1, double *pt2, double *mean)
{
    for (int i = 0; i < n_dims; i++)
    {
        mean[i] = (pt1[i] + pt2[i]) / 2;
    }
}

int cmpfunc(const void *a, const void *b)
{
    return ((sort_t *)a)->distance - ((sort_t *)b)->distance;
}

/* Computes the median point of a set of points in a line */
void median(double **pts, long pts_size, long a, double *median_pt)
{
    sort_t to_sort[pts_size];

    /* Calculates distances to point a */
    for (int i = 0; i < pts_size; i++)
    {
        to_sort[i].pt = pts[i];
        to_sort[i].distance = distance(pts[i], pts[a]);
    }

    /* Sorts array */
    qsort(to_sort, pts_size, sizeof(sort_t), cmpfunc);

    if (pts_size % 2 != 0)
    {
        median_pt = to_sort[pts_size / 2].pt;
    }
    else
    {
        mean(to_sort[pts_size / 2 - 1].pt, to_sort[pts_size / 2].pt, median_pt);
    }
}

void get_furthest_points(double **pts, long *current_set, long set_size, long *a, long *b)
{
    long i;
    double dist, max_distance = 0.0;

    /* Lock b as first point in set and find a */
    *b = *current_set;
    for (i = 1; i < set_size; i++)
    {
        if ((dist = distance(pts[*b], pts[current_set[i]])) > max_distance)
        {
            *a = current_set[i];
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (i = 0; i < set_size; i++)
    {
        if (i != *a && (dist = distance(pts[*a], pts[current_set[i]])) > max_distance)
        {
            *b = current_set[i];
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

double *project(double *p, double *a, double *b)
{
    double *auxiliary1 = (double *)malloc(n_dims * sizeof(double));
    double *auxiliary2 = (double *)malloc(n_dims * sizeof(double));
    assert(auxiliary1);
    assert(auxiliary2);
    double numerator, denominator, fraction;

    /* Numerator of formula */
    sub_points(p, a, auxiliary1);
    sub_points(b, a, auxiliary2);
    numerator = inner_product(auxiliary1, auxiliary2);

    /* Denominator of formula */
    denominator = inner_product(auxiliary2, auxiliary2);

    fraction = numerator / denominator;
    mul_point(auxiliary2, fraction);
    add_points(auxiliary2, a, auxiliary1);

    free(auxiliary2);
    return auxiliary1;
}

void build_tree(double **pts, long *current_set, long set_size)
{
    long a, b, i;
    double **projected = (double **)malloc(set_size * sizeof(double *));
    assert(projected);

    get_furthest_points(pts, current_set, set_size, &a, &b);

    /* Project points onto ab */
    for (i = 0; i < set_size; i++)
    {
        projected[i] = project(pts[current_set[i]], pts[a], pts[b]);
    }

#ifdef DEBUG
    if (n_dims != 2)
    {
        printf("Visualization only works in 2d, skipping...");
    }
    else
    {
        projOutputFile = fopen("projections.csv", "w");
        fprintf(projOutputFile, "x,y,projx,projy\n");
        fprintf(projOutputFile, "%f,%f,%f,%f\n", pts[a][0], pts[a][1], pts[b][0], pts[b][1]);
        for (i = 0; i < set_size; i++)
        {
            fprintf(projOutputFile, "%f,%f,%f,%f\n", pts[current_set[i]][0], pts[current_set[i]][1], projected[i][0], projected[i][1]);
        }
        fclose(projOutputFile);
    }
#endif

    /* Find median point */
    double median_pt[n_dims];

    median(projected, set_size, a, median_pt);
    print_point(median_pt, n_dims);

    /* Split */

    /* Ask professor if it's fine to have memory leaks
     * If it is, remove this part
     * */
    for (i = 0; i < set_size; i++)
    {
        free(projected[i]);
    }
    free(projected);

    printf("a = ");
    print_point(pts[a], n_dims);
    printf("b = ");
    print_point(pts[b], n_dims);
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

    /* node_t *root = build_tree(pts, current_set); */
    build_tree(pts, current_set, n_points);
    free(current_set);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.11f\n", exec_time);

    /* fastest way would be to print the tree as we build it */
    /* dump_tree(root); */

    free(*pts);
    free(pts);
}
