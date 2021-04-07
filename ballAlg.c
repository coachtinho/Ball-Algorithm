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

void get_furthest_points(double **pts, long l, long r, double **a, double **b)
{

    long i;
    double dist, max_distance = 0.0;

    /* finds first point relative to the original set */
    *b = pts[l];
    for (i = l + 1; i < r + 1; i++)
    {
        *b = pts[i] < *b ? pts[i] : *b;
    }

    /* Lock b as first point in set and find a */
    for (i = l; i < r + 1; i++)
    {
        if ((dist = distance(*b, pts[i])) > max_distance)
        {
            *a = pts[i];
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (i = l; i < r + 1; i++)
    {
        if ((dist = distance(*a, pts[i])) > max_distance)
        {
            *b = pts[i];
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
void mul_point(double *p1, double constant, double *result)
{
    long i;

    for (i = 0; i < n_dims; i++)
    {
        result[i] = p1[i] * constant;
    }
}

/* Projects p onto ab */
void project(double *p, double *a, double *b_a, double *result)
{
    double *auxiliary = (double *)malloc(n_dims * sizeof(double));
    assert(auxiliary);
    double numerator, denominator, fraction;

    /* Numerator of formula */
    sub_points(p, a, result);
    numerator = inner_product(result, b_a);

    /* Denominator of formula */
    denominator = inner_product(b_a, b_a);

    fraction = numerator / denominator;
    mul_point(b_a, fraction, auxiliary);
    add_points(auxiliary, a, result);

    free(auxiliary);
}

/**************************************************************************************
 * QUICKSELECT
 * ***********************************************************************************/

#define SWAP(x, y)         \
    {                      \
        double *temp1 = x; \
        x = y;             \
        y = temp1;         \
    }

int less_than(double *p1, double *p2)
{
    for (long i = 0; i < n_dims; i++)
    {
        if (p1[i] < p2[i])
        {
            return 1;
        }
        else if (p1[i] > p2[i])
        {
            return 0;
        }
    }
    return 0;
}

long partition(double **pts, double **projs, long l, long r, long pivotIndex, long offset)
{

    /* this is a little slower but it should be fiiiiiiiiiiiiiiiiiiiiiine */
    double *pivotValue = projs[pivotIndex];

    SWAP(projs[pivotIndex], projs[r]);
    SWAP(pts[pivotIndex + offset], pts[r + offset]);
    long storeIndex = l;

    for (long i = l; i < r; i++)
    {
        if (less_than(projs[i], pivotValue))
        {
            SWAP(projs[storeIndex], projs[i]);
            SWAP(pts[storeIndex + offset], pts[i + offset]);
            storeIndex++;
        }
    }

    SWAP(projs[r], projs[storeIndex]);
    SWAP(pts[r + offset], pts[storeIndex + offset]);

    return storeIndex;
}

double *qselect(double **pts, double **projs, long l, long r, long k, long offset)
{
    /* This way of doing it uses less stack */
    while (1)
    {
        if (l == r)
        {
            return projs[l];
        }
        long pivotIndex = l + rand() % (r - l + 1);
        pivotIndex = partition(pts, projs, l, r, pivotIndex, offset);
        if (k == pivotIndex)
        {
            return projs[k];
        }
        else if (k < pivotIndex)
        {
            r = pivotIndex - 1;
        }
        else
        {
            l = pivotIndex + 1;
        }
    }
}

/**************************************************************************************
 * END OF QUICKSELECT
 * ***********************************************************************************/

/* Computes the median point of a set of points in a line */
long median(double **pts, double **projs, long l, long r, double *center_pt)
{
    long projs_size = (r - l + 1);
    long k = projs_size / 2;

    if (projs_size % 2 != 0)
    {
        memcpy(center_pt, qselect(pts, projs, 0, projs_size - 1, k, l), sizeof(double) * n_dims);
    }
    else
    {
        qselect(pts, projs, 0, projs_size - 1, k, l);

        /* Finds point immediately before kth point */
        double *current = projs[0];
        for (long i = 1; i < k; i++)
        {
            if (less_than(current, projs[i]))
            {
                current = projs[i];
            }
        }
        mean(current, projs[k], center_pt);
    }
    k--;
    return k;
}

node_t *build_tree(double **pts, long l, long r)
{
#ifdef DEBUG
    printf("l = %ld; r = %ld\n", l, r);
    for (long i = l; i < r + 1; i++)
    {
        print_point(pts[i], n_dims);
    }
#endif
    node_t *node = (node_t *)malloc(sizeof(node_t));
    assert(node);

    node->id = current_id++;
    node->radius = 0.0;

    node->center = (double *)malloc(n_dims * sizeof(double));
    assert(node->center);

    if (r - l == 0)
    {
        memcpy(node->center, pts[l], sizeof(double) * n_dims);
        node->radius = 0.0;
        node->L = NULL;
        node->R = NULL;
        return node;
    }
    double *a, *b;
    double *b_a = (double *)malloc(n_dims * sizeof(double));
    assert(b_a);

    get_furthest_points(pts, l, r, &a, &b);

    /* Compute b - a */
    sub_points(b, a, b_a);

#ifdef DEBUG
    printf("a = ");
    print_point(a, n_dims);
    printf("b = ");
    print_point(b, n_dims);
#endif

    /* Project points onto ab */
    double **projections = (double **)malloc((r - l + 1) * sizeof(double *));
    double *proj = (double *)malloc((r - l + 1) * n_dims * sizeof(double));
    for (long i = l; i < r + 1; i++)
    {
        projections[i - l] = &proj[(i - l) * n_dims];
        project(pts[i], a, b_a, projections[i - l]);
    }
    free(b_a);

#ifdef DEBUG
    for (long i = 0; i < r - l + 1; i++)
    {
        printf("projection = ");
        print_point(projections[i], n_dims);
    }
#endif

    /* Find median point and split; 2 in 1 GIGA FAST */
    long split_index = median(pts, projections, l, r, node->center);

#ifdef DEBUG
    for (long i = 0; i < r - l + 1; i++)
    {
        printf("sorted projection = ");
        print_point(projections[i], n_dims);
    }
    printf("center = ");
    print_point(center, n_dims);
#endif

    free(projections);
    free(proj);

    /* Compute radius */
    for (long i = l; i < r + 1; i++)
    {
        double dist = distance(node->center, pts[i]);
        if (dist > node->radius)
        {
            node->radius = dist;
        }
    }

    node->L = build_tree(pts, l, l + split_index);
    node->R = build_tree(pts, l + split_index + 1, r);

    return node;
}

void print_node(node_t *node)
{
    if (node->L)
        print_node(node->L);
    if (node->R)
        print_node(node->R);

    printf("%ld %ld %ld %lf",
           node->id,
           node->L ? node->L->id : -1,
           node->R ? node->R->id : -1,
           node->radius);

    for (long i = 0; i < n_dims; i++)
    {
        printf(" %lf", node->center[i]);
    }
    printf(" \n");
}

void dump_tree(node_t *root)
{
    printf("%d %ld\n", n_dims, current_id);
    print_node(root);
}

void free_tree(node_t *root)
{
    if (root->L)
        free_tree(root->L);
    if (root->R)
        free_tree(root->R);

    free(root->center);
    free(root);
}

int main(int argc, char *argv[])
{
    double exec_time;

    exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);
    double *to_free = *pts;

    node_t *root = build_tree(pts, 0, n_points - 1);

    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1f\n", exec_time);

    /* fastest way would be to print the tree as we build it */
    dump_tree(root);
    free_tree(root);

    free(to_free);
    free(pts);
}
