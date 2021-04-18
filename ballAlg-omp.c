#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "gen_points.h"

int n_dims;
long n_points;
long max_depth = 0;

FILE *ptsOutputFile;

typedef struct _node
{
    long id;
    double *center;
    double radius;
    struct _node *L;
    struct _node *R;
} node_t;

#pragma region math

double quick_distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for (int d = 0; d < n_dims; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

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
        if ((dist = quick_distance(*b, pts[i])) > max_distance)
        {
            *a = pts[i];
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Find b */
    for (i = l; i < r + 1; i++)
    {
        if ((dist = quick_distance(*a, pts[i])) > max_distance)
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
void project(double *p, double *a, double *b_a, double *common_factor, double *result)
{
    double product;

    sub_points(p, a, result);
    product = inner_product(result, common_factor);

    mul_point(b_a, product, result);
}

#pragma endregion

#pragma region qselect

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

long median_of_three(double **pts, double **projs, long l, long r)
{
    long m = (l + r) / 2;
    if (less_than(projs[r], projs[l]))
    {
        SWAP(projs[r], projs[l]);
        SWAP(pts[r], pts[l]);
    }
    if (less_than(projs[m], projs[l]))
    {
        SWAP(projs[m], projs[l]);
        SWAP(pts[m], pts[l]);
    }
    if (less_than(projs[r], projs[m]))
    {
        SWAP(projs[r], projs[m]);
        SWAP(pts[r], pts[m]);
    }
    return m;
}

long partition(double **pts, double **projs, long l, long r, long pivotIndex)
{
    double *pivotValue = projs[pivotIndex];

    SWAP(projs[pivotIndex], projs[r]);
    SWAP(pts[pivotIndex], pts[r]);
    long storeIndex = l;

    for (long i = l; i < r; i++)
    {
        if (less_than(projs[i], pivotValue))
        {
            SWAP(projs[storeIndex], projs[i]);
            SWAP(pts[storeIndex], pts[i]);
            storeIndex++;
        }
    }

    SWAP(projs[r], projs[storeIndex]);
    SWAP(pts[r], pts[storeIndex]);

    return storeIndex;
}

double *qselect(double **pts, double **projs, long l, long r, long k)
{
    /* This way of doing it uses less stack */
    while (1)
    {
        if (l == r)
        {
            return projs[l];
        }
        long pivotIndex = median_of_three(pts, projs, l, r);
        pivotIndex = partition(pts, projs, l, r, pivotIndex);
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

/* Computes the median point of a set of points in a line */
long median(double **pts, double **projs, long l, long r, double *center_pt)
{
    long projs_size = (r - l + 1);
    long k = projs_size / 2;

    if (projs_size % 2 != 0)
    {
        memcpy(center_pt, qselect(pts, projs, l, r, k + l), sizeof(double) * n_dims);
    }
    else
    {
        qselect(pts, projs, l, r, k + l);

        /* Finds point immediately before kth point */
        double *current = projs[l];
        for (long i = l + 1; i < k + l; i++)
        {
            if (less_than(current, projs[i]))
            {
                current = projs[i];
            }
        }
        mean(current, projs[k + l], center_pt);
    }
    k--;
    return k;
}

#pragma endregion

node_t *build_tree(double **pts, double **projections, node_t *nodes, long l, long r, long depth, long id)
{

#ifdef DEBUG
    double exec_time = 0.0;
    if (depth < max_depth)
    {
        exec_time = -omp_get_wtime();
    }
#endif

    node_t *node = &nodes[id];

    node->id = id;
    node->radius = 0.0;

    if (r - l == 0)
    {
        memcpy(node->center, pts[l], sizeof(double) * n_dims);
        node->radius = 0.0;
        node->L = NULL;
        node->R = NULL;
        return node;
    }
    double *a, *b;

    get_furthest_points(pts, l, r, &a, &b);

    /* Compute common factors to all projections */
    double b_a[n_dims];
    sub_points(b, a, b_a);
    double denominator = inner_product(b_a, b_a);
    double common_factor[n_dims];
    mul_point(b_a, 1 / denominator, common_factor);

    /* Project points onto ab */
    for (long i = l; i < r + 1; i++)
    {
        project(pts[i], a, b_a, common_factor, projections[i]);
    }

    /* Find median point and split; 2 in 1 GIGA FAST */
    long split_index = median(pts, projections, l, r, node->center);

    add_points(node->center, a, node->center);

    /* Compute radius */
    for (long i = l; i < r + 1; i++)
    {
        double dist = distance(node->center, pts[i]);
        if (dist > node->radius)
        {
            node->radius = dist;
        }
    }

    if (depth >= max_depth)
    {
        node->L = build_tree(pts, projections, nodes, l, l + split_index, depth + 1, id + 1);
        node->R = build_tree(pts, projections, nodes, l + split_index + 1, r, depth + 1, id + 2 * (l + split_index - l + 1));
    }
    else
    {
#ifdef DEBUG
        exec_time += omp_get_wtime();
        printf("Hello from thread %d, launching two tasks at depth %ld. It took me %.1f to process this node\n", omp_get_thread_num(), depth, exec_time);
#endif
#pragma omp taskgroup
        {
#pragma omp task
            node->L = build_tree(pts, projections, nodes, l, l + split_index, depth + 1, id + 1);
#pragma omp task
            node->R = build_tree(pts, projections, nodes, l + split_index + 1, r, depth + 1, id + 2 * (l + split_index - l + 1));
        }
    }

    return node;
}

#pragma region print

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
    printf("%d %ld\n", n_dims, 2 * n_points - 1);
    print_node(root);
}

#pragma endregion print

int main(int argc, char *argv[])
{
    double exec_time = -omp_get_wtime();
    double **pts = get_points(argc, argv, &n_dims, &n_points);
    double *to_free = *pts;
    node_t *root;

    //TODO: FIX THIS FOR THE CASES WHERE OMP_NUM_THREADS IS NOT A POWER OF 2
    max_depth = (int)log2(omp_get_max_threads());

    /* Allocate memory for projections */
    double **projections = (double **)malloc(n_points * sizeof(double *));
    double *proj = (double *)malloc(n_points * n_dims * sizeof(double));
    for (long i = 0; i < n_points; i++)
    {
        projections[i] = &proj[i * n_dims];
    }

    /* Allocate memory for nodes */
    node_t *nodes = (node_t *)malloc((2 * n_points - 1) * sizeof(node_t));
    assert(nodes);
    double *centers = (double *)malloc((2 * n_points - 1) * n_dims * sizeof(node_t));
    assert(centers);

    for (long i = 0; i < 2 * n_points - 1; i++)
    {
        nodes[i].center = &centers[i * n_dims];
    }

#ifdef DEBUG
    printf("%ld\n", max_depth);
#endif

#pragma omp parallel
#pragma omp single
    {
#ifdef DEBUG
        printf("number of threads in team = %d\n", omp_get_num_threads());
#endif
#pragma omp task
        root = build_tree(pts, projections, nodes, 0, n_points - 1, 0, 0);
    }
    exec_time += omp_get_wtime();
    fprintf(stderr, "%.1f\n", exec_time);

    dump_tree(root);

    free(nodes);
    free(centers);
    free(projections);
    free(proj);
    free(to_free);
    free(pts);
}
