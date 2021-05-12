#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "gen_points.h"
#include <mpi.h>

int n_dims, n_procs, id;
long n_points;
MPI_Status status;

enum MESSAGES {
    TERMINATE = -1,
    PRINT = -2
};

enum TAGS {
    PTS = 1,
    ID = 2,
    DEPTH = 3
};

typedef struct _node
{
    long id;
    long left;
    long right;
    double *center;
    double radius;
    struct _node *next;
} node_t;

node_t *attach_node(node_t *head, node_t *new) {
    new->next = head;

    return new;
}

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
        for (long i = l + 1; i < l + k; i++)
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

node_t *create_node(long id) {
    node_t *node = (node_t*) malloc(sizeof(node_t));
    node->id = id;
    node->center = (double*) malloc(n_dims * sizeof(double));
    node->radius = 0.0;
    node->next = NULL;

    return node;
}

void free_node(node_t *node) {
    free(node->center);
    free(node);
}

void free_list(node_t *head) {
    node_t *aux = head;

    while (aux->next) {
        node_t *t = aux;
        aux = aux->next;
        free_node(t);
    }

    free_node(aux);
}

/**
 * Messaging protocol:
 * # of points
 * array of points
 * depth
 * id
 */
void send_workload(int target, double **pts, long size, long depth, long node_id) {
    // Send number of points
    MPI_Send(&size, 1, MPI_LONG, target, PTS, MPI_COMM_WORLD);

    // Send array of points
    for (long i = 0; i < size; i++) {
        MPI_Send(pts[i], n_dims, MPI_DOUBLE, target, i, MPI_COMM_WORLD);
    }

    // Send depth
    MPI_Send(&depth, 1, MPI_LONG, target, DEPTH, MPI_COMM_WORLD);

    // Send id
    MPI_Send(&node_id, 1, MPI_LONG, target, ID, MPI_COMM_WORLD);
}

int receive_workload(int sender, double ***pts, long *size, long *depth, long *node_id) {
    // Receive number of points
    MPI_Recv(size, 1, MPI_LONG, sender, PTS, MPI_COMM_WORLD, &status);
    // Check if termination message
    if (*size == TERMINATE) {
        return 1;
    } else {
        // Alocate space for points
        double *_p = (double*) malloc(n_dims * *size * sizeof(double));
        *pts = (double**) malloc(*size * sizeof(double*));
        for (long i = 0; i < *size; i++) {
            (*pts)[i] = &_p[i * n_dims];
        }
    }

    // Receive array of points
    for (long i = 0; i < *size; i++) {
        MPI_Recv((*pts)[i], n_dims, MPI_DOUBLE, sender, i, MPI_COMM_WORLD, &status);
    }

    // Receive depth
    MPI_Recv(depth, 1, MPI_LONG, sender, DEPTH, MPI_COMM_WORLD, &status);

    // Receive id
    MPI_Recv(node_id, 1, MPI_LONG, sender, ID, MPI_COMM_WORLD, &status);

    return 0;
}

void build_tree(double **pts, node_t **nodes, double **projections, long l, long r, long depth, long node_id)
{

    node_t *node = create_node(node_id);

    /* It's a leaf */
    if (r - l == 0)
    {
        memcpy(node->center, pts[l], n_dims * sizeof(double));
        node->left = -1;
        node->right = -1;
        *nodes = attach_node(*nodes, node);
        return;
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

    /* Since the projection skips summing a at the end it must be done here */
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

    node->left = node_id + 1;
    node->right = node_id + 2 * (split_index + 1);

    /* Add new node to list */
    *nodes = attach_node(*nodes, node);

    /* Calculate target processor */
    int target = id + (1 << depth);

    if (target < n_procs) {
        /* Send right child to other processor */ 
        long left = l + split_index + 1;
        long size = r - l - split_index;
        send_workload(target, pts + left, size, depth + 1, node->right);
    } else {
        build_tree(pts, nodes, projections, l + split_index + 1, r, depth + 1, node->right);
    }

    build_tree(pts, nodes, projections, l, l + split_index, depth + 1, node->left);
}

#pragma region print

void print_node(node_t *node)
{
    printf("%ld %ld %ld %lf",
           node->id,
           node->left,
           node->right,
           node->radius);

    for (long i = 0; i < n_dims; i++)
    {
        printf(" %lf", node->center[i]);
    }
    printf(" \n");
}

void dump_tree(node_t *nodes)
{
    for (node_t *aux = nodes; aux != NULL; aux = aux->next)
        print_node(aux);
}

#pragma endregion print

int main(int argc, char *argv[])
{
    double exec_time = -omp_get_wtime();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);
    /* long max_depth = (int)log2(n_procs); */
    double **pts = NULL;
    double *to_free = NULL;
    node_t *nodes = NULL;
    double **projections = NULL;
    double *proj = NULL;
    long depth, node_id;


    if (!id) {
        pts = get_points(argc, argv, &n_dims, &n_points);

        // Broadcast number of dimensions
        MPI_Bcast(&n_dims, 1, MPI_INT, 0, MPI_COMM_WORLD);

        to_free = *pts;
        depth = 0;
        node_id = 0;

        /* Allocate memory for projections */
        projections = (double **)malloc(n_points * sizeof(double *));
        proj = (double *)malloc(n_points * n_dims * sizeof(double));
        for (long i = 0; i < n_points; i++)
        {
            projections[i] = &proj[i * n_dims];
        }

        build_tree(pts, &nodes, projections, 0, n_points - 1, depth, node_id);
    } else {
        MPI_Bcast(&n_dims, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // Calculate processor that will send workload
        int sender = id - (1 << (int)log2(id));
        receive_workload(sender, &pts, &n_points, &depth, &node_id);
        to_free = *pts;
        /* Allocate memory for projections */
        projections = (double **)malloc(n_points * sizeof(double *));
        proj = (double *)malloc(n_points * n_dims * sizeof(double));
        for (long i = 0; i < n_points; i++)
        {
            projections[i] = &proj[i * n_dims];
        }

        build_tree(pts, &nodes, projections, 0, n_points - 1, depth, node_id);
    }

    if (!id) {
        // Send termination message if there are idle processors
        for (int i = n_points; i < n_procs; i++) {
            int message = TERMINATE;
            MPI_Send(&message, 1, MPI_INT, i, PTS, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    exec_time += omp_get_wtime();
    if (!id) {
        fprintf(stderr, "%.1f\n", exec_time);
        printf("%d %ld\n", n_dims, 2 * n_points - 1);
    }
    fflush(stdout);

    int send = PRINT, recv;
    /* print */
    if (n_procs < 2) {
        dump_tree(nodes);
    } else if (!id) {
        MPI_Send(&send, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
        MPI_Recv(&recv, 1, MPI_INT, n_procs - 1, 0, MPI_COMM_WORLD, &status);
        dump_tree(nodes);
    } else {
        MPI_Recv(&recv, 1, MPI_INT, id - 1, id, MPI_COMM_WORLD, &status);
        dump_tree(nodes);
        MPI_Send(&send, 1, MPI_INT, (id + 1) % n_procs, (id + 1) % n_procs, MPI_COMM_WORLD);
    }


    if (nodes) free_list(nodes);
    if (projections) free(projections);
    if (proj) free(proj);
    if (to_free) free(to_free);
    if (pts) free(pts);
    MPI_Finalize();
}
