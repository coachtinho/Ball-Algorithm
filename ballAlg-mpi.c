#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "gen_points.h"
#include <mpi.h>

int n_dims, n_procs, id;
long n_points, max_depth, diff;
MPI_Status status;

enum MESSAGES {
    PRINT = -1
};

enum TAGS {
    PTS = 1,
    ID = 2,
    CENTER = 3,
    INTER = 4
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

    for (int d = 0; d < n_dims - 1; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return dist;
}

double distance(double *pt1, double *pt2)
{
    double dist = 0.0;

    for (int d = 0; d < n_dims - 1; d++)
        dist += (pt1[d] - pt2[d]) * (pt1[d] - pt2[d]);
    return sqrt(dist);
}

void mean(double *pt1, double *pt2, double *mean)
{
    for (long i = 0; i < n_dims - 1; i++)
    {
        mean[i] = (pt1[i] + pt2[i]) / 2;
    }
}

void get_furthest_points(double **pts, long l, long r, double **a, double **b)
{
    long i;
    double dist, max_distance = 0.0;

    *b = pts[l];
    for (i = l + 1; i < r + 1; i++) {
        *b = pts[i][n_dims - 1] < (*b)[n_dims - 1] ? pts[i] : *b;
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

void distr_get_furthest_points(double **pts, MPI_Comm comm, long size, double *a, double *b) {
    long i;
    double dist, max_distance = 0.0;
    double possible_points[n_dims * n_procs]; 

    /* Find first point in initial set and send to leader */
    memcpy(b, pts[0], n_dims * sizeof(double));
    for (i = 1; i < size; i++) {
        if (pts[i][n_dims - 1] < b[n_dims - 1]) {
            memcpy(b, pts[i], n_dims * sizeof(double));
        }
    }
    MPI_Gather(b, n_dims, MPI_DOUBLE, possible_points, n_dims, MPI_DOUBLE, 0, comm);

    /* Lock b as first point in set and broadcast it */
    if (!id) {
        memcpy(b, possible_points, n_dims * sizeof(double));
        for (int p = 1; p < n_procs; p++) {
            if (possible_points[p * n_dims + n_dims - 1] < b[n_dims - 1]) {
                memcpy(b, &possible_points[p * n_dims], n_dims * sizeof(double));
            }
        }
    }
    MPI_Bcast(b, n_dims, MPI_DOUBLE, 0, comm);
    
    for (i = 0; i < size; i++)
    {
        if ((dist = quick_distance(b, pts[i])) >= max_distance)
        {
            memcpy(a, pts[i], n_dims * sizeof(double));
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Calculate real a at leader */
    MPI_Gather(a, n_dims, MPI_DOUBLE, possible_points, n_dims, MPI_DOUBLE, 0, comm);
    if (!id) {
        for (int p = 0; p < n_procs; p++) {
            if ((dist = quick_distance(b, &possible_points[p * n_dims])) >= max_distance)
            {
                memcpy(a, &possible_points[p * n_dims], n_dims * sizeof(double));
                max_distance = dist;
            }
        }
    }

    max_distance = 0.0;

    /* Broadcast a and calculate b */
    MPI_Bcast(a, n_dims, MPI_DOUBLE, 0, comm);
    
    for (i = 0; i < size; i++)
    {
        if ((dist = quick_distance(a, pts[i])) >= max_distance)
        {
            memcpy(b, pts[i], n_dims * sizeof(double));
            max_distance = dist;
        }
    }

    max_distance = 0.0;

    /* Calculate real b at leader */
    MPI_Gather(b, n_dims, MPI_DOUBLE, possible_points, n_dims, MPI_DOUBLE, 0, comm);
    if (!id) {
        for (int p = 0; p < n_procs; p++) {
            if ((dist = quick_distance(a, &possible_points[p * n_dims])) >= max_distance)
            {
                memcpy(b, &possible_points[p * n_dims], n_dims * sizeof(double));
                max_distance = dist;
            }
        }
    }

    /* Broadcast b */
    MPI_Bcast(b, n_dims, MPI_DOUBLE, 0, comm);
}

/* Subtracts p2 from p1 and saves in result */
void sub_points(double *p1, double *p2, double *result)
{
    long i;

    for (i = 0; i < n_dims - 1; i++)
    {
        result[i] = p1[i] - p2[i];
    }
}

/* Adds p2 to p1 and saves in result */
void add_points(double *p1, double *p2, double *result)
{
    long i;

    for (i = 0; i < n_dims - 1; i++)
    {
        result[i] = p1[i] + p2[i];
    }
}

/* Computes inner product of p1 and p2 */
double inner_product(double *p1, double *p2)
{
    long i;
    double result = 0.0;

    for (i = 0; i < n_dims - 1; i++)
    {
        result += p1[i] * p2[i];
    }

    return result;
}

/* Multiplies p1 with constant */
void mul_point(double *p1, double constant, double *result)
{
    long i;

    for (i = 0; i < n_dims - 1; i++)
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
    add_points(result, a, result);
}

#pragma endregion

#pragma region qsort

/* swap function from qsort implementation */
#define MEMSWAP(a, b, size)         \
{\
      size_t __size = (size);                                                      \
      char *__a = (char*) (a), *__b = (char*) (b);                                              \
      do                                                                      \
        {                                                                      \
          char __tmp = *__a;                                                      \
          *__a++ = *__b;                                                      \
          *__b++ = __tmp;                                                      \
        } while (--__size > 0);                                                      \
}

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

long qsort_partition(double *pts, double *projs, long l, long r) {
    double *pivot = &projs[(r + l) / 2 * n_dims];
    long i = l - 1;
    long j = r + 1;

    while (1) {
        do
        {
            i++;
        } while (less_than(&projs[i * n_dims], pivot));
        do
        {
            j--;
        } while (less_than(pivot, &projs[j * n_dims]));
        if (i >= j) {
            return j;
        }
        MEMSWAP(&projs[i * n_dims], &projs[j * n_dims], n_dims * sizeof(double));
        MEMSWAP(&pts[i * n_dims], &pts[j * n_dims], n_dims * sizeof(double));
    }
}

void quicksort(double *pts, double *projs, long l, long r) {
    if (l < r) {
        long p = qsort_partition(pts, projs, l, r);
        quicksort(pts, projs, l, p);
        quicksort(pts, projs, p + 1, r);
    }
}

#pragma endregion

#pragma region qselect

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

long q_select_partition(double **pts, double **projs, long l, long r, long pivotIndex)
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
        pivotIndex = q_select_partition(pts, projs, l, r, pivotIndex);
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

#pragma region distributed

#define SWAPDOUBLES(x, y)         \
    {                      \
        double temp; \
        temp = x; \
        x = y; \
        y = temp; \
    }

int cmpdoubles(const void *a, const void *b)
{
    double p1 = *(double*)a;
    double p2 = *(double*)b;
    if (p1 > p2) return 1;
    else if (p1 < p2) return -1;
    else return 0;
}

int cmppoints(const void *a, const void *b)
{
    double *p1 = (double*) a;
    double *p2 = (double*) b;
    if (less_than(p2, p1)) return 1;
    else return -1;
}

/* Parallel sorting by regular sampling */
long distr_sorting(double *pts, double *proj, long size, MPI_Comm comm, double **sort_proj, double **sort_pts) {
    double my_pivots[n_procs];
    double all_pivots[n_procs * n_procs];
    int counts[n_procs];
    int displacements[n_procs];
    int rec_counts[n_procs];
    int rec_displacements[n_procs];

    /* Sort own portion of the array */
    quicksort(pts, proj, 0, size - 1);

    /* Select pivots */
    for (long i = 0; i < n_procs; i++) {
        my_pivots[i] = proj[i * (size / n_procs) * n_dims];
    } 

    /* Collect pivots at leader */
    MPI_Gather(my_pivots, n_procs, MPI_DOUBLE, all_pivots, n_procs, MPI_DOUBLE, 0, comm);
    if (!id) {
        qsort(all_pivots, n_procs * n_procs, sizeof(double), cmpdoubles);
        for (long i = 1; i < n_procs; i++) {
            my_pivots[i - 1] = all_pivots[i * n_procs];
        }
    }

    /* Broadcast final pivots */
    MPI_Bcast(my_pivots, n_procs - 1, MPI_DOUBLE, 0, comm);

    /* Calculate counts and displacements and share them  */
    long current = 0, count = 0, sum = 0;
    for (long pivot = 0; pivot < n_procs - 1; pivot++) {
        while (current < size && proj[current * n_dims] < my_pivots[pivot]) {
            count += n_dims;
            current++;
        }
        displacements[pivot] = sum;
        counts[pivot] = count;
        sum += count;
        count = 0;
    }
    displacements[n_procs - 1] = sum;
    counts[n_procs - 1] = (size - current) * n_dims;
    MPI_Alltoall(counts, 1, MPI_INT, rec_counts, 1, MPI_INT, comm);
    sum = 0;
    for (long i = 0; i < n_procs; i++) {
        rec_displacements[i] = sum;
        sum += rec_counts[i];
    }

    /* Final distribution of sorted array */
    *sort_proj = (double*) malloc(sum * sizeof(double));
    assert(sort_proj);
    *sort_pts = (double*) malloc(sum * sizeof(double));
    assert(sort_pts);
    MPI_Alltoallv(proj, counts, displacements, MPI_DOUBLE, *sort_proj, rec_counts, rec_displacements, MPI_DOUBLE, comm);
    MPI_Alltoallv(pts, counts, displacements, MPI_DOUBLE, *sort_pts, rec_counts, rec_displacements, MPI_DOUBLE, comm);
    quicksort(*sort_pts, *sort_proj, 0, sum / n_dims - 1);

    return sum / n_dims;
}

long distr_find_center(double *projs, long sort_size, long distr_size, double *center, MPI_Comm comm) {
    /* Each processor searches its portion for right indexes and sends back to leader for broadcast */
    int n_centers = distr_size % 2 == 1 ? 1 : 2;
    long has_centers[n_centers];
    long center_indexes[n_centers];
    long base;

    /* Init values */
    for (int i = 0; i < n_centers; i++) {
        has_centers[i] = -1;
        center_indexes[i] = n_centers == 1 ? distr_size / 2 : distr_size / 2 - 1 + i;
    }

    /* Each processor looks for center points and sends them to leader if found */
    if (!id) {
        double centers[n_centers][n_dims];

        /* Tell next processor to start searching */
        MPI_Send(&sort_size, 1, MPI_LONG, 1, 1, comm);

        for (int j = 0; j < n_centers; j++) {
            if (center_indexes[j] >= 0 && center_indexes[j] < sort_size) {
                has_centers[j] = center_indexes[j];
                memcpy(centers[j], &projs[center_indexes[j] * n_dims], n_dims * sizeof(double));
            }
        }

        /* Accumulate centers at leader to calculate real center */
        memset(center, 0, n_dims * sizeof(double));
        for (int j = 0; j < n_centers; j++) {
            if (has_centers[j] == -1) MPI_Recv(centers[j], n_dims, MPI_DOUBLE, MPI_ANY_SOURCE, j, comm, &status);

            /* Fill center */
            for (int d = 0; d < n_dims; d++) {
                center[d] += centers[j][d];
            }
        }
        for (int d = 0; d < n_dims; d++) center[d] /= n_centers;
    } else {
        MPI_Recv(&base, 1, MPI_LONG, id - 1, id, comm, &status);
        long max = sort_size + base;

        /* Tell next processor to start searching */
        if (id < n_procs - 1) MPI_Send(&max, 1, MPI_LONG, id + 1, id + 1, comm);

        for (int j = 0; j < n_centers; j++) {
            if (center_indexes[j] >= base && center_indexes[j] < max) {
                has_centers[j] = center_indexes[j] - base;
                /* Send center to leader */
                MPI_Send(&projs[(center_indexes[j] - base) * n_dims], n_dims, MPI_DOUBLE, 0, j, comm);
            }
        }
    }
    MPI_Bcast(center, n_dims, MPI_DOUBLE, 0, comm);
    return has_centers[n_centers - 1];
}

#pragma endregion

node_t *create_node(long id) {
    node_t *node = (node_t*) malloc(sizeof(node_t));
    assert(node);
    node->id = id;
    node->center = (double*) malloc(n_dims * sizeof(double));
    assert(node->center);
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

void finish_tree(double **pts, node_t *nodes, double **projections, long l, long r, long node_id, long depth, long base_id)
{
    node_t *node = &nodes[node_id - base_id];
    
    node->id = node_id;
    node->radius = 0.0;

    /* It's a leaf */
    if (r - l == 0)
    {
        memcpy(node->center, pts[l], n_dims * sizeof(double));
        node->left = -1;
        node->right = -1;
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

    if (depth < max_depth || (depth == max_depth && omp_get_thread_num() < diff)) {
#pragma omp taskgroup
        {
#pragma omp task
            finish_tree(pts, nodes, projections, l, l + split_index, node->left, depth + 1, base_id);
#pragma omp task
            finish_tree(pts, nodes, projections, l + split_index + 1, r, node->right, depth + 1, base_id);
        }
    } else {
        finish_tree(pts, nodes, projections, l, l + split_index, node->left, depth + 1, base_id);
        finish_tree(pts, nodes, projections, l + split_index + 1, r, node->right, depth + 1, base_id);
    }
}

long build_tree(double **pts, MPI_Comm team, node_t **nodes, long my_set, long team_set, long node_id) {
    MPI_Comm_size(team, &n_procs);
    MPI_Comm_rank(team, &id);

    /* Alone in team, finish sequentially */
    if (n_procs == 1) {
        double *to_free = *pts;
        /* Allocate memory for projections */
        double **projections = (double **)malloc(my_set * sizeof(double *));
        assert(projections);
        double *proj = (double *)malloc(my_set * n_dims * sizeof(double));
        assert(proj);
        for (long i = 0; i < my_set; i++)
        {
            projections[i] = &proj[i * n_dims];
        }

        max_depth = (int)log2(omp_get_max_threads());
        diff = omp_get_max_threads() - (1 << max_depth);

        /* Allocate memory for nodes */
        node_t *node_arr = (node_t *)malloc((2 * my_set - 1) * sizeof(node_t));
        assert(nodes);
        double *centers = (double *)malloc((2 * my_set - 1) * n_dims * sizeof(double));
        assert(centers);

        for (long i = 0; i < 2 * my_set - 1; i++)
        {
            node_arr[i].center = &centers[i * n_dims];
        }

#pragma omp parallel
#pragma omp single
        {
#pragma omp task
            finish_tree(pts, node_arr, projections, 0, my_set - 1, node_id, 0, node_id);
        }
        *nodes = attach_node(*nodes, node_arr);
        free(projections);
        free(proj);
        free(to_free);
        free(pts);
        return 2 * my_set - 1;
    }

    /* Find a and b */
    double a[n_dims], b[n_dims];

    distr_get_furthest_points(pts, team, my_set, a, b);
    
    /* Compute common factors to all projections */
    double b_a[n_dims];
    sub_points(b, a, b_a);
    double denominator = inner_product(b_a, b_a);
    double common_factor[n_dims];
    mul_point(b_a, 1 / denominator, common_factor);

    /* Allocate memory for projections */
    double *proj = (double*) malloc(my_set * n_dims * sizeof(double));
    assert(proj);
    for (long i = 0; i < my_set; i++)
    {
        project(pts[i], a, b_a, common_factor, &proj[i * n_dims]);
    }

    /* Sort first coordinate of projections */
    double *sorted_projs, *sorted_pts;
    long sorted_set = distr_sorting(*pts, proj, my_set, team, &sorted_projs, &sorted_pts);
    free(*pts);
    free(proj);

    /* Find center projection */
    /* split is -1 if I don't have the center */
    /* index of first R point otherwise */
    double center[n_dims];
    long split = distr_find_center(sorted_projs, sorted_set, team_set, center, team);
    free(sorted_projs);
    
    /* Calculate radius */
    double max_distance = 0.0;
    double radius;
    for (long i = 0; i < sorted_set; i++)
    {
        double dist = distance(center, &sorted_pts[i * n_dims]);
        if (dist > max_distance)
        {
            max_distance = dist;
        }
    }
    MPI_Reduce(&max_distance, &radius, 1, MPI_DOUBLE, MPI_MAX, 0, team);

    /* Processor with center sends its id to leader so it can broadcast it */
    int center_proc = 0;
    if (!id && split == -1) {
        MPI_Recv(&center_proc, 1, MPI_INT, MPI_ANY_SOURCE, CENTER, team, &status); 
    } else if (split != -1) {
        MPI_Send(&id, 1, MPI_INT, 0, CENTER, team);
    }
    MPI_Bcast(&center_proc, 1, MPI_INT, 0, team);

    long left = node_id + 1;
    long right = node_id + 2 * (team_set / 2);

    /* Add new node to leader's list */
    if (!id) {
        node_t *node = create_node(node_id);
        node->left = left;
        node->right = right;
        node->radius = radius;
        memcpy(node->center, center, n_dims * sizeof(double));
        *nodes = attach_node(*nodes, node);
    }

    /* Split communicator */
    /* Processor with center is part of R group unless it is processor 0 */
    MPI_Comm new_team, inter;
    int new_id, new_procs;
    long new_team_set;
    long new_node_id;
    if (center_proc) {
        /* Center processor is on R group */
        new_team_set = team_set / 2 + (team_set % 2 && id >= center_proc);
        new_node_id = id < center_proc ? left : right;
        MPI_Comm_split(team, id >= center_proc, id, &new_team);
    } else {
        /* Center processor is on L group */
        new_team_set = team_set / 2 + (team_set % 2 && id > center_proc);
        new_node_id = id <= center_proc ? left : right;
        MPI_Comm_split(team, id > center_proc, id, &new_team);
    }
    MPI_Comm_rank(new_team, &new_id);
    MPI_Comm_size(new_team, &new_procs);
    /* Create intercommunicator */
    if (center_proc) {
        /* Center processor is on R group */
        if (id < center_proc) {
            MPI_Intercomm_create(new_team, 0, team, center_proc, INTER, &inter);
        } else {
            MPI_Intercomm_create(new_team, 0, team, 0, INTER, &inter);
        }
    } else {
        /* Center processor is on L group */
        if (id <= center_proc) {
            MPI_Intercomm_create(new_team, 0, team, center_proc + 1, INTER, &inter);
        } else {
            MPI_Intercomm_create(new_team, 0, team, 0, INTER, &inter);
        }
    }

    /* Calculate counts and displacements */
    int rec_count = 0;
    int send_counts[n_procs - new_procs];
    int send_displacements[n_procs - new_procs];
    if (id == center_proc) {
        if (center_proc) {
            /* Center processor is on R group */
            long remainder = split % (n_procs - new_procs);
            int sum = 0;
            for (long i = 0; i < n_procs - new_procs; i++) {
                send_counts[i] = (split / (n_procs - new_procs) + (i < remainder)) * n_dims;
                send_displacements[i] = sum;
                sum += send_counts[i];
            }
        } else {
            /* Center processor is on L group */
            long remainder = (sorted_set - split) % (n_procs - new_procs);
            int sum = 0;
            for (long i = 0; i < n_procs - new_procs; i++) {
                send_counts[i] = ((sorted_set - split) / (n_procs - new_procs) + (i < remainder)) * n_dims;
                send_displacements[i] = sum;
                sum += send_counts[i];
            }
        }
    }
    if (center_proc) {
        /* Center processor is on R group */
        if (id < center_proc) {
            MPI_Scatter(NULL, 0, MPI_INT, &rec_count, 1, MPI_INT, 0, inter);
        } else if (id == center_proc) {
            MPI_Scatter(send_counts, 1, MPI_INT, NULL, 0, MPI_INT, MPI_ROOT, inter);
        } else {
            MPI_Scatter(NULL, 0, MPI_INT, NULL, 0, MPI_INT, MPI_PROC_NULL, inter);
        }
    } else {
        /* Center processor is on L group */
        if (id > center_proc) {
            MPI_Scatter(NULL, 0, MPI_INT, &rec_count, 1, MPI_INT, 0, inter);
        } else if (id == center_proc) {
            MPI_Scatter(send_counts, 1, MPI_INT, NULL, 0, MPI_INT, MPI_ROOT, inter);
        } else {
            MPI_Scatter(NULL, 0, MPI_INT, NULL, 0, MPI_INT, MPI_PROC_NULL, inter);
        }
    }

    /* Allocate memory for new points */
    if (rec_count > 0) {
        sorted_pts = (double*) realloc(sorted_pts, (sorted_set * n_dims + rec_count) * sizeof(double));
        assert(sorted_pts);
    }

    /* Center processor scatters its points */
    if (center_proc) {
        /* Center processor is on R group */
        if (id < center_proc) {
            MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, &sorted_pts[sorted_set * n_dims], rec_count, MPI_DOUBLE, 0, inter);
        } else if (id == center_proc) {
            MPI_Scatterv(sorted_pts, send_counts, send_displacements, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MPI_ROOT, inter);
        } else {
            MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, inter);
        }
    } else {
        /* Center processor is on L group */
        if (id > center_proc) {
            MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, &sorted_pts[sorted_set * n_dims], rec_count, MPI_DOUBLE, 0, inter);
        } else if (id == center_proc) {
            MPI_Scatterv(&sorted_pts[split * n_dims], send_counts, send_displacements, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MPI_ROOT, inter);
        } else {
            MPI_Scatterv(NULL, NULL, NULL, MPI_DOUBLE, NULL, 0, MPI_DOUBLE, MPI_PROC_NULL, inter);
        }
    }

    if (center_proc) {
        /* Get rid of L points at center processor */
        if (id == center_proc && split > 0) {
            double *new = (double*) malloc((sorted_set - split) * n_dims * sizeof(double));
            assert(new);
            memcpy(new, &sorted_pts[split * n_dims], (sorted_set - split) * n_dims * sizeof(double));
            double *t = sorted_pts;
            sorted_pts = new;
            free(t);
            sorted_set -= split;
        }
    } else {
        /* Get rid of R points at center processor */
        if (id == center_proc) {
            double *new = (double*) malloc(split * n_dims * sizeof(double));
            assert(new);
            memcpy(new, sorted_pts, split * n_dims * sizeof(double));
            double *t = sorted_pts;
            sorted_pts = new;
            free(t);
            sorted_set = split;
        }
    }

    if (rec_count > 0) sorted_set += rec_count / n_dims;

    /* Reconstruct pointer array */
    pts = (double**) realloc(pts, sorted_set * sizeof(double*));
    assert(pts);
    for (long i = 0; i < sorted_set; i++) {
        pts[i] = &sorted_pts[i * n_dims];
    }


    MPI_Barrier(new_team);

    return build_tree(pts, new_team, nodes, sorted_set, new_team_set, new_node_id);
}

#pragma region print

void print_node(node_t *node)
{
    printf("%ld %ld %ld %lf",
           node->id,
           node->left,
           node->right,
           node->radius);

    for (long i = 0; i < n_dims - 1; i++)
    {
        printf(" %lf", node->center[i]);
    }
    printf(" \n");
    fflush(stdout);
}

void dump_tree(long n_nodes, node_t *nodes)
{
    for (long i = 0; i < n_nodes; i++)
        print_node(&nodes[i]);
    for (node_t *aux = nodes->next; aux != NULL; aux = aux->next)
        print_node(aux);
}

#pragma endregion print

int main(int argc, char *argv[])
{
    double exec_time = -omp_get_wtime();

    MPI_Init(&argc, &argv);

    if (argc != 4) {
        printf("Usage: %s <n_dims> <n_points> <seed>\n", argv[0]);
        exit(1);
    }

    n_dims = atoi(argv[1]);
    if (n_dims < 2) {
        printf("Illegal number of dimensions (%d), must be above 1.\n", n_dims);
        exit(2);
    }

    n_points = atol(argv[2]);
    if (n_points < 1) {
        printf("Illegal number of points (%ld), must be above 0.\n", n_points);
        exit(3);
    }

    /* Create communicator without excess processors */
    MPI_Comm comm;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_split(MPI_COMM_WORLD, id < n_points, id, &comm);

    if (!id) {
        printf("%d %ld\n", n_dims, 2 * n_points - 1);
    }

    double **pts = NULL;
    node_t *nodes = NULL;
    long n_nodes = 0;
    if (id < n_points) {
        MPI_Comm_rank(comm, &id);
        MPI_Comm_size(comm, &n_procs);

        /* Generate points */
        long pts_per_proc = n_points / n_procs;
        long remainder = n_points % n_procs;
        long to_consume = n_dims;
        long my_set;

        if (id < remainder) {
            my_set = pts_per_proc + 1;
            to_consume *= (pts_per_proc + 1) * id;
        }
        else {
            my_set = pts_per_proc;
            to_consume *= (pts_per_proc + 1) * remainder + (pts_per_proc) * (id - remainder);
        }

        sprintf(argv[2], "%ld", my_set);
        pts = get_points(argc, argv, &n_dims, &my_set, to_consume, 1);

        /* Build tree */
        n_nodes = build_tree(pts, comm, &nodes, my_set, n_points, 0);
    }
    
    MPI_Barrier(MPI_COMM_WORLD);

    exec_time += omp_get_wtime();

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &n_procs);

    if (!id) {
        fprintf(stderr, "%.1f\n", exec_time);
    }

    int send = PRINT, recv;
    /* print */
    if (n_procs < 2) {
        dump_tree(n_nodes, nodes);
    } else if (!id) {
        if (nodes) dump_tree(n_nodes, nodes);
        fflush(stdout);
        MPI_Send(&send, 1, MPI_INT, 1, 1, MPI_COMM_WORLD);
        MPI_Recv(&recv, 1, MPI_INT, n_procs - 1, 0, MPI_COMM_WORLD, &status);
    } else {
        MPI_Recv(&recv, 1, MPI_INT, id - 1, id, MPI_COMM_WORLD, &status);
        if (nodes) dump_tree(n_nodes, nodes);
        fflush(stdout);
        MPI_Send(&send, 1, MPI_INT, (id + 1) % n_procs, (id + 1) % n_procs, MPI_COMM_WORLD);
    }


    if (nodes) free_list(nodes);
    MPI_Finalize();
}
