#define GASPI_EXT 1
#include <assert.h>
#include <string.h>
#include <GASPI.h>
#include <GASPI_Ext.h>
#include "graph.h"

#define SEGMENT_ID_N 0
#define SEGMENT_ID_GRAPH 1
#define SEGMENT_ID_LOCAL 2
#define SEGMENT_ID_DIST 3
#define SEGMENT_ID_PATH 4
#define SEGMENT_ID_HAS_CHANGE 5

int N; //number of vertices
Graph* graph; // adjacency list

void abort_with_error_message(const char* msg) {
    fprintf(stderr, "%s\n", msg);
    abort();
}

int print_result(bool has_negative_cycle, int *dist, int* path, int dest) {
    FILE *outputf = fopen("output.txt", "w");
    if (!outputf) {
        perror("Failed to open output file");
        return -1;
    }

    int* temp = (int *) malloc(N * sizeof(int));

    int i;
    int index;
    for (i=0, index=dest; index != -1; index=path[index], i++) {
        temp[i] = index;
    }

    fprintf(outputf, "Node/Parent:");
    int j;
    for (j=0; j<N; j++) {
        fprintf(outputf, " (%d, %d)", j, path[j]);
    }
    fprintf(outputf, "\n");

    int z;
    fprintf(outputf, "path: ");
    for (z=i-1; z >= 0; z--) {
        fprintf(outputf, "%d, ", temp[z]);
    }
    fprintf(outputf, "\n");

    if (!has_negative_cycle) {
        for (int i = 0; i < N; i++) {
            if (dist[i] > INF)
                dist[i] = INF;
            fprintf(outputf, "%d\n", dist[i]);
        }
        fflush(outputf);
    } else {
        fprintf(outputf, "FOUND NEGATIVE CYCLE!\n");
    }
    
    fclose(outputf);
    return 0;
}

int read_file(const char* filename) {
    char temp[100];
    FILE* inputf = fopen(filename, "r");
    if (inputf == NULL) {
        abort_with_error_message("ERROR OCCURRED WHILE READING INPUT FILE");
    }
    
    if (fscanf(inputf, "%d", &N) != 1) {
        abort_with_error_message("ERROR READING MATRIX SIZE FROM FILE");
    }
    
    assert(N < (1024 * 1024 * 20));
    
    graph = createGraph(N);
    if (graph == NULL) {
        abort_with_error_message("MEMORY ALLOCATION FAILED");
    }

    int o;

    int dest;
    int length;
    
    for (int i = 0; i < N; i++) {
        if (fscanf(inputf, "%d", &o) != 1) {
            abort_with_error_message("ERROR READING MATRIX SIZE FROM FILE");
        }

        for (int j = 0; j < o; j++) {
            if (fscanf(inputf, "%d,%d", &dest, &length) != 2) {
                sprintf(temp, "ERROR READING MATRIX ELEMENT FROM FILE, i: %d, j: %d", i, j);
                abort_with_error_message(temp);
            }

            addEdge(graph, i, dest, length);
        }
    }
    
    fclose(inputf);
    return 0;
}

void bellman_ford(gaspi_rank_t my_rank, gaspi_rank_t num_ranks, gaspi_group_t group, int n, Graph* graph, int *dist, int* path,int source) {
if (my_rank == 0) {
  printf(">>> Entered bellman_ford with N=%d, ranks=%d, source=%d\n",
         n, num_ranks, source);
  fflush(stdout);
}
    int loc_n = n;
    int loc_start, loc_end;
    Graph* loc_graph;
    int* loc_dist;
    int* loc_path;
    int* all_dists;
    int* all_paths;
    bool* all_has_change;

    // Broadcast n and total_edges
    gaspi_pointer_t seg0_ptr;
    gaspi_segment_ptr(SEGMENT_ID_N, &seg0_ptr);
    int* seg0_data = (int*)seg0_ptr;
    
    if (my_rank == 0) {
        seg0_data[0] = n;
        int total_edges = 0;
        for (int i = 0; i < n; i++) {
            total_edges += graph->array[i].size;
        }
        seg0_data[1] = total_edges;
    }

    gaspi_barrier(group, GASPI_BLOCK);

    if (my_rank != 0) {
        gaspi_read(SEGMENT_ID_N, 0, 0, SEGMENT_ID_N, 0, 2 * sizeof(int), 0, GASPI_BLOCK);
        gaspi_wait(0, GASPI_BLOCK);  // Wait for read to complete
    }
    gaspi_barrier(group, GASPI_BLOCK);

    loc_n = seg0_data[0];
    int total_edges = seg0_data[1];

    // Allocate and build graph
    gaspi_pointer_t graph_ptr;
    gaspi_segment_ptr(SEGMENT_ID_GRAPH, &graph_ptr);
    int* graph_data = (int*)graph_ptr;

    if (my_rank == 0) {
        int idx = 0;
        for (int i = 0; i < loc_n; i++) {
            graph_data[idx++] = graph->array[i].size;
            AdjListNode* node = graph->array[i].head;
            for (int j = 0; j < graph->array[i].size; j++) {
                graph_data[idx++] = node->dest;
                graph_data[idx++] = node->length;
                node = node->next;
            }
        }
    }

    gaspi_barrier(group, GASPI_BLOCK);
    if (my_rank != 0) {
        gaspi_read(SEGMENT_ID_GRAPH, 0, 0, SEGMENT_ID_GRAPH, 0, (loc_n + 2 * total_edges) * sizeof(int), 0, GASPI_BLOCK);
        gaspi_wait(0, GASPI_BLOCK);  // Wait for read to complete
    }
    gaspi_barrier(group, GASPI_BLOCK);

    loc_graph = createGraph(loc_n);
    int idx = 0;
    for (int i = 0; i < loc_n; i++) {
        int num_edges = graph_data[idx++];
        for (int j = 0; j < num_edges; j++) {
            int dest = graph_data[idx++];
            int weight = graph_data[idx++];
            addEdge(loc_graph, i, dest, weight);
        }
    }

    // Determine local range
    int ave = loc_n / num_ranks;
    loc_start = ave * my_rank;
    loc_end = (my_rank == num_ranks - 1) ? loc_n : ave * (my_rank + 1);
    printf("rank %d: loc_start=%d, loc_end=%d\n", my_rank, loc_start, loc_end);
fflush(stdout);

    // Allocate memory
    loc_dist = (int*)malloc(loc_n * sizeof(int));
    loc_path = (int*)malloc(loc_n * sizeof(int));
    all_dists = (int*)malloc(num_ranks * loc_n * sizeof(int));
    all_paths = (int*)malloc(num_ranks * loc_n * sizeof(int));
    all_has_change = (bool*)malloc(num_ranks * sizeof(bool));

    // Initialize arrays
    for (int i = 0; i < loc_n; i++) {
        loc_dist[i] = (i == source) ? 0 : INF;
        loc_path[i] = -1;
    }

    gaspi_barrier(group, GASPI_BLOCK);

    // Main Bellman-Ford loop
    for (int iter = 0; iter < loc_n - 1; iter++) {
        bool loc_has_change = false;

        // Relax edges in local range
        for (int u = loc_start; u < loc_end; u++) {
            AdjListNode* pCrawl = loc_graph->array[u].head;
            while (pCrawl != NULL) {
                int dest = pCrawl->dest;
                int length = pCrawl->length;
                
                // Only relax if we have a valid path to u
                if (loc_dist[u] != INF && loc_dist[u] + length < loc_dist[dest]) {
                    loc_dist[dest] = loc_dist[u] + length;
                    loc_path[dest] = u;
                    loc_has_change = true;
                }
                pCrawl = pCrawl->next;
            }
        }

        // 2) do a global sum of that 0/1 flag => if sum>0 then someone changed
        int  loc_flag    = loc_has_change ? 1 : 0;
        int  global_flag = 0;
        printf("rank %d before allreduce: loc_flag=%d\n", my_rank, loc_flag);
fflush(stdout);
        gaspi_allreduce(
            &loc_flag,         // send buffer
            &global_flag,      // recv buffer
           1,                 // count = 1 element
            GASPI_OP_SUM,      // sum across all ranks
            GASPI_TYPE_INT,    // 32-bit integer type
            group,
            GASPI_BLOCK
        );
        printf("rank %d after allreduce: global_flag=%d\n", my_rank, global_flag);
fflush(stdout);
        bool global_has_change = (global_flag > 0);

        if (!global_has_change) {
            if (my_rank == 0) {
               printf("Converged at iteration %d, stopping.\n", iter);
            }
            break;
        }

                
// 3) First, copy our new loc_dist[]/loc_path[] (and the flag) into LOCAL segment
        gaspi_pointer_t local_ptr;
        gaspi_segment_ptr(SEGMENT_ID_LOCAL, &local_ptr);
        int   *seg_dist = (int*) local_ptr;
        int   *seg_path = seg_dist + loc_n;
        bool  *seg_flag = (bool*)(seg_path + loc_n);

        // copy working arrays into the segment
        memcpy(seg_dist, loc_dist, loc_n * sizeof(int));
        memcpy(seg_path, loc_path, loc_n * sizeof(int));
        *seg_flag = loc_has_change;

        // now push LOCAL→DIST and LOCAL→PATH exactly as before
        gaspi_queue_id_t que = 0;
        for (gaspi_rank_t r = 0; r < num_ranks; ++r) {
          // distances
          gaspi_write(SEGMENT_ID_LOCAL,
                      0,
                      r,
                      SEGMENT_ID_DIST,
                      my_rank * loc_n * sizeof(int),
                      loc_n * sizeof(int),
                      que, GASPI_BLOCK);
          // paths
          gaspi_write(SEGMENT_ID_LOCAL,
                      loc_n * sizeof(int),
                      r,
                     SEGMENT_ID_PATH,
                      my_rank * loc_n * sizeof(int),
                      loc_n * sizeof(int),
                      que, GASPI_BLOCK);
        }
        gaspi_wait(que, GASPI_BLOCK);
        gaspi_barrier(group, GASPI_BLOCK);

        // 4) Read global distances and paths (now correctly populated)
        gaspi_pointer_t dist_ptr, path_ptr;
        gaspi_segment_ptr(SEGMENT_ID_DIST, &dist_ptr);
        gaspi_segment_ptr(SEGMENT_ID_PATH, &path_ptr);
        
        memcpy(all_dists, dist_ptr, num_ranks * n * sizeof(int));
        memcpy(all_paths, path_ptr, num_ranks * n * sizeof(int));

        // Update local arrays
        for (int j = 0; j < n; j++) {
            int min_val = INF;
            int min_rank = -1;
            
            for (int r = 0; r < num_ranks; r++) {
                int idx = r * n + j;
                if (all_dists[idx] < min_val) {
                    min_val = all_dists[idx];
                    min_rank = r;
                }
            }
            
            loc_dist[j] = min_val;
            loc_path[j] = all_paths[min_rank * n + j];
        }

        if (my_rank == 0) {
            printf("Iteration %d complete. Changes: %s\n", 
                   iter, global_has_change ? "YES" : "NO");
        }
        
        gaspi_barrier(group, GASPI_BLOCK);
    }

    // Copy results to root
    if (my_rank == 0) {
        memcpy(dist, loc_dist, loc_n * sizeof(int));
        memcpy(path, loc_path, loc_n * sizeof(int));
    }

    // Free memory
    free(loc_dist);
    free(loc_path);
    free(all_dists);
    free(all_paths);
    free(all_has_change);
    freeGraph(loc_graph);
}

int main(int argc, char **argv) {
    if (argc <= 3) {
        abort_with_error_message("not enough number of input arguments!");
    }
    fprintf(stderr, "rank %d: entering main()\n", getpid());  
    fflush(stderr);
    

    char* filename = argv[1];
    int source = atoi(argv[2]);
    int dest_main = atoi(argv[3]);

    int *dist;
    int* path;

    gaspi_config_t config;
    gaspi_config_get(&config);
    config.queue_size_max = 2048;
    gaspi_config_set(config);

    fprintf(stderr, "[%d] before gaspi_proc_init\n", getpid()); fflush(stderr);
    gaspi_proc_init(GASPI_BLOCK);
    fprintf(stderr, "[rank?] after gaspi_proc_init\n"); fflush(stderr);
    
    gaspi_rank_t my_rank;
    gaspi_rank_t num_ranks;
    gaspi_proc_rank(&my_rank);
    
    fprintf(stderr, "[%d] got my_rank=%d\n", getpid(), my_rank);
fflush(stderr);

    gaspi_proc_num(&num_ranks);
    gaspi_group_t group;
    gaspi_group_create(&group);
    fprintf(stderr, "[%d] after group_create\n", my_rank);
fflush(stderr);
    gaspi_group_add(group, GASPI_GROUP_ALL);
    fprintf(stderr, "[%d] after group_add\n", my_rank);
fflush(stderr);
    
    fprintf(stderr, "rank %d: initialized GASPI, now at segment-create\n", my_rank);
fflush(stderr);

    // Create segments
    gaspi_segment_create(SEGMENT_ID_N, 2 * sizeof(int), group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);
    // On rank 0, read the file and fill seg0_data[0]=N, seg0_data[1]=total_edges
    gaspi_pointer_t seg0_ptr;
    gaspi_segment_ptr(SEGMENT_ID_N, &seg0_ptr);
    int* seg0_data = (int*)seg0_ptr;
    if (my_rank == 0) {
        assert(read_file(filename) == 0);
        int total_edges = 0;
        for (int i = 0; i < N; i++) {
            total_edges += graph->array[i].size;
        }
        seg0_data[0] = N;
        seg0_data[1] = total_edges;
        dist = malloc(N * sizeof(int));
       path = malloc(N * sizeof(int));
    }

    // Barrier, then everyone reads both N and total_edges
    gaspi_barrier(group, GASPI_BLOCK);
    if (my_rank != 0) {
        gaspi_read(SEGMENT_ID_N, 0, 0,
                  SEGMENT_ID_N, 0, 2*sizeof(int),
                  0, GASPI_BLOCK);
       gaspi_wait(0, GASPI_BLOCK);
       // Now update local N and total_edges
        N = seg0_data[0];
        int total_edges = seg0_data[1];
        (void)total_edges;  // we'll use it below
    }
    gaspi_barrier(group, GASPI_BLOCK);

    // Now every rank has a correct N and total_edges in seg0_data[]
    int total_edges = seg0_data[1];

    // Calculate segment sizes
    gaspi_size_t graph_segment_size = (N + 2 * total_edges) * sizeof(int);
    gaspi_size_t local_segment_size = (2 * N * sizeof(int)) + sizeof(bool);
    gaspi_size_t dist_segment_size = num_ranks * N * sizeof(int);
    gaspi_size_t path_segment_size = num_ranks * N * sizeof(int);
    gaspi_size_t has_change_segment_size = num_ranks * sizeof(bool);

    gaspi_segment_create(SEGMENT_ID_GRAPH, graph_segment_size, group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);
    gaspi_segment_create(SEGMENT_ID_LOCAL, local_segment_size, group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);
    gaspi_segment_create(SEGMENT_ID_DIST, dist_segment_size, group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);
    gaspi_segment_create(SEGMENT_ID_PATH, path_segment_size, group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);
    gaspi_segment_create(SEGMENT_ID_HAS_CHANGE, has_change_segment_size, group, GASPI_BLOCK, GASPI_MEM_INITIALIZED);

    gaspi_time_t t1, t2;
    gaspi_time_get(&t1);

    bellman_ford(my_rank, num_ranks, group, N, graph, dist, path, source);

    gaspi_time_get(&t2);
    if (my_rank == 0) {
        printf("Computation time: %f seconds\n", (double)(t2 - t1));
        print_result(false, dist, path, dest_main);
        free(dist);
        free(path);
        freeGraph(graph);
    }

    gaspi_barrier(group, GASPI_BLOCK);
    gaspi_proc_term(GASPI_BLOCK);
    return 0;
}