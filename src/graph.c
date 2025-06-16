#include "graph.h"

void freeGraph(Graph* graph) {
    if (!graph) return;
    
    if (graph->array) {
        for (int v = 0; v < graph->V; v++) {
            AdjListNode* curr = graph->array[v].head;
            while (curr) {
                AdjListNode* next = curr->next;
                free(curr);
                curr = next;
            }
        }
        free(graph->array);
    }
    free(graph);
}

AdjListNode* newAdjListNode(int dest, int length) {
    AdjListNode* newNode = (AdjListNode*) malloc(sizeof(AdjListNode));
    newNode->dest = dest;
    newNode->length = length;
    newNode->next = NULL;

    return newNode;
}

Graph* createGraph(int V) {
    Graph* graph = (Graph*) malloc(sizeof(Graph));
    graph->V = V;
    graph->array = (AdjList*) malloc(V * sizeof(AdjList));

    int i;
    for (i = 0; i < V; ++i) {
        graph->array[i].head = NULL;
        graph->array[i].size = 0;
    }

    return graph;
}

void addEdge(Graph* graph, int src, int dest, int length) {
    AdjListNode* newNode = newAdjListNode(dest, length);
    newNode->next = graph->array[src].head;
    graph->array[src].head = newNode;
    graph->array[src].size++;
}

Graph* copyGraph(Graph* graph) {
    int V = graph->V;

    Graph* result = createGraph(V);
    for (int v = 0; v < V; v++) {
        int size = graph->array[v].size;
        AdjListNode* pCrawl = graph->array[v].head;

        while (pCrawl) {
            addEdge(result, v, pCrawl->dest, pCrawl->length);
            pCrawl = pCrawl->next;
        }
    }

    return result;
}

void printGraph(Graph* graph) {
    int v;
    for (v = 0; v < graph->V; v++) {
        int size = graph->array[v].size;
        AdjListNode* pCrawl = graph->array[v].head;
        printf("\n Adjacency list of vertex %d, size: %d \n head ", v, size);
        while (pCrawl) {
            printf("-> dest: %d, length: %d \n", pCrawl->dest, pCrawl->length);
            pCrawl = pCrawl->next;
        }
        printf("\n");
    }
}
