#ifndef GRAPH_H
#define GRAPH_H

#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

#define INF INT_MAX / 2

typedef struct AdjListNode AdjListNode;

struct AdjListNode {
    int dest;
    int length;
    AdjListNode* next;
};

typedef struct {
    AdjListNode *head;
    int size;
} AdjList;

typedef struct {
    int V;
    AdjList* array;
} Graph;

void freeGraph(Graph* graph);
AdjListNode* newAdjListNode(int dest, int length);
Graph* createGraph(int V);
void addEdge(Graph* graph, int src, int dest, int length);
Graph* copyGraph(Graph* graph);
void printGraph(Graph* graph);

#endif