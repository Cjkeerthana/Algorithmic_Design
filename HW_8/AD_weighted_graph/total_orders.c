#include <stdlib.h>
#include "weighted_graph.h"

int leq_float(const void *a, const void *b)
{
  return *((float*)a)<=*((float*)b);
}

int leq_int(const void *a, const void *b)
{
  return *((int*)a)<=*((int*)b);
}

int geq_int(const void *a, const void *b)
{
  return *((int*)a)>=*((int*)b);
}

int leq_node(const void *a, const void* b)
{
  return ((NODE*)a)->dist <= ((NODE*)b)->dist;
}