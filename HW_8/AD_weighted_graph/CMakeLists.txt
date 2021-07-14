cmake_minimum_required (VERSION 2.6)
project (dijkstra)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O4")

add_executable (test_dijkstra main.c dijkstra.c weighted_graph.c priority_queue.c binheap.c total_orders.c)