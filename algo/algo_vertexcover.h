#ifndef ALGO_VERTEXCOVER_H
#define ALGO_VERTEXCOVER_H

#include "../utils/tools.h"
struct Adjlist;
struct Bheap;

ul matchLB(const Adjlist &g);
ul cliqueLB(const Adjlist &g);
ul bounding_heuristics_comparison(const Adjlist &g);
std::vector<ul> greedyMatching(const Adjlist &g);
std::vector<ul> greedyVC(const Adjlist &g);
std::vector<ul> solution_heuristics_comparison(const Adjlist &g);
double certifVC(const Adjlist &g, const ul &higestLB, const std::vector<ul> &smallestVC);

#endif
