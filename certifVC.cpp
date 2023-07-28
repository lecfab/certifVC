// certifVC: Quality certification for Vertex Cover
// Copyright (C) 2023 Fabrice Lécuyer (fabrice.lecuyer@lip6.fr)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 or later.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License
// for more details: https://www.gnu.org/licenses/

#include <iomanip>
#include <iostream>
#include <string>

#include <cmath>
#include <vector>

#include <algorithm> // std::random_shuffle
#include <cstdlib>   // std::rand, std::srand
#include <ctime>     // std::time

#include "algo/algo_vertexcover.h"
#include "utils/CLI11.h" // options parser
#include "utils/adjlist.h"
#include "utils/edgelist.h"
#include "utils/inout.h"
#include "utils/tools.h"

using namespace std;

int main(int argc, char **argv) {
  TimeBegin()

      CLI::App app{"Quality certification for Vertex Cover: finds a small "
                   "cover and certifies its quality with a lower-bound on the "
                   "unknown minimum. See README.md for more information."};

  string filename, sol_heur = "greedyVC", bound_heur = "cliqueLB";
  bool output_cover = false;
  app.add_option(
         "dataset", filename,
         "Text file: list of `a b` edges with nodes IDs ranging from 0 to N-1")
      ->required();
  app.add_option("-s,--solution", sol_heur,
                 "Solution-heuristic, the algorithm used to obtain a small vertex cover")
      ->check(CLI::IsMember({"none", "greedyVC", "greedyMatching", "comparison"}, CLI::ignore_case))
      ->capture_default_str();
  app.add_option("-b,--bounding", bound_heur,
                 "Lower-bounding-heuristic, the algorithm used to obtain a "
                 "lower-bound on the minimum cover size")
      ->check(CLI::IsMember({"none", "matchLB", "cliqueLB", "comparison"}, CLI::ignore_case))
      ->capture_default_str();
  app.add_flag("-c,--cover", output_cover,
               "Print nodes of the vertex cover in standard output (add >f to redirect into file f)")
      ->capture_default_str();

  CLI11_PARSE(app, argc, argv);

  ifstream file(filename);
  srand(unsigned(time(0))); // initialise randomness

  // --------------------------------------------------
  // -------- Store the graph in an edgelist ----------
  // --------------------------------------------------
  Info("Reading edgelist from file " << filename)
  Edgelist h = Edgelist(file);

  Info("Number of nodes: " << h.n)
  Info("Number of edges: " << h.e)
  TimeStep("Read")

  // --------------------------------------------------
  // -- Convert edgelist file into Adjlist structure --
  // --------------------------------------------------
  Info("Converting to adjacency list");
  Uadjlist g(h);
  TimeStep("Adjlist")

  // --------------------------------------------------
  // -------- Execute vertex cover algorithms ---------
  // --------------------------------------------------
  vector<ul> smallestVC;
  if (sol_heur == "greedyMatching") { smallestVC = greedyMatching(g); }
  if (sol_heur == "greedyVC") { smallestVC = greedyVC(g); }
  if (sol_heur == "comparison") { smallestVC = solution_heuristics_comparison(g); }
  TimeStep("SolHeuristic")

  ul higestLB = 0;
  if (bound_heur == "matchLB") { higestLB = matchLB(g); }
  if (bound_heur == "cliqueLB") { higestLB = cliqueLB(g); }
  if (bound_heur == "comparison") { higestLB = bounding_heuristics_comparison(g); }
  TimeStep("BoundHeuristic")

  if(bound_heur != "none" and sol_heur != "none") {
    certifVC(g, higestLB, smallestVC);
    TimeStep("Certif")
  }

  if(sol_heur != "none" and output_cover) {
    for(ul &u: smallestVC)
      cout << u << endl;
    TimeStep("Output")
  }

  TimeTotal()
  return 0;
}
