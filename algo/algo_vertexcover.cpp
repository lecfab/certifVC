// See header for general documentation

#include "algo_vertexcover.h"
#include <cmath>
#include <map>
#include <set>
#include <iterator>     // std::insert_iterator
#include <random>
#include <algorithm>
#include "../utils/adjlist.h"
#include "../utils/heap.h"

using namespace std;

// --------------------------------------------------
// ----------- nodeGreedy largest-first -------------
// --------------------------------------------------
vector<ul> nodeGreedy_largestFirst(const Adjlist &g) {
  Bheap h(g.n);
	for (ul u = 0; u < g.n; ++u)
		h.insert(Keyvalue(u, g.n-g.get_degree(u)));

  vector<ul> cover; cover.reserve(g.n / 2);
	while (h.n > 0) {
     // node with highest degree is added to the cover
		Keyvalue kv = h.popmin();
    if(kv.val == g.n) continue; // node without edges
		ul u = kv.key;
    Debug("Found "<<u)
		cover.push_back(u);

    // its neighbours loose a degree and are deleted if they have no more uncovered edge
    for (auto &v : g.neigh_iter(u)) {
      if(!h.contains(v)) continue;
      h.update_increment(v);
      if(h.key_to_val(v) >= g.n) h.remove(v);
    }
	}

  Info("nodeGreedy_largestFirst: size " << cover.size())
	return cover;
}

// --------------------------------------------------
// ----------- edgeGreedy largest-first -------------
// --------------------------------------------------
vector<ul> edgeGreedy_largestFirst(const Adjlist &g) {
  Bheap h(g.n);
	for (ul u = 0; u < g.n; ++u)
		h.insert(Keyvalue(u, g.n-g.get_degree(u)));

  vector<ul> cover; cover.reserve(g.n / 2);
	while (h.n > 0) {
     // node with highest degree is added to the cover
		Keyvalue kv = h.popmin();
    if(kv.val == g.n) continue; // node without edges
		ul u1 = kv.key;
		cover.push_back(u1);

    // its neighbour with highest degree is added to the cover
    ul u2 = u1, deg_ranking_u2 = g.n;
    for (auto &v : g.neigh_iter(u1)) {
      if(!h.contains(v)) continue;

      if(h.key_to_val(v) < deg_ranking_u2) {
        u2 = v;
        deg_ranking_u2 = h.key_to_val(v);
      }
    }

    cover.push_back(u2);
    h.remove(u2);

    // their neighbours loose a degree and are deleted if they have no more uncovered edge
    for (auto &v : g.neigh_iter(u1)) {
      if(!h.contains(v)) continue;
      h.update_increment(v);
      if(h.key_to_val(v) >= g.n) h.remove(v);
    }
    for (auto &v : g.neigh_iter(u2)) {
      if(!h.contains(v)) continue;
      h.update_increment(v);
      if(h.key_to_val(v) >= g.n) h.remove(v);
    }
	}
  Info("edgeGreedy_largestFirst: size " << cover.size())

	return cover;
}


// --------------------------------------------------
// ----------- nodeGreedy largest-first -------------
// --------------------------------------------------
vector<ul> edgeGreedy_smallestFirst(const Adjlist &g) {
  Bheap h(g.n);
	for (ul u = 0; u < g.n; ++u)
		h.insert(Keyvalue(u, g.get_degree(u)));

  vector<ul> cover; cover.reserve(g.n / 2);
	while (h.n > 0) {
    // node with lowest degree is added to the cover
		Keyvalue kv = h.popmin();
    if(kv.val == 0) continue; // node without edges
		ul u1 = kv.key;
    Debug("Found "<<u1)
		cover.push_back(u1);

    // its neighbour with lowest degree is added to the cover
    ul u2 = u1, deg_ranking_u2 = g.n;
    for (auto &v : g.neigh_iter(u1)) {
      if(!h.contains(v)) continue;

      if(h.key_to_val(v) < deg_ranking_u2) {
        u2 = v;
        deg_ranking_u2 = h.key_to_val(v);
      }
    }
    cover.push_back(u2);
    h.remove(u2);

    // their neighbours loose a degree and are deleted if they have no more uncovered edge
    for (auto &v : g.neigh_iter(u1)) {
      if(!h.contains(v)) continue;
      h.update_decrement(v);
      if(h.key_to_val(v) == 0) h.remove(v);
    }
    for (auto &v : g.neigh_iter(u2)) {
      if(!h.contains(v)) continue;
      h.update_decrement(v);
      if(h.key_to_val(v) == 0) h.remove(v);
    }
	}

  Info("edgeGreedy_smallestFirst: size " << cover.size())
	return cover;
}




/*
// --------------------------------------------------
// ------------- nodeGreedy generalised -------------
// --------------------------------------------------
Parameters allow to choose between different priority functions and update rules

  For greedy vertex cover algorithms, one can choose:
    - the first_node strategy, or priority function to choose a node
    - the second_node strategy, which consists in replacing the selected node by one of its neighbours
    - whether update_deg, whereby the degrees are updated after each deletion

  The following parameters are available:
    - VC_ASC: selection by ascending degree (low-degre first)
    - VC_DESC: selection by descending degree (high-degre first)
    - VC_RAND: random selection generated at every request
    - VC_NONE: no selection of a second node

  For greedy clique cover, the same parameters apply except that there is no second_node selection.
  For pruning, the first_node strategies are the same, and an extra score parameter defines how nodes are selected:
    - VC_DEG: depending on their total degree
    - VC_NEIGH: depending on their number of neighbours in the cover
    - VC_RNEIGH: depending on their number of neighbours outside of the cover
*/

int VC_ASC=0, VC_DESC=1, VC_RAND=2, VC_NONE=3, VC_DEG=4, VC_NEIGH=5, VC_RNEIGH=6;

// ----------------------- greedy_cover ----------------------- //
vector<ul> nodeGreedy_general(const Adjlist &g, int first_node=VC_ASC, int second_node=VC_ASC, bool update_deg=false) {
  Bheap h(g.n);
  function<ul(ul)> priority; // define the priority order to select the next node: ASC
  priority = [&g](ul u) { return rand(); };
  if(first_node == VC_ASC)
    priority = [&g](ul u) { return g.get_degree(u); };
  if(first_node == VC_DESC)
    priority = [&g](ul u) { return g.n - g.get_degree(u); };

  function<ul(vector<pair<ul,ul>>&)> select_neigh; // define if and how to select a neighbour instead
  select_neigh = [](vector<pair<ul,ul>> &a) { return a[rand() % a.size()].second; };
  if(second_node == VC_ASC)
    select_neigh = [](vector<pair<ul,ul>> &a) { return (*(min_element(begin(a), end(a)))).second; };
  if(second_node == VC_DESC)
    select_neigh = [](vector<pair<ul,ul>> &a) { return (*(max_element(begin(a), end(a)))).second; };

  function<void(ul)> update; // define how to update priority of neighbours of inserted node
  update = [&h](ul u) { };
  if(update_deg) {
    if(first_node == VC_ASC)
    update = [&h](ul v) { h.update_decrement(v); };
    if(first_node == VC_DESC)
    update = [&h](ul v) { h.update_increment(v); };
  }

  vector<ul> degrees; degrees.reserve(g.n); // fill the priority queue
	for (ul u = 0; u < g.n; ++u) {
		degrees.push_back(g.get_degree(u));
    h.insert(Keyvalue(u, priority(u)));
  }

  vector<ul> cover; cover.reserve(g.n);
  for (ul i = 0; i < g.n; ++i) { // for each node w in the priority order
    ul w = h.getmin().key;
    if(degrees[w] == 0) {
      h.remove(w);
      continue;
    }

    ul u = w; // select u neighbour of w with second_node criteria
    if(second_node != VC_NONE) {
      vector<pair<ul,ul>> neighs;
      for (auto &v : g.neigh_iter(w)) {
        if(!h.contains(v) or degrees[v]==0) continue;
        neighs.push_back({degrees[v], v});
      }
      if(neighs.size()) u = select_neigh(neighs);
    }

    cover.push_back(u); // cover u
    h.remove(u);

    for (auto &v : g.neigh_iter(u)) { // its neighbours loose a degree
      if(degrees[v] == 0) continue;
      degrees[v] --;
      update(v);
    }
  }

  Info("nodeGreedy_general("<< first_node <<","<< second_node <<","<< update_deg <<"): size " << cover.size())
	return cover;
}


// --------------------------------------------------
// ----------- greedyPruning generalised ------------
// --------------------------------------------------
// Takes a cover and makes it minimal
vector<ul> greedyPruning_general(const Adjlist &g, const vector<ul> &cover, int first_node=VC_ASC, int score=VC_DEG) {
  // if a node has all its neighbours in the cover, remove it from the cover

  vector<ul> neighbours_in_cover(g.n, 0);
  for(auto &u : cover) {
    for(auto &v : g.neigh_iter(u))
      neighbours_in_cover[v] ++;
  }

  // store the nodes with fully covered neighbours sorted according to first_node rule
  vector<bool> has_all_neighbours_covered(g.n, false);
  vector<ul> removable; removable.reserve(cover.size());
  vector<ul> pruned_cover; pruned_cover.reserve(cover.size());
  for(auto &u : cover) {
    if(neighbours_in_cover[u] == g.get_degree(u)) {
      has_all_neighbours_covered[u] = true;
      removable.push_back(u);
    }
    else pruned_cover.push_back(u);
  }

  function<ul(ul)> priority;
  if(score == VC_DEG) priority = [&g](ul u) { return g.get_degree(u); };
  if(score == VC_NEIGH) priority = [&neighbours_in_cover](ul u) { return neighbours_in_cover[u]; };
  if(score == VC_RNEIGH) priority = [&g, &neighbours_in_cover](ul u) { return g.get_degree(u)-neighbours_in_cover[u]; };

  if(first_node == VC_ASC)  sort(removable.begin(), removable.end(), [&priority,&g, &neighbours_in_cover](ul u, ul v) { return priority(u) < priority(v); });
  if(first_node == VC_DESC) sort(removable.begin(), removable.end(), [&priority,&g, &neighbours_in_cover](ul u, ul v) { return priority(u) > priority(v); });
  if(first_node == VC_RAND) random_shuffle(removable.begin(), removable.end());

  for(auto &u : removable) {
    if(!has_all_neighbours_covered[u]) pruned_cover.push_back(u);
    else // u is not added to pruned_cover so all its neighbours have to be in it
      for(auto &v : g.neigh_iter(u)) has_all_neighbours_covered[v] = false;
  }

  Info("greedyPruning_general(" << first_node <<","<< score <<"): size "<< pruned_cover.size() <<", reduction: "<< (100-100*pruned_cover.size()/cover.size()) <<"%")
  return pruned_cover;
}


// --------------------------------------------------
// --------- greedyCliqueCover generalised ----------
// --------------------------------------------------
ul greedyCliqueCover_general(const Adjlist &g, int first_node=VC_ASC, bool update_deg=false) {
  Bheap h(g.n);

  function<ul(ul)> priority; // define the priority order to select the next node
  priority = [&g](ul u) { return rand(); };
  if(first_node == VC_ASC)
    priority = [&g](ul u) { return (g.get_degree(u)); };
  if(first_node == VC_DESC)
    priority = [&g](ul u) { return (g.n - g.get_degree(u)); };

  function<ul(vector<pair<ul,ul>>&)> select_neigh; // define if and how to select a neighbour instead
  select_neigh = [](vector<pair<ul,ul>> &a) { return a[rand() % a.size()].second; };
  if(first_node == VC_ASC)
    select_neigh = [](vector<pair<ul,ul>> &a) { return (*(min_element(begin(a), end(a)))).second; };
  if(first_node == VC_DESC)
    select_neigh = [](vector<pair<ul,ul>> &a) { return (*(max_element(begin(a), end(a)))).second; };

  function<void(ul)> update; // define how to update priority of neighbours of inserted node
  update = [&h](ul u) { };
  if(update_deg) {
    if(first_node == VC_ASC)
    update = [&h](ul v) { h.update_decrement(v); };
    if(first_node == VC_DESC)
    update = [&h](ul v) { h.update_increment(v); };
  }

  vector<ul> degrees; degrees.reserve(g.n); // fill the priority queue
	for (ul u = 0; u < g.n; ++u) {
		degrees.push_back(g.get_degree(u));
    h.insert(Keyvalue(u, priority(u)));
  }

  ul lower_bound = 0;
  vector<ul> cover; cover.reserve(g.n);
  for (ul i = 0; i < g.n; ++i) { // for each node u1 in the priority order
    if(h.n == 0) break;

    auto tmp = h.popmin();
    ul u1 = tmp.key;
    if(degrees[u1] == 0) continue;

    ul current_clique = 1;

    // if a node v has all its remaining neighbours in the clique, it does not need to be added
    ul fully_covered = u1, d_fully_covered = degrees[u1];

    set<ul> neigh;
    for (auto &v : g.neigh_iter(u1))
      if(h.contains(v) and degrees[v] > 0)
        neigh.insert(v);

    while(!neigh.empty()) {
      current_clique ++;

      // select neighbour of best priority
      ul u2 = u1;
      vector<pair<ul,ul>> neighs_deg;
      for (auto &v : neigh) neighs_deg.push_back({degrees[v], v});
      if(neighs_deg.size()) u2 = select_neigh(neighs_deg);
      else { Alert("No neighbours left for node "<< u1) }
      if(u2 == u1) { Alert("Identical selection "<< u1)}


      // add it to current clique
      if(degrees[u2] + current_clique - 2 < d_fully_covered) {
        cover.push_back(fully_covered);
        fully_covered = u2;
        d_fully_covered = degrees[u2] + current_clique - 2;
      }
      else cover.push_back(u2);
      if(!h.contains(u2)) { Alert("not in heap "<< u2) }
      h.remove(u2);

      // neighbours of u2 loose a degree and are deleted if they have no more uncovered edge
      for (auto &v : g.neigh_iter(u2)) {
        if(!h.contains(v) or degrees[v] == 0) continue;
        degrees[v] --;
        update(v);
      }

      // compute intersection with neighbours of u2
      set<ul> neigh2(g.neigh_beg(u2), g.neigh_end(u2));
      set<ul> neigh_tmp;
      for(auto &v : neigh)
        if(neigh2.count(v) and degrees[v] > 0)
          neigh_tmp.insert(v);
      neigh.swap(neigh_tmp);
    }

    // neighbours of u1 loose a degree and are deleted if they have no more uncovered edge
    for (auto &v : g.neigh_iter(u1)) {
      if(!h.contains(v) or degrees[v] == 0) continue;
      degrees[v] --;
      update(v);
    }

    // if u1 has all its remaining neighbours in clique, no need to select it
    if(d_fully_covered >= current_clique)
      cover.push_back(fully_covered);

    lower_bound += current_clique - 1;
	}

  Info("greedyCliqueCover_general("<< first_node <<","<< update_deg <<"): size " << cover.size() <<", bound "<< lower_bound)
  return lower_bound;
}


// --------------------------------------------------
// -- Check the coverage and minimality of a cover --
// --------------------------------------------------
bool check_cover(const Adjlist &g, const vector<ul> &cover, const bool require_minimal=false) {
  vector<bool> is_in_cover(g.n, false);
  for(ul u: cover)
    is_in_cover[u] = true;
  bool coverage = true, minimality = true;
  for(ul u=0; u<g.n; ++u) {
    if(is_in_cover[u]) { // check minimality: not all neighbours should be covered
      bool has_all_neighbours_covered = true;
      for (auto &v : g.neigh_iter(u))
        if(!is_in_cover[v]) has_all_neighbours_covered = false;
      if(has_all_neighbours_covered) { minimality = false; break; }
    }
    else { // node is not covered; check coverage: all neighbours should be covered
      for (auto &v : g.neigh_iter(u))
        if(!is_in_cover[v]) { coverage = false; break; }
    }
  }
  if(!coverage) { Alert("Invalid vertex cover") }
  if(require_minimal and !minimality) { Alert("Not minimal cover") }
  return minimality;
}



// --------------------------------------------------
// ------------ Aliasses for algorithms -------------
// --------------------------------------------------
ul matchLB(const Adjlist &g) {
  ul matching_lower_bound = edgeGreedy_smallestFirst(g).size() / 2;
  Info("matchLB: lower bound "<< matching_lower_bound)
  return matching_lower_bound;
}

ul cliqueLB(const Adjlist &g) {
  ul clique_lower_bound = greedyCliqueCover_general(g, VC_ASC, true);
  Info("cliqueLB: lower bound "<< clique_lower_bound)
  return clique_lower_bound;
}

ul bounding_heuristics_comparison(const Adjlist &g) {
  ul highestLB = matchLB(g);

  for(auto first_node: {VC_ASC, VC_DESC, VC_RAND})
  for(auto update_deg: {true, false}) {
    highestLB = max(highestLB, greedyCliqueCover_general(g, first_node, update_deg));
  }

  Info("highestLB: lower bound "<< highestLB)
  return highestLB;
}


vector<ul> greedyMatching(const Adjlist &g) {
  vector<ul> cover = edgeGreedy_largestFirst(g);
  vector<ul> post_cover = greedyPruning_general(g, cover);
  Info("greedyMatching: cover size "<< post_cover.size())
  return post_cover;
}

vector<ul> greedyVC(const Adjlist &g) {
  vector<ul> cover = nodeGreedy_largestFirst(g);
  vector<ul> post_cover = greedyPruning_general(g, cover);
  Info("greedyVC: cover size "<< post_cover.size())
  return post_cover;
}

vector<ul> solution_heuristics_comparison(const Adjlist &g) {
  vector<ul> smallestVC = greedyMatching(g);
  vector<ul> cover = greedyVC(g);
  vector<ul> post_cover;

  if(cover.size() < smallestVC.size()) smallestVC.swap(cover);

  for(auto first_node: {VC_ASC, VC_DESC, VC_RAND})
  for(auto second_node: {VC_NONE, VC_ASC, VC_DESC, VC_RAND})
  for(auto update_deg: {true, false}) {
    cover = nodeGreedy_general(g, first_node, second_node, update_deg);
    if(check_cover(g, cover)) {
      Info("is minimal")
      if(cover.size() < smallestVC.size()) smallestVC.swap(cover);
    }
    else {
      for(auto postproc: {VC_ASC, VC_DESC, VC_RAND})
      for(auto score: {VC_DEG, VC_NEIGH, VC_RNEIGH}) {
        post_cover = greedyPruning_general(g, cover, postproc, score);
        check_cover(g, post_cover, true);
        if(post_cover.size() < smallestVC.size()) smallestVC.swap(post_cover);
      }
    }
  }

  cover.clear();
  for(ul u=0; u<g.n; ++u) cover.push_back(u);
  Info("all nodes: cover size " << cover.size())
  for(auto postproc: {VC_ASC, VC_DESC, VC_RAND})
  for(auto score: {VC_DEG, VC_NEIGH, VC_RNEIGH}) {
    post_cover = greedyPruning_general(g, cover, postproc, score);
    check_cover(g, post_cover, true);
    if(post_cover.size() < smallestVC.size()) smallestVC.swap(post_cover);
  }

  Info("smallestVC: size "<< smallestVC.size())
  return smallestVC;
}



double certifVC(const Adjlist &g, const ul &higestLB, const vector<ul> &smallestVC) {
  double certified_quality = ((double) smallestVC.size()) / higestLB;
  Info("Certified quality: "<< certified_quality)
  return certified_quality;
}
