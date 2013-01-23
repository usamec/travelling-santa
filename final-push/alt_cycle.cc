#include "alt_cycle.h"
#include <cassert>
#include <queue>
#include <tr1/unordered_map>
#include <map>

vector<int> Graph4::NormalizePath(vector<int>& p) {
  int mp = 0;
  for (int i = 0; i < p.size(); i++)
    if (p[i] < p[mp]) mp = i;
  vector<int> ret;
  for (int i = 0; i < p.size(); i++) {
    int ii = (i + mp) % p.size();
    ret.push_back(p[ii]);
  }
  sort(ret.begin(), ret.end());
  return ret;
}

void Graph4::Alternate(vector<int>& cycle) {
  vector<pair<int, int> > in, out;
  in.push_back(make_pair(cycle[0], cycle.back()));
  for (int i = 1; i < cycle.size(); i++) {
    if (i % 2 == 1) out.push_back(make_pair(cycle[i], cycle[i-1]));
    if (i % 2 == 0) in.push_back(make_pair(cycle[i], cycle[i-1]));
  }
  for (int i = 0; i < out.size(); i++) {
    int a = out[i].first, b = out[i].second;
    next_[a].erase(b);
    next_[b].erase(a);
  }
  for (int i = 0; i < in.size(); i++) {
    int a = in[i].first, b = in[i].second;
    next_[a].insert(b);
    next_[b].insert(a);
  }
  // check
  map<int, int> nd;
  for (int i = 0; i < next_.size(); i++)
    nd[next_[i].size()]++;
  for (auto it = nd.begin(); it != nd.end(); ++it)
    printf("nd %d %d\n", it->first, it->second);
  if (nd[3] != 4 || nd[4] != points_.size() - 4) {
    printf("bad nd\n");
    exit(1);
  }
}

vector<int> Graph4::FindAlternatingCycleBellman(int start, int limit,
    unordered_set<vector<int> >& banned_cycles) const {
  int bc = 0;
  for (int i = 0; i < next_.size(); i++) {
    int ii = (i + start) % next_.size();
    if (ii % 10000 == 0)
      printf("alter dfs %d %d %d\n", ii, bc, limit);

    
    vector<int> altDFS;
    unordered_set<pair<int, int> > cycle_edge_set_dfs;
    altDFS.push_back(ii);
    if (FindAlternatingCycleDFSExpand(altDFS, cycle_edge_set_dfs, 0, 300000, banned_cycles)) {
      printf("find dfs %d %d\n", altDFS.size(), ii);
      return altDFS;
    }
  }

  for (int i = 0; i < next_.size(); i++) {
    int ii = (i + start) % next_.size();
    if (ii % 1000 == 0)
      printf("alter %d %d %d\n", ii, bc, limit);

    
    unordered_map<pair<int, int>, int> prev;
    unordered_map<int, double> dd;
    dd[ii] = 0;
    for (int j = 0; j < limit; j++) {
      if (j % 2 == 0) {
        unordered_map<int, double> dd2;
        for (auto it = dd.begin(); it != dd.end(); ++it) {
          int a = it->first;
          unordered_set<pair<int, int> > prev2;
          int cur = a;
          for (int k = j-1; k >= 1 && k >= j-24; k--) {
            int pr = cur;
            cur = prev[make_pair(cur, k)];
            prev2.insert(make_pair(min(cur, pr), max(cur, pr)));
          }
          double d = it->second;
          for (auto it2 = next_[a].begin(); it2 != next_[a].end(); ++it2) {
            int b = *it2;
            if (prev2.count(make_pair(min(a,b), max(a,b)))) continue;
            double d2 = d - CalcDist(points_[a], points_[b]);
            if (dd2.count(b) == 0 || dd2[b] > d2) {
              dd2[b] = d2;
              prev[make_pair(b, j)] = a;
            }
          }
        }
        dd = dd2;
      } else {
        unordered_map<int, double> dd2;
        for (auto it = dd.begin(); it != dd.end(); ++it) {
          int a = it->first;
          unordered_set<pair<int, int> > prev2;
          int cur = a;
          for (int k = j-1; k >= 1 && k >= j-24; k--) {
            int pr = cur;
            cur = prev[make_pair(cur, k)];
            prev2.insert(make_pair(min(cur, pr), max(cur, pr)));
          }
          double d = it->second;

          // try to close cycle
          if (a != ii && next_[a].count(ii) == 0) {
            double dc = d + CalcDist(points_[ii], points_[a]);
            if (dc < 0) {
              vector<int> cycle;
              cycle.reserve(j+1);
              unordered_set<pair<int, int> > edge_set;
              int cur = a;
              bool bad = false;
              for (int k = j-1; k >= 0; k--) {
                cycle.push_back(cur);
                cur = prev[make_pair(cur, k)];
              }
              cycle.push_back(cur);
              edge_set.insert(make_pair(min(cycle[0], cycle.back()),
                                        max(cycle[0], cycle.back())));
              for (int k = 1; k < cycle.size(); k++) {
                if (edge_set.count(make_pair(min(cycle[k], cycle[k-1]),
                                             max(cycle[k], cycle[k-1])))) {
                  bad = true;
                  break;
                }
                edge_set.insert(make_pair(min(cycle[k], cycle[k-1]),
                                           max(cycle[k], cycle[k-1]))); 
              }

              if (bad) {
/*                printf("bad %lf\n", dc);
                for (int i = 0; i < cycle.size(); i++)
                  printf("%d ", cycle[i]);
                printf("\n");*/
                bc++;
              }
              if (cur == ii && bad == false && banned_cycles.count(NormalizePath(cycle)) == 0) {
                printf("improv %lf %d bc %d\n", dc, cycle.size(), bc);
                reverse(cycle.begin(), cycle.end());
                return cycle;
              }
            }
          }

          for (int k = 0; k < closest_[a].size(); k++) {
            int b = closest_[a][k];
            if (prev2.count(make_pair(min(a,b), max(a,b)))) continue;
            if (next_[a].count(b)) continue;
            double d2 = d + CalcDist(points_[a], points_[b]);
            if (d2 < 0 && (dd2.count(b) == 0 || dd2[b] > d2)) {
              dd2[b] = d2;
              prev[make_pair(b, j)] = a;
            }
          }
        }
        dd = dd2;
      }
    }
  }
  return vector<int>();
}

vector<int> Graph4::FindAlternatingCycleDFS(int start, int dfs_limit,
          unordered_set<vector<int> >& banned_cycles) const {
  for (int i = 0; i < next_.size(); i++) {
    int ii = (i + start) % next_.size();
    if (ii % 1000 == 0)
      printf("alter %d\n", ii);
    vector<int> cycle;
    unordered_set<pair<int, int> > cycle_edge_set;
    cycle.push_back(ii);
    if (FindAlternatingCycleDFSExpand(cycle, cycle_edge_set, 0, dfs_limit, banned_cycles))
      return cycle;
  }
  return vector<int>();
}

bool Graph4::FindAlternatingCycleDFSExpand(
    vector<int>& cycle, unordered_set<pair<int, int> >& cycle_edge_set,
    double cur_cost, int dfs_limit,
    unordered_set<vector<int> >& banned_cycles) const {
  if (cycle.size() == dfs_limit) return false;
  if (cycle.size() > 1)
    if (cycle.back() < cycle[0]) return false;
  double cur_cost_backup = cur_cost;
  if (cycle.size() % 2 == 0) {
    // try to close_cycle 
    if (banned_cycles.count(NormalizePath(cycle)) == 0 &&
        next_[cycle.back()].count(cycle[0]) == 0 && cycle.back() != cycle[0]) {
      if (cycle_edge_set.count(make_pair(
              min(cycle[0], cycle.back()), max(cycle[0], cycle.back())))==0) {
        cur_cost += CalcDist(points_[cycle.back()], points_[cycle[0]]);
        if (cur_cost < 0) {
          printf("find improv %lf\n", cur_cost);
          return true;
        }
        cur_cost = cur_cost_backup;
      }
    }
    // insert edge from good_next_
    int cur = cycle.back();
    for (auto it = good_next_[cur].begin(); it != good_next_[cur].end(); ++it) {
      int next = *it;
      if (next_[cur].count(next) != 0) continue;
      if (cycle_edge_set.count(make_pair(min(next, cur), max(next, cur))) != 0) continue;
      cycle.push_back(next);
      cycle_edge_set.insert(make_pair(min(next, cur), max(next, cur)));
      cur_cost += CalcDist(points_[cur], points_[next]);
      if (FindAlternatingCycleDFSExpand(cycle, cycle_edge_set, cur_cost, dfs_limit, banned_cycles))
        return true;
      cur_cost = cur_cost_backup;
      cycle_edge_set.erase(make_pair(min(next, cur), max(next, cur)));
      cycle.pop_back();
    }
  } else {
    // insert (erase) edge from next_
    int cur = cycle.back();
    for (auto it = next_[cur].begin(); it != next_[cur].end(); ++it) {
      int next = *it;
      if (good_next_[cur].count(next) != 0) continue;
      if (cycle_edge_set.count(make_pair(min(next, cur), max(next, cur))) != 0) continue;
      cycle.push_back(next);
      cycle_edge_set.insert(make_pair(min(next, cur), max(next, cur)));
      cur_cost -= CalcDist(points_[cur], points_[next]);
      if (FindAlternatingCycleDFSExpand(cycle, cycle_edge_set, cur_cost, dfs_limit, banned_cycles))
        return true;
      cur_cost = cur_cost_backup;
      cycle_edge_set.erase(make_pair(min(next, cur), max(next, cur)));
      cycle.pop_back();
    }
  }
  return false;
}

vector<int> Graph4::GetNearPoints(vector<int>&start, int dist_lim) const {
  unordered_map<int, int> dist;
  queue<int> fr;
  vector<int> ret;
  for (int i = 0; i < start.size(); i++) {
    if (dist.count(start[i])) continue;
    dist[start[i]] = 0;
    fr.push(start[i]);
  }
  while (!fr.empty()) {
    int x = fr.front(); fr.pop();
    int d = dist[x];
    ret.push_back(x);
    if (d == dist_lim) continue;
    for (auto it = next_[x].begin(); it != next_[x].end(); ++it) {
      int nn = *it;
      if (dist.count(nn) != 0) continue;
      dist[nn] = d + 1;
      fr.push(nn);
    }
  }
  return ret;
}

/*int main(int argc, char **argv) {
  FILE *fpoints = fopen(argv[2], "r");

  int id, x, y;
  vector<pair<double, double> > points;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points.push_back(make_pair(x, y));
  }
  FILE *fpath = fopen("xo", "r");

  vector<int> p1, p2, ip1, ip2;
  for (int i = 0; i < points.size(); i++) {
    int a, b; fscanf(fpath, "%d,%d", &a, &b);
    p1.push_back(a); p2.push_back(b);
  }

  FILE *fclosest = fopen("closest.dat", "r");
  FILE *fedges = fopen("edge.dat", "r");
  vector<vector<int> > closest, edge4(points.size());

  for (int i = 0; i < points.size(); i++) {
    vector<int> x;
    unordered_set<int> ss;
    for (int j = 0; j < 4; j++) {
      int ii; int a; fscanf(fedges, "%d %d", &ii, &a);
      assert(ii == i);
      x.push_back(a);
      ss.insert(a);
      edge4[i].push_back(a);
    }
    for (int j = 0; j < 100; j++) {
      int a; fscanf(fclosest, "%d", &a);
      if (x.size() < 7 && ss.count(a) == 0)
        x.push_back(a);
    }
    closest.push_back(x);
  }

  Graph4 g4(p1, p2, edge4, closest, points);

  int dfs_limit = 20;
  vector<int> cycle = g4.FindAlternatingCycleDFS(dfs_limit);
  vector<int> near = g4.GetNearPoints(cycle, 2);
  for (int i = 0; i < cycle.size(); i++)
    printf("%d ", cycle[i]);
  printf("\n");
  for (int i = 0; i < near.size(); i++)
    printf("%d ", near[i]);
  printf("\n");
}*/
