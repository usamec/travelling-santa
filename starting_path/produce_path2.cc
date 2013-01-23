#include "utility.h"
#include <cstdio>
#include <map>
#include <set>

int test_count = 200;
int N = 150000;

vector<int> ProducePath(map<int, pair<double, double> >& points,
    set<pair<int, int> >& banned_edges) {
  while (true) {
    vector<int> ids;
    for (auto it = points.begin(); it != points.end(); ++it) {
      ids.push_back(it->first);
    }
    int f = rand()%ids.size();
    vector<int> path;
    path.push_back(ids[f]);
    RemoveFromVector(ids, f);

    bool bad = false;
    while (ids.size() > 0) {
      int best_id = -1;
      double best_dist = 0;
      for (int i = 0; i < test_count && i < ids.size(); i++) {
        int p;
        if (test_count < ids.size())
          p = rand()%ids.size();
        else
          p = i;
        double dist = CalcDist(
            points[path.back()], points[ids[p]]);
        if (banned_edges.count(make_pair(path.back(), ids[p])))
          continue;
        if (best_id == -1 || dist < best_dist) {
          best_id = p;
          best_dist = dist;
        }
      }
      if (best_id == -1) {
        bad = true;
        break;
      }
      path.push_back(ids[best_id]);
      RemoveFromVector(ids, best_id);
      if (ids.size() % 10000 == 0)
        fprintf(stderr, "ids %d\n", ids.size());
    }
    if (bad) continue;
    for (int i = 1; i < path.size(); i++) {
      banned_edges.insert(make_pair(path[i], path[i-1]));
      banned_edges.insert(make_pair(path[i-1], path[i]));
    }
    return path;
  }
}

int main(int argc, char **argv) {
  srand(time(NULL));
  FILE *fpoints = fopen(argv[1], "r");

  int id, x, y;
  map<int, pair<double, double> > points;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }
  set<pair<int, int> > banned_edges;
  vector<int> p1 = ProducePath(points, banned_edges); 
  vector<int> p2 = ProducePath(points, banned_edges); 

  for (int i = 0; i < p1.size(); i++) {
    printf("%d,%d\n", p1[i], p2[i]);
  }
}
