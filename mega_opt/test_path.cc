#include <cstdio>
#include <map>

#include "utility.h"

using namespace std;

int TestPath(const vector<int>& p1, const vector<int>& p2,
             map<int, pair<double, double> >& points) {
  if (p1.size() != points.size()) {
    printf("too few points");
    return 1;
  }
  set<pair<int, int> > bad_edges;
  double dist1 = 0, dist2 = 0;
  int previd = p1[0];
  set<int> visited;
  visited.insert(previd);
  for (int i = 1; i < p1.size(); i++) {
    int id = p1[i];
    visited.insert(id);
    if (points.count(id) == 0 || points.count(previd) == 0) {
      printf("unknown point\n");
      return 1;
    }
    dist1 += CalcDist(points[id], points[previd]);
    bad_edges.insert(make_pair(id, previd));
    bad_edges.insert(make_pair(previd, id));
    previd = id;
  }
  if (visited.size() != points.size()) {
    printf("didn't visit all points on path 1\n");
    return 1;
  }
  printf("path 1 length: %lf\n", dist1);

  if (p2.size() != points.size()) {
    printf("too few points");
    return 1;
  }
  previd = p2[0];

  visited.clear();
  visited.insert(previd);
  for (int i = 1; i < points.size(); i++) {
    int id = p2[i];
    visited.insert(id);
    if (points.count(id) == 0 || points.count(previd) == 0) {
      printf("unknown point\n");
      return 1;
    }
    dist2 += CalcDist(points[id], points[previd]);
    if (bad_edges.count(make_pair(id, previd))) {
      printf("double edge %d %d\n", id, previd);
    }
    previd = id;
  }
  if (visited.size() != points.size()) {
    printf("didn't visit all points on path 2\n");
    return 1;
  }
  printf("path 2 length: %lf\n", dist2);
  printf("score: %lf\n", max(dist1, dist2));
  printf("avg: %lf\n", (dist1+dist2)/2);
  return 0;
}

int main(int argc, char **argv) {
  FILE *fpath = fopen(argv[1], "r");
  FILE *fpoints = fopen(argv[2], "r");

  int id, x, y;
  map<int, pair<double, double> > points;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }

  vector<int> p1, p2;
  for (int i = 0; i < points.size(); i++) {
    int a, b; fscanf(fpath, "%d,%d", &a, &b);
    p1.push_back(a); p2.push_back(b);
  }
  return TestPath(p1, p2, points);
}
