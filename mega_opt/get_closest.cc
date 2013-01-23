#include <cstdio>
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char**argv) {
  FILE *fpoints = fopen(argv[1], "r");

  int id, x, y;
  vector<pair<long long, long long> > points(150000);
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }

  fprintf(stderr, "loaded\n");
  int n = 100;

  for (int i = 0; i < 150000; i++) {
    vector<pair<long long, int> > closest;
    for (int j = 0; j < 150000; j++) {
      if (i == j) continue;
      long long d = (points[i].first - points[j].first)*(points[i].first - points[j].first) +
                    (points[i].second - points[j].second)*(points[i].second - points[j].second);
      closest.push_back(make_pair(d, j));
      if (closest.size() > 2*n) {
        sort(closest.begin(), closest.end());
        closest.resize(n);
      }
    }
    sort(closest.begin(), closest.end());
    for (int j = 0; j < n; j++)
      printf("%d ", closest[j].second);
    printf("\n");
    fprintf(stderr, "%d\n", i);
  }
}
