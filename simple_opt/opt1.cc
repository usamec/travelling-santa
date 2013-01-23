#include <cstdio>
#include <unordered_map>
#include <map>

#include "utility.h"

#include "opts.h"

using namespace std;


int main(int argc, char **argv) {
  srand(time(NULL));

  // Start with this opt for couple of times
  FILE *fpath = fopen(argv[1], "r");
  FILE *fpoints = fopen(argv[2], "r");

  int id, x, y;
  unordered_map<int, pair<double, double> > points;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }

  vector<int> p1, p2;
  for (int i = 0; i < points.size(); i++) {
    int a, b; fscanf(fpath, "%d,%d", &a, &b);
    p1.push_back(a); p2.push_back(b);
  }
  fclose(fpath);
  unordered_set<pair<int, int> > used_edges;
  double d1 = 0, d2 = 0;
  for (int i = 1; i < p1.size(); i++) {
    used_edges.insert(make_pair(p1[i-1], p1[i]));
    used_edges.insert(make_pair(p1[i], p1[i-1]));
    d1 += CalcDist(points[p1[i-1]], points[p1[i]]);
    used_edges.insert(make_pair(p2[i-1], p2[i]));
    used_edges.insert(make_pair(p2[i], p2[i-1]));
    d2 += CalcDist(points[p2[i-1]], points[p2[i]]);
  }
  printf("d1 %lf d2 %lf\n", d1, d2);

  for (int i = 0; i < 10; i++) {
    if (d1 > d2) {
      d1 = OptPath(p1, d1, used_edges, points);
      printf("d1 %lf d2 %lf\n", d1, d2);
      d2 = OptPath(p2, d2, used_edges, points);
    } else {
      d2 = OptPath(p2, d2, used_edges, points);
      printf("d1 %lf d2 %lf\n", d1, d2);
      d1 = OptPath(p1, d1, used_edges, points);
    }
    printf("d1 %lf d2 %lf\n", d1, d2);
  }
  FILE *fo = fopen(argv[1], "w");
  for (int i = 0; i < p1.size(); i++) {
    fprintf(fo, "%d,%d\n", p1[i], p2[i]);
  }
  fclose(fo);
  while (true) {
    FILE *fpath = fopen(argv[1], "r");
    FILE *fpoints = fopen(argv[2], "r");

    int id, x, y;
    unordered_map<int, pair<double, double> > points;
    while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
      points[id] = make_pair(x, y);
    }

    vector<int> p1, p2;
    for (int i = 0; i < points.size(); i++) {
      int a, b; fscanf(fpath, "%d,%d", &a, &b);
      p1.push_back(a); p2.push_back(b);
    }
    fclose(fpath);
    unordered_set<pair<int, int> > used_edges;
    double d1 = 0, d2 = 0;
    for (int i = 1; i < p1.size(); i++) {
      used_edges.insert(make_pair(p1[i-1], p1[i]));
      used_edges.insert(make_pair(p1[i], p1[i-1]));
      d1 += CalcDist(points[p1[i-1]], points[p1[i]]);
      used_edges.insert(make_pair(p2[i-1], p2[i]));
      used_edges.insert(make_pair(p2[i], p2[i-1]));
      d2 += CalcDist(points[p2[i-1]], points[p2[i]]);
    }
    printf("d1 %lf d2 %lf\n", d1, d2);
    double best_score = max(d1, d2);
    double last_score = best_score;
    bool rozdrb = false;

    int gi = 0;
    last_score = max(d1, d2);
    int rozdrb_count = 0;
    while (rozdrb_count < 4) {
      gi++;

      for (int iiii = 0; iiii < 2; iiii++) {
        sector_div = 500 + rand()%100;
        printf("sector div %lf\n", sector_div);
        if (rozdrb) {
          printf("unstuck\n");
          kill_coef = 2;
          kill_limit = 20000;
          if (d1 > d2) {
            LocalOptPathRHarder(p1, d1, used_edges, points, p2, d2, 1);
          } else {
            LocalOptPathRHarder(p2, d2, used_edges, points, p1, d1, 1);
          }
          kill_coef = 1.5;
          for (int i = 0; i < 10; i++) {
            printf("d1 %lf d2 %lf\n", d1, d2);
            if (d1 > d2) {
              KillCrossingsHarder(p1, d1, used_edges, points, p2, d2);
            } else {
              KillCrossingsHarder(p2, d2, used_edges, points, p1, d1);
            }
          }
          kill_coef = 1.5;
          kill_limit = 2000;
          double lm = max(d1, d2);
          double lim = 900;
          while (fabs(d1 - d2) > lim || lim < 1000) {
            printf("d1 %lf d2 %lf %lf\n", d1, d2, lim);
            if (d1 > d2) {
              KillCrossingsHarder(p1, d1, used_edges, points, p2, d2);
            } else {
              KillCrossingsHarder(p2, d2, used_edges, points, p1, d1);
            }
            if (fabs(max(d1, d2) - lm) < 5)
              break;
            lm = max(d1, d2);
            lim += 10; kill_limit -= 10;
          }
        }
        printf("d1 %lf d2 %lf\n", d1, d2);
        if (d1 > d2) {
          d1 = OptPath(p1, d1, used_edges, points);
        } else {
          d2 = OptPath(p2, d2, used_edges, points);
        } 
        printf("d1 %lf d2 %lf\n", d1, d2);
        for (int i = 0; i < 2; i++) {
          if (d2 > d1) {
            printf("d1 %lf d2 %lf\n", d1, d2);
            d2 = KillCrossings(p2, d2, used_edges, points);
          } else {
            printf("d1 %lf d2 %lf\n", d1, d2);
            d1 = KillCrossings(p1, d1, used_edges, points);
          }
        }
        printf("d1 %lf d2 %lf\n", d1, d2);
        if (d1 > d2) 
          d1 = LocalOptPath(p1, d1, used_edges, points);
        else
          d2 = LocalOptPath(p2, d2, used_edges, points);
        printf("d1 %lf d2 %lf\n", d1, d2);
        if (d1 > d2) {
          d1 = LocalPermsC(p1, d1, used_edges, points, 10);
        } else {
          d2 = LocalPermsC(p2, d2, used_edges, points, 10);
        }
        printf("d1 %lf d2 %lf\n", d1, d2);
        if (iiii == 0) { 
          if (d1 > d2) {
            d1 = LocalPermsD(p1, d1, used_edges, points, 12);
          } else {
            d2 = LocalPermsD(p2, d2, used_edges, points, 12);
          }
        }
  //      LocalOptPaths(p1, d1, p2, d2, used_edges, points);
        printf("d1 %lf d2 %lf\n", d1, d2);
        double dd1 = 0, dd2 = 0;
        for (int i = 1; i < p1.size(); i++) {
          dd1 += CalcDist(points[p1[i-1]], points[p1[i]]);
          dd2 += CalcDist(points[p2[i-1]], points[p2[i]]);
        }
        d1 = dd1; d2 = dd2;
        double score = max(d1, d2);
        if (last_score - score < 10 && rozdrb == false) {
          rozdrb_count++;
          rozdrb = true;
        } else
          rozdrb = false;
        last_score = score;
        printf("dd1 %lf dd2 %lf cur %lf best %lf\n", dd1, dd2, score, best_score);
        FILE *flog = fopen("opt.log", "a");
        fprintf(flog, "%lld %lf %lf %lf %lf\n", time(NULL)-1355852120, best_score, score, dd1, dd2);
        fclose(flog);
        if (score < best_score) {
          best_score = score;
          FILE *fo = fopen("best.dat", "w");
          for (int i = 0; i < p1.size(); i++) {
            fprintf(fo, "%d,%d\n", p1[i], p2[i]);
          }
          fclose(fo);
        }
        FILE *fo = fopen(argv[1], "w");
        for (int i = 0; i < p1.size(); i++) {
          fprintf(fo, "%d,%d\n", p1[i], p2[i]);
        }
        fclose(fo);
        reverse(p1.begin(), p1.end());
        reverse(p2.begin(), p2.end());
      }
    }
  }
}
