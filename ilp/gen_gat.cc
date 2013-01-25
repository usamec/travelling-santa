#include <cstdio>
#include <map>
#include <queue>

#include "utility.h"
#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

using namespace std;
boost::mutex mutex;

double bd;
double bm;
map<int, pair<double, double> > points;

void Optimize() {
  queue<pair<int, int> > buffer;
  int id = rand()%1000;
  char inname[20];
  char outname[20];
  char logname[20];
  sprintf(inname, "gr%d.tsp", id);
  sprintf(outname, "out%d.tt", id);
  sprintf(logname, "x%d.log", id);

  while(true) {
    mutex.lock();
    FILE *fpath = fopen("xo", "r");

    vector<int> p1, p2, ip1, ip2;
    for (int i = 0; i < points.size(); i++) {
      int a, b; fscanf(fpath, "%d,%d", &a, &b);
      p1.push_back(a); p2.push_back(b);
    }
    fclose(fpath);
    mutex.unlock();
    ip1.resize(p1.size()); ip2.resize(p2.size());

    double d1 = 0, d2 = 0;
    for (int i = 1; i < p1.size(); i++) {
      d1 += CalcDist(points[p1[i]], points[p1[i-1]]);
      d2 += CalcDist(points[p2[i]], points[p2[i-1]]);
    }
    for (int i = 0; i < p1.size(); i++) {
      ip1[p1[i]] = i;
      ip2[p2[i]] = i;
    }

    map<int, int> points_in_r;
    vector<int> points_in;
    int mi = 100, ma = 150;
    FILE *fin = fopen(inname, "w");

    //printf("loaded\n");
    double a, c;
    if (buffer.empty()) {
      a = rand() % 20000;
      c = rand() % 20000;
    } else {
      a = buffer.front().first;
      c = buffer.front().second;
      buffer.pop();
    }
    double bmi = 0;
    double bma = 20000;
    double l, r, u, dow;
    while (bma - bmi > 0.001) {
      double bm = (bma + bmi) / 2;
      double b = a + bm;
      double d = c + bm;
      l = min(a,b), r = max(a, b);
      u = min(c,d), dow = max(c, d);
      points_in.clear();
      for (int i = 0; i < p1.size(); i++) {
        if (points[p1[i]].first < r && points[p1[i]].second < dow && 
            points[p1[i]].first > l && points[p1[i]].second > u) {
          points_in.push_back(p1[i]);
        }
        if (points_in.size() > ma) break;
      }
      if (points_in.size() >= mi && points_in.size() <= ma)
        break;
      if (points_in.size() < mi) bmi = bm;
      else bma = bm;
    }
    if (points_in.size() < mi || points_in.size() > ma+3) continue;
/*    while (points_in.size() < mi || points_in.size() > ma) {
  //    printf("ps %d\n", points_in.size());
      points_in.clear();
      points_in_r.clear();
      bool last_in = false;
      int a = rand() % 20000;
      int b = rand() % 1000;
      int c = rand() % 20000;
      int d = rand() % 20000;
      d = c + b;
      b = a + b;
      int l = min(a,b), r = max(a, b);
      int u = min(c,d), dow = max(c, d);
      for (int i = 0; i < p1.size(); i++) {
        if (points[p1[i]].first < r && points[p1[i]].second < dow && 
            points[p1[i]].first > l && points[p1[i]].second > u) {
          points_in.push_back(p1[i]);
        }
        if (points_in.size() > ma) break;
      }
//      fprintf(stderr, "pi %d\n", points_in.size());
    }*/
    for (int i = 0; i < points_in.size(); i++)
      points_in_r[points_in[i]] = i;
    int N = points_in.size();
    //printf("N %d %lf %lf %lf %lf\n", N, l, r, u, dow);
    fprintf(fin, "%lf\n", (d1+d2)/2.0); //max(d1, d2));
    fprintf(fin, "%d\n", N);
    for (int i = 0; i < points_in.size(); i++) {
      for (int j = 0; j < points_in.size(); j++) {
        fprintf(fin, "%lf ", CalcDist(points[points_in[i]], points[points_in[j]]));
      }
      fprintf(fin, "\n");
    }

    vector<vector<int> > nei1(p1.size()), nei2(p2.size());
    set<pair<int, int> > fixed1, fixed2;
    int last = -1;
    int first = -1;
    bool last_in = true;
    int last_pos = -1;
    vector<pair<int, int> > fixed;
    vector<double> fixed_l;
    double acc = 0;
    double first_acc = 0;
    for (int i = 0; i < p1.size(); i++) {
      if (i > 0) acc += CalcDist(points[p1[i-1]], points[p1[i]]);
      if (first == -1 && points_in_r.count(p1[i])) {
        first = p1[i];
        first_acc = acc;
      }

      if (last_in == false && last != -1 && points_in_r.count(p1[i])) {
        fixed.push_back(make_pair(points_in_r[last], points_in_r[p1[i]]));
        fixed_l.push_back(acc);
        fixed1.insert(make_pair(last, p1[i]));
        fixed1.insert(make_pair(p1[i], last));
//        fprintf(stderr, "fix %d %d\n", i, last_pos);
      }

      bool in = points_in_r.count(p1[i]);
      if (i > 0){
        if (!last_in || !in) {
          nei1[p1[i-1]].push_back(p1[i]);
          nei1[p1[i]].push_back(p1[i-1]);
        }
      }

      last_in = false;
      if (points_in_r.count(p1[i])) {
        last_in = true;
        last = p1[i];
        last_pos = i;
        acc = 0;
      }
    }
    fixed.push_back(make_pair(points_in_r[last], points_in_r[first]));
    fixed1.insert(make_pair(last, first));
    fixed1.insert(make_pair(first, last));
    fixed_l.push_back(first_acc + acc);
    fprintf(fin, "%d\n", fixed.size());
    fprintf(stderr, "%d %d\n", N, fixed.size());
    for (int i = 0; i < fixed.size(); i++) {
      fprintf(fin, "%d %d %lf\n", fixed[i].first, fixed[i].second, fixed_l[i]);
    }
    last = -1;
    first = -1;
    last_in = true;
    fixed.clear();
    fixed_l.clear();
    acc = 0;
    first_acc = 0;
    for (int i = 0; i < p2.size(); i++) {
      if (i > 0) acc += CalcDist(points[p2[i-1]], points[p2[i]]);
      if (first == -1 && points_in_r.count(p2[i])) {
        first = p2[i];
        first_acc = acc;
      }

      if (last_in == false && last != -1 && points_in_r.count(p2[i])) {
        fixed.push_back(make_pair(points_in_r[last], points_in_r[p2[i]]));
        fixed_l.push_back(acc);
        fixed2.insert(make_pair(last, p2[i]));
        fixed2.insert(make_pair(p2[i], last));
//        fprintf(stderr, "fix %d %d\n", i, last_pos);
      }
      bool in = points_in_r.count(p2[i]);
      if (i > 0) {
        if (!last_in || !in) {
          nei2[p2[i-1]].push_back(p2[i]);
          nei2[p2[i]].push_back(p2[i-1]);
        }
      }

      last_in = false;
      if (points_in_r.count(p2[i])) {
        last_in = true;
        last = p2[i];
        last_pos = i;
        acc = 0;
      }
    }
    fixed.push_back(make_pair(points_in_r[last], points_in_r[first]));
    fixed2.insert(make_pair(last, first));
    fixed2.insert(make_pair(first, last));
    fixed_l.push_back(first_acc + acc);
    fprintf(fin, "%d\n", fixed.size());
    fprintf(stderr, "%d %d\n", N, fixed.size());
    for (int i = 0; i < fixed.size(); i++) {
      fprintf(fin, "%d %d %lf\n", fixed[i].first, fixed[i].second, fixed_l[i]);
    }
    fclose(fin);
    char command[100];
    // presolving
    sprintf(command, "./presolve.py %s %s >%s", inname, outname, logname);
    system(command);
    FILE *fout = fopen(outname, "r");
    char ttt[100];
    fscanf(fout, "%s", ttt);
    printf("** %s\n", ttt);
    fclose(fout);
    if (ttt[0] == 'w') {
        continue; // was_ok
    }
    
    sprintf(command, "./usama.py %s %s >%s", inname, outname, logname);
    system(command);
    
    fout = fopen(outname, "r");
    fscanf(fout, "%s", ttt);
    printf("** %s\n", ttt);
    assert(ttt[0] == 'n'); // new_solution
    printf("good N %d %lf %lf %lf %lf\n", N, l, r, u, dow);

    double smx[] = {0.5, 0.5, -0.5, -0.5};
    double smy[] = {0.5, -0.5, 0.5, -0.5};
    for (int i = 0; i < 4; i++) {
      buffer.push(make_pair(l + smx[i]*(r-l), u + smy[i]*(dow-u)));
    }

    vector<int> pp1, pp2;
    for (int i = 0; i < N; i++) {
      int a; fscanf(fout, "%d", &a);
      pp1.push_back(points_in[a]);
    }
    for (int i = 0; i < N; i++) {
      int a; fscanf(fout, "%d", &a);
      pp2.push_back(points_in[a]);
    }

    for (int i = 0; i < pp1.size(); i++) {
      int nn = (i+1)%N;
      if (fixed1.count(make_pair(pp1[i], pp1[nn]))==0) {
        nei1[pp1[i]].push_back(pp1[nn]);
        nei1[pp1[nn]].push_back(pp1[i]);
      }
      if (fixed2.count(make_pair(pp2[i], pp2[nn]))==0) {
        nei2[pp2[i]].push_back(pp2[nn]);
        nei2[pp2[nn]].push_back(pp2[i]);
      }
    }

    vector<bool> vis(p1.size(), false);
    vector<int> np1, np2;

    int cur = p1[0];
    while (true) {
      np1.push_back(cur);
      int nn = -1;
      vis[cur] = true;
      for (int i = 0; i < nei1[cur].size(); i++) {
        if (vis[nei1[cur][i]]) continue;
        nn = nei1[cur][i];
        break;
      }
      if (nn == -1) break;
      cur = nn;
    }
    if (np1.size() != p1.size()) {
      printf("fuck 1\n"); continue;
    }
    vis = vector<bool>(p2.size(), false);
    cur = p2[0];
    while (true) {
      np2.push_back(cur);
      int nn = -1;
      vis[cur] = true;
      for (int i = 0; i < nei2[cur].size(); i++) {
        if (vis[nei2[cur][i]]) continue;
        nn = nei2[cur][i];
        break;
      }
      if (nn == -1) break;
      cur = nn;
    }
    if (np2.size() != p2.size()) {
      printf("fuck 2 %d %d\n", np2.size(), p2.size());
      map<int, int> nd;
      for (int i = 0; i < p2.size(); i++)
        nd[nei2[i].size()]++;
      for (auto it = nd.begin(); it != nd.end(); ++it)
        printf("%d %d\n", it->first, it->second);
      continue;
    }
    double nd1 = 0, nd2 = 0;
    for (int i = 1; i < p1.size(); i++) {
      nd1 += CalcDist(points[np1[i]], points[np1[i-1]]);
      nd2 += CalcDist(points[np2[i]], points[np2[i-1]]);
    }
    if ((nd1+ nd2)/2 > (d1+ d2)/2+1) {
      printf("fuck 3 %lf %lf %lf %lf\n", nd1, nd2, d1, d2); continue;
    }
    printf("old %lf %lf %lf %lf\n", d1, d2, max(d1, d2), (d1+d2)/2);
    printf("new %lf %lf %lf %lf\n", nd1, nd2, max(nd1, nd2), (nd1+nd2)/2);
    mutex.lock();
    if ((nd1+ nd2)/2 < bd) {
      FILE *fo = fopen("xo", "w");
      for (int i = 0; i < p1.size(); i++) {
        fprintf(fo, "%d,%d\n", np1[i], np2[i]);
      }
      fclose(fo);
      bd = (nd1+ nd2)/2;
      printf("new bd %lf\n", bd);
    }
    if (max(nd1, nd2) < bm) {
      bm = max(nd1, nd2);
      FILE *fo = fopen("xom", "w");
      for (int i = 0; i < p1.size(); i++) {
        fprintf(fo, "%d,%d\n", np1[i], np2[i]);
      }
      fclose(fo);
    }
    mutex.unlock();
  }
}

int main(int argc, char **argv) {
  srand(atoi(argv[3]));
  FILE *fpoints = fopen(argv[2], "r");

  int id, x, y;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }
  FILE *fpath = fopen("xo", "r");

  vector<int> p1, p2, ip1, ip2;
  for (int i = 0; i < points.size(); i++) {
    int a, b; fscanf(fpath, "%d,%d", &a, &b);
    p1.push_back(a); p2.push_back(b);
  }

  double d1 = 0, d2 = 0;
  for (int i = 1; i < p1.size(); i++) {
    d1 += CalcDist(points[p1[i]], points[p1[i-1]]);
    d2 += CalcDist(points[p2[i]], points[p2[i-1]]);
  }

  bd = (d1+ d2)/2;
  bm = max(d1, d2);
  boost::thread th1(Optimize);
  boost::thread th2(Optimize);
  boost::thread th3(Optimize);
  boost::thread th4(Optimize);
  boost::thread th5(Optimize);
  boost::thread th6(Optimize);
  boost::thread th7(Optimize);
  th1.join();
  return 0;
}
