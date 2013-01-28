#include <stdio.h>
#include <cmath>
#include <cassert>
#include <utility>
#include <set>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <map>

using namespace std;

#define __assert(cond, ...) assert(cond || !fprintf(stderr, __VA_ARGS__))
#define FOR(q,n) for(int q=0; q<n; q++)
#define FOREACH(it, t) for(__typeof(t.begin()) it=t.begin(); it!=t.end(); ++it)

#define mp make_pair
#define fi first
#define se second
#define pb push_back


bool file_exists(const char *filename) {
  struct stat buffer;
  return stat(filename, &buffer) == 0;
}

typedef pair<int, int> PII;
typedef pair<double, double> PDD;
#define N 150000

class Points {
public:
  double posx[N];
  double posy[N];
  double adjust[N];

  double dist_e2(PDD a, PDD b) {
    double dx = a.fi - b.fi;
    double dy = a.se - b.se;
    return sqrt(dx*dx + dy*dy + 1e-12);
    
  }

  double dist_e2(int a, int b) {
    double dx = posx[a] - posx[b];
    double dy = posy[a] - posy[b];
    return sqrt(dx*dx + dy*dy + 1e-12);
  }

  double dist(int a, int b) {
    return dist_e2(a, b) + adjust[a] + adjust[b];
  }

  static const int BUCKETS_X = 20;
  static const int BUCKETS_Y = 20;
  static const int MAX_X = 20000;
  static const int MAX_Y = 20000;
  static const int STEP_X = (MAX_X / BUCKETS_X + 1);
  static const int STEP_Y = (MAX_Y / BUCKETS_Y + 1);

  vector<int> buckets[BUCKETS_X][BUCKETS_Y];
  double min_adjust[BUCKETS_X][BUCKETS_Y];

  void load() {
    FILE* f= fopen("santa_cities.csv", "r");
    assert(f != NULL);

    FOR(q, N) {
      int id,x,y;
      fscanf(f, "%d,%d,%d", &id, &x, &y);
    //printf("%d %d %d\n", id, x, y);
      posx[q] = x;
      posy[q] = y;
      assert(id == q);
    }
    fill_buckets();
  }

  PII get_bucket_of_point(int q) {
      int bx = posx[q] / STEP_X;
      int by = posy[q] / STEP_Y;
      assert (0 <= bx && bx < BUCKETS_X);
      assert (0 <= by && by < BUCKETS_Y);
      return mp(bx, by);
  }

  void fill_buckets() {
    FOR(q, N) {
      PII b = get_bucket_of_point(q);
      buckets[b.fi][b.se].pb(q);
    }
  }

  double min_distance_point_bucket(int q, PII b) {
    int tmpx = posx[q];
    int tmpy = posy[q];
    if (tmpx < b.fi * STEP_X) tmpx = b.fi * STEP_X;
    if (tmpx > (b.fi + 1) * STEP_X) tmpx = (b.fi + 1) * STEP_X;
    
    if (tmpy < b.se * STEP_Y) tmpy = b.se * STEP_Y;
    if (tmpy > (b.se + 1) * STEP_Y) tmpy = (b.se + 1) * STEP_Y;

    double d = dist_e2(mp(posx[q], posy[q]), mp(tmpx, tmpy));
    d += adjust[q];
    d += min_adjust[b.fi][b.se];
    return d;
  }


  void update_minadjust() {
    FOR(q, BUCKETS_X)
      FOR(w, BUCKETS_Y)
        min_adjust[q][w] = min_adjust_in_bucket(mp(q, w));
  }
private:
  double min_adjust_in_bucket(PII b) {
    double m = 1e30;
    FOREACH(it, buckets[b.fi][b.se])
      if (adjust[*it] < m) m = adjust[*it];
    return m;
  }

};

Points points;


typedef pair<double, int> PDI;

void load_adjust() {
  FILE *f = fopen("adjust.dat", "r");
  FOR(q, N) fscanf(f, "%lf", &points.adjust[q]);
  fclose(f);
}

void normalize() {
  double sum = 0;
  FOR(q, N)
    sum += points.adjust[q];
  double delta = -sum / N;
  if (abs(delta) > 1e-6) printf("NEEDED ADJUST %lf!!!!!!!!\n", delta);
  FOR(q, N)
    points.adjust[q] += delta;
}

// **************
const int SOLVE_SIZE = 120;


int no_processed  = 0;

vector<double> merge_results(const vector<double>& res1, const vector<double>& res2) {
  assert(res1.size() == 5);
  assert(res2.size() == 5);
  vector<double> res;
  FOR(q, 5) {
    double best = 1e30;
    FOR(w, q + 1) {
      double t = res1[w] + res2[q-w];
      if (t < best) best = t;
    }
    res.push_back(best);
  }
  return res;
}

vector<double> recursive_lowerbound(const vector<int>& pts,
    double xlo, double xhi, double ylo, double yhi) {
  FOR(q, (int) pts.size()) {
    int v = pts[q];
    assert(points.posx[v] >= xlo);
    assert(points.posx[v] <= xhi);
    assert(points.posy[v] >= ylo);
    assert(points.posy[v] <= yhi);
  }

  if ((int) pts.size() > SOLVE_SIZE) {
    double xmid = (xlo + xhi) / 2;
    double ymid = (ylo + yhi) / 2;
    vector<int> tmp[2][2];
    FOR(q, (int) pts.size()) {
      int v = pts[q];
      tmp[points.posx[v] > xmid][points.posy[v] > ymid].push_back(v);
    }

    vector<double> res(5,0);
    vector<double> t;
    t = recursive_lowerbound(tmp[0][0], xlo, xmid, ylo, ymid);
    res = merge_results(res, t);
    
    t = recursive_lowerbound(tmp[0][1], xlo, xmid, ymid, yhi);
    res = merge_results(res, t);
    
    t = recursive_lowerbound(tmp[1][0], xmid, xhi, ylo, ymid);
    res = merge_results(res, t);
    
    t = recursive_lowerbound(tmp[1][1], xmid, xhi, ymid, yhi);
    res = merge_results(res, t);

    //printf("**************************************** %lf %lf %lf %lf %lf\n",
    //      res[0], res[1], res[2], res[3], res[4]);
    return res;
  }

  if (pts.size() == 0) {
    return vector<double>(5,0);
  }

  // Now let us solve LP on all points in pts
  set<int> pts_set;
  FOREACH(it, pts) pts_set.insert(*it);
  FILE *fin = fopen("solver.in", "w");
  assert(fin != NULL);

  fprintf(fin, "%lf %lf %lf %lf\n", xlo, xhi, ylo, yhi);
  fprintf(fin, "%d\n", (int) pts.size());

  FOR(q, (int) pts.size()) {
    int v = pts[q];
    PII b = points.get_bucket_of_point(v);
    
    vector<pair<double, int> > edge4;
    FOR(__tmp, 4) edge4.pb(mp(1e30, -47));

    FOREACH(it2, points.buckets[b.fi][b.se]) {
      if (pts_set.count(*it2) == 0) {
        edge4.push_back(mp(points.dist(v, *it2), *it2));
      }
    }
    sort(edge4.begin(), edge4.end());
    while (edge4.size() > 4) edge4.pop_back();
    assert(edge4.size() == 4);

    FOR(bbx, points.BUCKETS_X) {
      FOR(bby, points.BUCKETS_Y) {
        if (b.fi == bbx && b.se == bby) {
          continue; // bucket already processed
        }
        if (points.min_distance_point_bucket(v, mp(bbx, bby) ) > edge4[3].fi)  {
          continue; // bucket too far
        }

        FOREACH(it2, points.buckets[bbx][bby]) {
          if (pts_set.count(*it2) == 0) {
            edge4.push_back(mp(points.dist(v, *it2), *it2));
          }
        }
        sort(edge4.begin(), edge4.end());
        while (edge4.size() > 4) edge4.pop_back();
        assert(edge4.size() == 4);
      }
    }
    fprintf(fin, "%lf %lf %lf %lf\n", edge4[0].fi, edge4[1].fi, edge4[2].fi, edge4[3].fi);
  }

  FOR(q, (int) pts.size()) {
    FOR(w, (int) pts.size()) {
      fprintf(fin, "%lf ", points.dist(pts[q], pts[w]));
    }
    fprintf(fin, "\n");
  }
  fclose(fin);
  int ret = system("./solver.py");
  assert(ret == 0);
  FILE* fout = fopen("solver.out", "r");
  double c0, c1, c2, c3, c4;
  fscanf(fout, "%lf %lf %lf %lf %lf", &c0, &c1, &c2, &c3, &c4);
  fclose(fout);
  no_processed += pts.size();
 

  printf(">>>>>>>>>>>>>>> %lf%% done\n", no_processed * 100.0 / N);
  printf("----------------------------------------------------------------------------");
  printf("\n\n\n\n");
  //printf(" %lf %lf %lf %lf %lf\n", c0, c1, c2, c3, c4);
  vector<double> cost;
  cost.pb(c0);
  cost.pb(c1);
  cost.pb(c2);
  cost.pb(c3);
  cost.pb(c4);
  return cost;
}

void solve_lowerbound() {
  no_processed = 0;
  vector<int> pts;
  FOR(q, N) pts.push_back(q);
  vector<double> len = recursive_lowerbound(pts, -0.1, 20000.7, -0.1, 20000.7);
  
  double min4 = 0;
  vector<double> tmp;
  FOR(q, N) tmp.pb(points.adjust[q]);
  sort(tmp.begin(), tmp.end());
  
  FOR(q, 4)
    min4 += tmp[q];
  printf("cycles: %lf \n paths: %lf (raw %lf)\n", len[0] / 4, len[4] / 4 + min4, len[4] / 4);
}

int main() {
  printf("loading points\n");
  points.load();
  printf("points loaded\n");
  load_adjust();
  normalize();
  solve_lowerbound();
}
