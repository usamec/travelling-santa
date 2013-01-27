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
    return sqrt(dx*dx + dy*dy + 1e-20);
    
  }

  double dist_e2(int a, int b) {
    double dx = posx[a] - posx[b];
    double dy = posy[a] - posy[b];
    return sqrt(dx*dx + dy*dy + 1e-20);
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


double gradient[N];
double oldgradient[N];
typedef pair<double, int> PDI;

vector<PDI> incoming[N];
vector<PDI> outgoing[N];

/* Proof:
 - Let 2TSP be any solution for 2 edge-disjoint TSP paths.
   Cost of the solution is the length of the longer path.

 - Similarly, let 2HAM be a solution for 2 edge-disjoint hamiltonian cycles
   with the cost of longer one.
 
 - Let 4EDGE be any solution for the problem
  "Find a set of edges such that each vertex is incident to 4 edges".
   Cost of the solution is the sum of all edge lengths.
 
 - Let 4DIR be any solution for the problem
   "Find a set of 4 directed edges of from each vertex" with the cost being the sum
   of all lengths.

 Lemma 1: OPTCOST(4DIR) / 2 <= OPTCOST(4EDGE).
   Let 4EDGE be any solution of the problem.
   Then each edge in 4EDGE can be decomposed into two directed edges.
   Note that this transformation produces 4DIR solution of twice the cost.

 Note 1: OPT(4DIR) can be trivially calculated (calculate for each vertex separately).

 Lemma 2: OPTCOST(4EDGE) / 2 <= OPTCOST(2HAM).
   Similarly to lemma 1, each solution for 2HAM is a solution for 4EDGE.
   Moreover, the cost of 2HAM (maximum path) is at least the cost of average 
   length of these cycles which is a cost of 4EDGE/2.
   
 Trick 1: Clever transformation.
  If we would like to estimate the lower bound on OPTCOST(2HAM), we can use following trick:

 - Let D_ij = C_ij + P_i + P_j be a new length matrix for edges.
  we may conclude that
  - each solution of 2HAM on C is also a solution for 2HAM on D and vice versa.
  - moreover, assuming that SUM(P_i) = 0, their cost is exactly the same!
    (the difference is 4 * SUM(P_i)).
  Thus OPTCOST(2HAM, C) = OPTCOST(2HAM, D).

 Therefore, we have a following chain of inequalities: 
  for each cost matrix D (vector P): 
    OPTCOST(4DIR, D) / 4  <= OPTCOST(4EDGE, D) /2 <= OPTCOST(2HAM, D)
  and OPTCOST(2HAM, D) = OPTCOST(2HAM, C). 
  Thus, in order to estimate lower bound of OPTCOST(2HAM, C), we may estimate
  MAX(OPTCOST(4DIR, D) for all P). This can be done by a subgradient optimization.

 Now, let us worry about 2TSP:
  The trick which worked nicely for Hamiltonian cycles is a bit harder to use for paths.
  The reason is that the transformation between C and D will now change the cost of the
  solution (in fact, cost(2TSP, D) = cost(2TSP, C) - P_i - P_j - P_k - P_l).
  This also means that OPT(2TSP, D) != OPT(2TSP, C).

  Nevertheless, we can overcome these problems in a following way:
  cost(2TSP, D) + P_i + P_j + P_k + P_l = cost(2TSP, C) for each solution 2TSP (note that the
  solution is same on both sides of the equation!)
  cost(2TSP, D) + sum(min4(P)) <= cost(2TSP, C)
  cost(OPT(2TSP, C), D) + sum(min4(P)) <= cost(OPT(2TSP, C), C) = OPTCOST(2TSP, C)
  OPTCOST(2TSP, D) + sum(min4(P)) <= OPTCOST(2TSP, C).

  Assuming that 4 minimal values in P are of reasonable size (they will be negative), this can be a good estimate.

  Now, we can continue similarly to what we did for hamiltonian cycles.
  Formally, let us define two new problems:

  Let 4EDGE2MISS be a set of 4 edges from each vertext with the exception of two missing edges.
  Similarly, Let 4DIR4MISS be a set of 4 directed edges from each vertex with the exception that there might be
  4 directed edges missing. Then
  OPTCOST(4DIR4MISS, D) / 4 <= OPT(4EDGE2MISS, D) /2 <= OPT(2TSP, D)

  Therefore, in order to obtain a good lower bound on 2TSP, we can maximize the following quantity:
  OPT(4DIR4MISS, D) / 4 + sum(min4(P)) / 4.

  Note that calculating 4DIR4MISS is easy, we just use 4 smallest edges from each vertext and then
  globally erase 4 longest edges.

*/
double ascent(double step) {
  FOR(q, N) oldgradient[q] = gradient[q];
  FOR(q, N) gradient[q] *= .4; // dampen oscillations + boost single direction
  FOR(q, N) incoming[q].clear();
  FOR(q, N) outgoing[q].clear();

  double sum = 0;
  vector<PII> edges;

  points.update_minadjust();

  FOR(bx, points.BUCKETS_X) {
    FOR(by, points.BUCKETS_Y) {
      FOREACH(it, points.buckets[bx][by]) {

        int v = *it;
        vector<pair<double, int> > min4;
        FOREACH(it2, points.buckets[bx][by]) {
          if (*it2 == v) continue;
          min4.push_back(mp(points.dist(v, *it2), *it2));
        }

        sort(min4.begin(), min4.end());
        while (min4.size() > 5) min4.pop_back();
        assert(min4.size() == 5);

        FOR(bbx, points.BUCKETS_X) {
          FOR(bby, points.BUCKETS_Y) {
            if (bx == bbx && by == bby)  continue;
            if (points.min_distance_point_bucket(v, mp(bbx, bby) ) > min4[3].fi)  continue;
            
            FOREACH(it2, points.buckets[bbx][bby]) {
              min4.push_back(mp(points.dist(v, *it2), *it2));
            }
            sort(min4.begin(), min4.end());
            while (min4.size() > 5) min4.pop_back();
            assert(min4.size() == 5);
          }
        }
        assert(min4.size() == 5);
        FOR(q, 4) {
          edges.push_back(mp(v, min4[q].se));
          sum += min4[q].fi;
          incoming[min4[q].se].pb(mp(min4[q].fi, v));
        }
        outgoing[v] = min4;
      }
    }
  }

  map<int, int> distrib;
 
  FOR(q, N) {
    distrib[incoming[q].size()]++;
  }
  printf("Indegree distribution:\n");
  FOREACH(it, distrib) {
    printf("%d: %d\n", it->fi, it->se);
  }

  vector<int> badness(N);

  int total_badness = 0;
  FOR(q, N) {
    int extra_cnt = 0;
    FOR(tmp, (int) incoming[q].size()) {
      int w = incoming[q][tmp].se;
      int in_out = 0;
      FOR(tmp2, 4)
        if (w == outgoing[q][tmp2].se)
          in_out = 1;
      extra_cnt += !in_out;
    }
    badness[q] = extra_cnt;
    total_badness += badness[q];
    // extra_cnt = ((int) incoming[q].size() - 4);
    gradient[q] += step * extra_cnt; //(incoming[q].size() - 4.0);
    FOR(w, 4)
      gradient[outgoing[q][w].se] += step * 0.2 * extra_cnt; //0.1 * (incoming[q].size() - 4.0);
  }

  printf("Total badness: %d\n", total_badness);
  /*
  FILE *edgefile = fopen("edge.dat", "w");
  assert(edgefile != NULL);
  FOR(q, N) {
    fprintf(edgefile, "%d %d %d %d\n",
          outgoing[q][0].se,
          outgoing[q][1].se,
          outgoing[q][2].se,
          outgoing[q][3].se
        );
  }
  fclose(edgefile);
  
  FILE *plotfile = fopen("plot.dat", "w");
  assert(plotfile != NULL);
  FOREACH(it, edges) {
    fprintf(plotfile, "%lf %lf %lf %lf\n", points.posx[it->fi], points.posy[it->fi],
                points.posx[it->se], points.posy[it->se]);
  }
  fclose(plotfile);
  
  FILE *badfile = fopen("bad1.dat", "w");
  assert(badfile != NULL);
  FOR(q, N) {
    fprintf(badfile, "%lf %lf %d\n", points.posx[q], points.posy[q], badness[q]);
  }
  fclose(badfile);

  FILE *adjplotfile = fopen("adjplot.dat", "w");
  assert(adjplotfile != NULL);
  FOR(q, N) {
    fprintf(adjplotfile, "%lf %lf %lf\n", points.posx[q], points.posy[q], points.adjust[q]);
  }
  fclose(adjplotfile);

  FILE *gradientfile = fopen("gradient.dat", "w");
  assert(gradientfile != NULL);
  FOR(q, N) {
    fprintf(adjplotfile, "%lf %lf %lf\n", points.posx[q], points.posy[q], gradient[q]);
  }
  fclose(gradientfile);
  
  FILE *badfile = fopen("bad2.dat", "w");
  assert(badfile != NULL);
  FOR(q, N) {
    if (incoming[q].size() != 4) {
      fprintf(badfile, "%lf %lf %d\n", points.posx[q], points.posy[q], (int)incoming[q].size() - 4);
    }
  }
  fclose(badfile);
 */
  { 
    vector<double> tmp;
    FOR(q, N) tmp.pb(points.adjust[q]);
    sort(tmp.begin(), tmp.end());
    FOR(q, N) if (points.adjust[q] < tmp[10]) gradient[q] += 2 * step;
  }

  double gsum = 0;
  FOR(q, N)
    gsum += gradient[q];
  double gadjust = -gsum / N;
  FOR(q, N)
    points.adjust[q] += gradient[q] + gadjust;

  vector<double> all_outgoing_dist;
  FOR(q, N) {
    FOREACH(it, outgoing[q])
      all_outgoing_dist.push_back(it->fi);
  }
  sort(all_outgoing_dist.begin(), all_outgoing_dist.end());
  reverse(all_outgoing_dist.begin(), all_outgoing_dist.end());
  double longest4 = 0;
  FOR(q, 4)
    longest4 += all_outgoing_dist[q];

  //printf("circle->path: %lf\n", circle_to_path);
  double min4 = 0;
  vector<double> tmp;
  FOR(q, N) tmp.pb(points.adjust[q]);
  sort(tmp.begin(), tmp.end());
  FOR(q, 4)
    min4 += tmp[q];
 
  printf("raw: %lf\n", sum/4);
  printf("circle->path  %lf\n", -longest4/4);
  printf("transfo patch %lf\n", min4/4);
  return (sum - longest4) / 4 + min4/4;
}

void save_adjust() {
  FILE *f = fopen("adjust.dat", "w");
  FOR(q, N) fprintf(f, "%lf\n", points.adjust[q]);
  fclose(f);
}

void load_adjust() {
  FILE *f = fopen("adjust.dat", "r");
  if (f != NULL) {
    FOR(q, N) fscanf(f, "%lf", &points.adjust[q]);
    fclose(f);
  } else {
    // empty adjust
    FOR(q, N) points.adjust[q] = 0;
  }
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

int main() {
  srand(time(NULL));
  printf("loading points\n");
  points.load();
  printf("points loaded\n");
  load_adjust();
  FOR(q, N) gradient[q] = 0;

  double step = 0.5;
  int period = 500;
  while (period) {
    FOR(i, period) {
      normalize();

      double res = ascent(step);
      FILE *flog = fopen("progress.log", "a");
      assert(flog != NULL);
      fprintf(flog, "%lf\n", res - 6500000);
      fclose(flog);

      printf("current lower bound for paths: %lf\n", res);
      save_adjust();
    }
  period = period/2;
  step /= 5;
  }
  
}
