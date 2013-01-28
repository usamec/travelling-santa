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
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

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

class LowerBoundTask {
  public:
    LowerBoundTask(const vector<int>& pts) {
      this->pts = pts;
      this->id = __id++;
      this->pid = -1;
      sprintf(this->in_file, "task%d.in", this->id);
      sprintf(this->out_file, "task%d.out", this->id);
    }

    int size() {
      return pts.size();
    }

    pair<vector<double>, vector<PII> > run_sync() {
      generate_in_file();
      char cmdline[1000];
      sprintf(cmdline, "./edges.py %s %s", this->in_file, this->out_file);

      int ret = system(cmdline);
      assert(ret == 0);
      pair<vector<double>, vector<PII> > res = read_result();
      cleanup();
      return res;
    }

    void run_async() {
      pid_t cpid;
      cpid = fork();
      assert(cpid != -1); // failed fork?
      if (cpid == 0) {
        generate_in_file();
        // exec
        char* program = "./edges.py";
        char in_file[1000]; sprintf(in_file, "task%d.in", this->id);
        char out_file[1000]; sprintf(out_file, "task%d.out", this->id);
        char* args[] = { program, in_file, out_file, NULL};
        char* const env[] = {NULL};
        execve(program, args, env);
        assert(false); // should not be here if program launched
      } else {
        // in parent
        this->pid = cpid;
      }
    }

    int is_running() {
      if (this->pid == -2) return false;
      assert(this->pid != -1);
      int status;
      int t = waitpid(pid, &status, WNOHANG);
      if (t == 0) return true; // still running
      assert(t == pid);
      this->pid = -2;

      if (WIFSIGNALED(status)) {
        // if the child was killed, fail miserably
        printf("STATUS: %x\n", status);
        assert(false);
      }
      if (WIFEXITED(status)) {
        if (WEXITSTATUS(status)) {
          assert(false); // child not exitted cleanly
        }
        return false; // dinished
      }
      assert(false); // wtf?
    }

    void cleanup() {
      unlink(in_file);
      unlink(out_file);
    }
   
    pair<vector<double>, vector<PII> > finish_async() {
      assert(!is_running());
      pair<vector<double>, vector<PII> > res = read_result();
      cleanup();
      return res;
    }
  private:
    pair<vector<double>, vector<PII> > read_result() {
      FILE* fout = fopen(this->out_file, "r");
      printf("read %s\n", this->out_file);
      double c0, c1, c2, c3, c4;
      fscanf(fout, "%lf %lf %lf %lf %lf", &c0, &c1, &c2, &c3, &c4);

      vector<double> cost;
      cost.pb(c0);
      cost.pb(c1);
      cost.pb(c2);
      cost.pb(c3);
      cost.pb(c4);

      vector<PII> edges;
      int a,b;
      FOR(q, (int) pts.size() * 4) {
        fscanf(fout, "%d %d", &a, &b);
        edges.pb(mp(a,b));
      }
      fscanf(fout, "%d %d", &a, &b);
      assert(a==-1 && b==-1);

      fclose(fout);
      return mp(cost, edges);
    }
    
    void generate_in_file() {
      FILE *fin = fopen(this->in_file, "w");
      assert(fin != NULL);
      fprintf(fin, "task %d\n", id);
      fprintf(fin, "%d\n", (int) pts.size());

      set<int> pts_set;
      FOREACH(it, pts) pts_set.insert(*it);

      FOR(q, (int) pts.size()) {
        fprintf(fin, "%d\n", pts[q]);
      }

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
        fprintf(fin, "%d %d %d %d\n", edge4[0].se, edge4[1].se, edge4[2].se, edge4[3].se);
      }

      FOR(q, (int) pts.size()) {
        FOR(w, (int) pts.size()) {
          fprintf(fin, "%lf ", points.dist(pts[q], pts[w]));
        }
        fprintf(fin, "\n");
      }
      FOR(q, (int) pts.size()) {
        fprintf(fin, "%lf %lf\n", points.posx[pts[q]], points.posy[pts[q]]);
      }
      fclose(fin);

    }

  private:
    pid_t pid;
    int id;
    char in_file[1000];
    char out_file[1000];
    static int __id;
    vector<int> pts;
};

int LowerBoundTask::__id = 0;

const int MAX_TASK_SIZE = 180;
vector<LowerBoundTask*> tasks;

void recursive_create_tasks(const vector<int>& pts, double xlo, double xhi, double ylo, double yhi) {
  // sanity checks
  FOR(q, (int) pts.size()) {
    int v = pts[q];
    assert(points.posx[v] >= xlo);
    assert(points.posx[v] <= xhi);
    assert(points.posy[v] >= ylo);
    assert(points.posy[v] <= yhi);
  }

  if ((int) pts.size() > MAX_TASK_SIZE) {
    if (xhi - xlo >= yhi - ylo) {
      vector<double> x;
      FOR(q, (int) pts.size()) x.pb(points.posx[pts[q]]);
      sort(x.begin(), x.end());
      int m = x.size() / 3 + rand() % (1 * x.size() / 3);
      double xmid = x[m] + 0.47;

      vector<int> tmp[2];
      FOR(q, (int) pts.size()) {
        int v = pts[q];
        tmp[points.posx[v] > xmid].push_back(v);
      }
      recursive_create_tasks(tmp[0], xlo, xmid, ylo, yhi);
      recursive_create_tasks(tmp[1], xmid, xhi, ylo, yhi);

    } else {
      vector<double> y;
      FOR(q, (int) pts.size()) y.pb(points.posy[pts[q]]);
      sort(y.begin(), y.end());
      int m = y.size() / 3 + rand() % (1 * y.size() / 3);
      double ymid = y[m] + 0.47;

      vector<int> tmp[2];
      FOR(q, (int) pts.size()) {
        int v = pts[q];
        tmp[points.posy[v] > ymid].push_back(v);
      }
      recursive_create_tasks(tmp[0], xlo, xhi, ylo, ymid);
      recursive_create_tasks(tmp[1], xlo, xhi, ymid, yhi);
    }
  } else if (pts.size() != 0) {
    LowerBoundTask* task = new LowerBoundTask(pts);
    tasks.push_back(task);
  }
}


const int MAX_RUNNING_TASKS = 4; // tasks will be big

void solve_lowerbound() {
  vector<int> pts;
  vector<PII> edges;
  FOR(q, N) pts.push_back(q);
  printf("preparing tasks\n");
  recursive_create_tasks(pts, -0.1, 20000.7, -0.1, 20000.7);
  random_shuffle(tasks.begin(), tasks.end());
  printf("preparing done, running tasks\n");
  vector<double> res = vector<double>(5, 0); // initial 
  int no_processed = 0;
  
  // synchronous version
  /*
  while (!tasks.empty()) {
    LowerBoundTask task = tasks.back(); tasks.pop_back();
    vector<double> tmp = task.run_sync();
    res = merge_results(res, tmp);
    no_processed += task.size();
    double _done = no_processed * 1.0 / N;
    printf(">>>>> %lf%% done, current size %lf, expected %lf\n",
        _done * 100, res[0] / 4, res[0] / 4 / (_done + 1e-6));
  }
  */

  FILE* fedge = fopen("edge.dat", "w");

  vector<LowerBoundTask*> running_tasks;

  while (!tasks.empty() || !running_tasks.empty()) {
    while ((int)running_tasks.size() < MAX_RUNNING_TASKS && !tasks.empty()) {
      LowerBoundTask* task = tasks.back(); tasks.pop_back();
      task->run_async();
      running_tasks.push_back(task);
    }

    sleep(1);
    FOR(q, (int) running_tasks.size()) {
      LowerBoundTask* task = running_tasks[q];
      if (task->is_running()) {
        continue;
      }
    
      no_processed += task->size();
      pair<vector<double>, vector<PII> > tmp = task->finish_async();
      res = merge_results(res, tmp.fi);
      FOREACH(it, tmp.se) {
        fprintf(fedge, "%d %d\n", it->fi, it->se);
        edges.pb(*it);
      }

      double _done = no_processed * 1.0 / N;
      printf(">>>>> %lf%% done, current size %lf, expected %lf\n",
          _done * 100, res[0] / 4, res[0] / 4 / (_done + 1e-6));

      running_tasks.erase(running_tasks.begin() + q);
      delete task;
      q--; // go back one step as there will be immediate q++
    }
  }

  fclose(fedge);
  double min4 = 0;
  vector<double> tmp;
  FOR(q, N) tmp.pb(points.adjust[q]);
  sort(tmp.begin(), tmp.end());
  
  FOR(q, 4)
    min4 += tmp[q];
  printf("RESULT: cycles: %lf paths: %lf (raw paths %lf)\n",
                    res[0] / 4, res[4] / 4 + min4, res[4] / 4);
}

int main() {
  printf("loading points\n");
  points.load();
  printf("points loaded\n");
  load_adjust();
  normalize();
  solve_lowerbound();
}
