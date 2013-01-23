#include "utility.h"

#include <cassert>
#include <queue>
#include <tr1/unordered_map>


vector<vector<int> > EdgeToNext(set<pair<int, int> >& e, int N) {
  vector<vector<int> > next(N);
  for (auto it = e.begin(); it != e.end(); ++it) {
    next[it->first].push_back(it->second);
    next[it->second].push_back(it->first);
  }
  return next;
}

vector<vector<int> > EdgeToNext(unordered_set<pair<int, int> >& e, int N) {
  vector<vector<int> > next(N);
  for (int i = 0; i < N; i++)
    next[i].reserve(2);
  for (auto it = e.begin(); it != e.end(); ++it) {
    next[it->first].push_back(it->second);
    next[it->second].push_back(it->first);
  }
  return next;
}

vector<vector<int> > GetCycles(set<pair<int, int> >&e, int N) {
  vector<bool> used(N, false);
  vector<vector<int> > next = EdgeToNext(e, N);
  vector<vector<int> > ret;
  for (int i = 0; i < N; i++) {
    if (used[i]) continue;
    queue<int> fr;
    fr.push(i);
    used[i] = true;
    vector<int> cycle;
    while (!fr.empty()) {
      int x = fr.front(); fr.pop();
      cycle.push_back(x);
      for (int i = 0; i < next[x].size(); i++) {
        int nn = next[x][i];
        if (used[nn]) continue;
        used[nn] = true;
        fr.push(nn);
        break;
      }
    }
    ret.push_back(cycle);
  }
  return ret;
}

int nextc[100000][3];
int nextcc[100000];

vector<vector<int> > GetCycles(unordered_set<pair<int, int> >&e, int N, int lim = 100000) {
  vector<bool> used(N, false);
//  vector<vector<int> > next = EdgeToNext(e, N);
  for (int i = 0; i < N; i++)
    nextcc[i] = 0;
  for (auto it = e.begin(); it != e.end(); ++it) {
    nextc[it->first][nextcc[it->first]] = it->second;
    nextcc[it->first]++;
    nextc[it->second][nextcc[it->second]] = it->first;
    nextcc[it->second]++;
  }
  vector<vector<int> > ret;
  for (int i = 0; i < N; i++) {
    if (used[i]) continue;
    queue<int> fr;
    fr.push(i);
    used[i] = true;
    vector<int> cycle;
    while (!fr.empty()) {
      int x = fr.front(); fr.pop();
      cycle.push_back(x);
      for (int i = 0; i < nextcc[x]; i++) {
        int nn = nextc[x][i];
        if (used[nn]) continue;
        used[nn] = true;
        fr.push(nn);
        break;
      }
    }
    ret.push_back(cycle);
    if (ret.size() > lim) return ret;
  }
  return ret;
}
unordered_map<int, int> GetCycleunordered_map(vector<vector<int> >& cycles) {
  unordered_map<int, int> ret;
  for (int i = 0; i < cycles.size(); i++) {
    for (int j = 0; j < cycles[i].size(); j++) {
      ret[cycles[i][j]] = i;
    }
  }
  return ret;
}

bool Try(unordered_set<pair<int, int> > paths[2], vector<vector<int> > cycles[2],
    vector<vector<int> > next[2], int start, int cur, unordered_set<pair<int, int> > &erased,
    unordered_set<pair<int, int> > &added, int left, int N,
    unordered_set<pair<int, int> > fix[2],
    unordered_map<int, int> cycle_unordered_map[2], int round,
    bool rozdrb = false, int ms = -1000) {
  if (cur == start && !erased.empty() && erased.size() == added.size()) {
    unordered_set<pair<int, int> > p2[2];
    CopyUnordered(paths[0], &p2[0]);
    CopyUnordered(paths[1], &p2[1]);

    for (auto it = erased.begin(); it != erased.end(); ++it) {
      p2[0].erase(*it);
      p2[1].insert(*it);
    }
    for (auto it = added.begin(); it != added.end(); ++it) {
      p2[1].erase(*it);
      p2[0].insert(*it);
    }

    vector<vector<int> > c2[2];
    c2[0] = GetCycles(p2[0], N, cycles[0].size());
    c2[1] = GetCycles(p2[1], N, cycles[1].size());

    if (c2[0].size() <= cycles[0].size() && c2[1].size() <= cycles[1].size() &&
        (c2[0].size() + c2[1].size() < cycles[0].size() + cycles[1].size() || rozdrb)) {
//      printf("improv %d %d\n", c2[0].size(), c2[1].size());
      CopyUnordered(p2[0], &paths[0]);
      CopyUnordered(p2[1], &paths[1]);
      cycles[0] = c2[0]; cycles[1] = c2[1];
/*      for (auto it = erased.begin(); it != erased.end(); ++it) {
        printf("erase %d %d\n", it->first, it->second);
      }
      for (auto it = added.begin(); it != added.end(); ++it) {
        printf("insert %d %d\n", it->first, it->second);
      }*/
      return true;
    }
    return false;
  }

  if (left == 0) return false;

  int pp = left%2;
  for (int i = 0; i < next[pp][cur].size(); i++) {
    int nn = next[pp][cur][i];
    if (fix[pp].count(make_pair(min(cur,nn), max(cur,nn)))) continue;
    if (round == 1) {
      if (cycle_unordered_map[1-pp][cur] == cycle_unordered_map[1-pp][nn]) {
        continue;
      }
      ms = cur + nn; 
    }
    if (round > 1 && round % 2 == 1) {
      if (cycle_unordered_map[1-pp][cur] != cycle_unordered_map[1-pp][nn]) {
        if (cur + nn < ms) continue;
      }
    }

    if (pp == 0) {
      if (erased.count(make_pair(min(cur,nn), max(cur,nn)))) continue;
      erased.insert(make_pair(min(cur,nn), max(cur,nn)));
      if (Try(paths, cycles, next, start, nn, erased, added, left-1, N, fix,
            cycle_unordered_map,
            round+1, rozdrb, ms))
        return true;
      erased.erase(make_pair(min(cur,nn), max(cur,nn)));
    }
    if (pp == 1) {
      if (added.count(make_pair(min(cur,nn), max(cur,nn)))) continue;
      added.insert(make_pair(min(cur,nn), max(cur,nn)));
      if (Try(paths, cycles, next, start, nn, erased, added, left-1, N, fix,
            cycle_unordered_map, round+1, rozdrb, ms))
        return true;
      added.erase(make_pair(min(cur,nn), max(cur,nn)));
    }
  }
  return false;
}

void Die(char* filename) {
  FILE *f = fopen(filename, "w");
  fprintf(f, "was_optimal\n");
  fclose(f);
}

int main(int argc, char** argv) {
  for (int tr = 0; tr < 10; tr++) {
    FILE *fproblem = fopen(argv[1], "r");
    FILE *ffactor4 = fopen(argv[2], "r");
    double cost; fscanf(fproblem, "%lf", &cost);
    int N; fscanf(fproblem, "%d", &N);
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        double x; fscanf(fproblem, "%lf", &x);
      }
    }
    set<pair<int, int> > fixed1, fixed2;
    int nfixed1; fscanf(fproblem, "%d", &nfixed1);
    for (int i = 0; i < nfixed1; i++) {
      int a, b; double c;
      fscanf(fproblem, "%d %d %lf", &a, &b, &c);
      fixed1.insert(make_pair(min(a, b), max(a, b)));
    }
    int nfixed2; fscanf(fproblem, "%d", &nfixed2);
    for (int i = 0; i < nfixed2; i++) {
      int a, b; double c;
      fscanf(fproblem, "%d %d %lf", &a, &b, &c);
      fixed2.insert(make_pair(min(a, b), max(a, b)));
    }

    fclose(fproblem);

    unordered_set<pair<int, int> > fix[2];
    unordered_set<pair<int, int> > paths[2];
    set<pair<int, int> > free_edges;
    vector<int> deg[2];
    deg[0] = vector<int>(N);
    deg[1] = vector<int>(N);
    int tries = 0;
    for (int i = 0; i < N; i++) {
      int a, b; fscanf(ffactor4, "%d %d", &a, &b);
      deg[0][a]++;
      deg[0][b]++;
      if (fixed1.count(make_pair(b, a))) {
        fixed1.erase(make_pair(b, a));
        fix[0].insert(make_pair(b, a));
        paths[0].insert(make_pair(b, a));
        continue;
      }
      paths[0].insert(make_pair(b, a));
      free_edges.insert(make_pair(b, a));
    }
    for (int i = 0; i < N; i++) {
      int a, b; fscanf(ffactor4, "%d %d", &a, &b);
      deg[1][a]++;
      deg[1][b]++;
      if (fixed2.count(make_pair(b, a))) {
        fixed2.erase(make_pair(b, a));
        fix[1].insert(make_pair(b, a));
        paths[1].insert(make_pair(b, a));
        continue;
      }
      paths[1].insert(make_pair(b, a));
      free_edges.insert(make_pair(b, a));
    }
    assert(fixed1.empty()); 
    assert(fixed2.empty());

    printf("loaded\n");

    vector<vector<int> > cycles[2];
    unordered_map<int, int> cycle_unordered_map[2];
    int limit = 10;
    bool bad = false;
    while (true) {
      for (int i = 0; i < 2; i++) {
        cycles[i] = GetCycles(paths[i], N);
        cycle_unordered_map[i] = GetCycleunordered_map(cycles[i]);
      }
      printf("cycles %d %d\n", cycles[0].size(), cycles[1].size());
      if (cycles[0].size() == 1 && cycles[1].size() == 1)
        break;

      for (int c = 0; c < cycles[0].size(); c++) {
        for (int i = 0; i < cycles[0][c].size(); i++) {
          int ii = (i+1) % cycles[0][c].size();
          int a = cycles[0][c][i], b = cycles[0][c][ii];
          if (fix[0].count(make_pair(min(a,b), max(a,b)))) continue;
          if (free_edges.count(make_pair(min(a,b), max(a,b)))) continue;
          printf("bad edge x0 %d %d %d %d %d %d\n", a, b, i, ii, c, cycles[0][c].size());
          Die(argv[3]);
          return 1;
        }
      }
      for (int c = 0; c < cycles[1].size(); c++) {
        for (int i = 0; i < cycles[1][c].size(); i++) {
          int ii = (i+1) % cycles[1][c].size();
          int a = cycles[1][c][i], b = cycles[1][c][ii];
          if (fix[1].count(make_pair(min(a,b), max(a,b)))) continue;
          if (free_edges.count(make_pair(min(a,b), max(a,b)))) continue;
          printf("bad edge x2 %d %d %d %d %d %d\n", a, b, i, ii, c, cycles[1][c].size());
          Die(argv[3]);
          return 1;
        }
      }

      vector<vector<int> > next[2];
      for (int i = 0; i < 2; i++)
        next[i] = EdgeToNext(paths[i], N);

      bool improv = false;
      vector<int> order;
      for (int i = 0; i < N; i++)
        order.push_back(i);
      random_shuffle(order.begin(), order.end());
      for (int i = 0; i < N; i++) {
        int ii = order[i];
        unordered_set<pair<int, int> > erased, added;
        if (Try(paths, cycles, next, ii, ii, erased, added, limit, N, fix, cycle_unordered_map, 0)) {
          improv = true;
          break;
        }
        if (Try(paths, cycles, next, ii, ii, erased, added, limit, N, fix, cycle_unordered_map, 1)) {
          improv = true;
          break;
        }
      }
      if (!improv) {
        printf("no improv\n");
        if (limit < 24) {
          limit += 2;
          printf("limit increase %d\n", limit);
        } else if (tries < 0) {
          printf("rozdrb\n");
          for (int i = 0; i < 10; i++) {
            vector<vector<int> > next[2];
            for (int i = 0; i < 2; i++)
              next[i] = EdgeToNext(paths[i], N);
            unordered_set<pair<int, int> > erased, added;
            int ss = rand() % N;
            Try(paths, cycles, next, ss, ss, erased, added, 20, N, fix, cycle_unordered_map, 0, true);
          }
          limit = 4;
          tries++;
        } else {
          bad = true;
          break;
        }
      }
    }
    if (bad) continue;
    // Test paths
    assert(cycles[0][0].size() == N);
    assert(cycles[1][0].size() == N);
    int fixed_c = 0;
    for (int c = 0; c < cycles[0].size(); c++) {
      for (int i = 0; i < cycles[0][c].size(); i++) {
        int ii = (i+1) % cycles[0][c].size();
        int a = cycles[0][c][i], b = cycles[0][c][ii];
        if (fix[0].count(make_pair(min(a,b), max(a,b)))) {
          fixed_c++;
          continue;
        }
        if (free_edges.count(make_pair(min(a,b), max(a,b)))) continue;
        printf("bad edge %d %d %d %d %d %d\n", a, b, i, ii, c, cycles[0][c].size());
        Die(argv[3]);
        return 1;
      }
    }
    if (fixed_c != fix[0].size()) {
      Die(argv[3]);
      printf("bad fixed\n"); return 1;
    }
    fixed_c = 0;
    for (int c = 0; c < cycles[1].size(); c++) {
      for (int i = 0; i < cycles[1][c].size(); i++) {
        int ii = (i+1) % cycles[1][c].size();
        int a = cycles[1][c][i], b = cycles[1][c][ii];
        if (fix[1].count(make_pair(min(a,b), max(a,b)))) {
          fixed_c++;
          continue;
        }
        if (free_edges.count(make_pair(min(a,b), max(a,b)))) continue;
        printf("bad edge %d %d %d %d %d %d\n", a, b, i, ii, c, cycles[0][c].size());
        Die(argv[3]);
        return 1;
      }
    }
    if (fixed_c != fix[1].size()) {
      Die(argv[3]);
      printf("bad fixed\n"); return 1;
    }
    printf("output\n");

    FILE *fout = fopen(argv[3], "w");
    fprintf(fout, "new_solution\n");
    for (int i = 0; i < cycles[0][0].size(); i++)
      fprintf(fout, "%d\n", cycles[0][0][i]);
    for (int i = 0; i < cycles[1][0].size(); i++)
      fprintf(fout, "%d\n", cycles[1][0][i]);
    fclose(fout);
    return 0;
  }
  Die(argv[3]);
  return 1;
}
