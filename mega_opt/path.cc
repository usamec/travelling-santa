#include "path.h"

#include <tr1/unordered_map>
#include <queue>

const int kMaxKicks = 5;
vector<double> good_non_seq;
int visits;

bool IsPairIn(pair<int, int> p, const vector<pair<int, int> >& v) {
  for (int i = 0; i < v.size(); i++) {
    if (v[i] == p) return true;
    if (v[i] == make_pair(p.second, p.first)) return true;
  }
  return false;
}

static unordered_map<pair<int, int>, double> cache;

inline double GetDistance(int p1, int p2, const vector<pair<double, double> >& points) {
/*  if (cache.count(make_pair(p1, p2))) {
    return cache[make_pair(p1, p2)];
  }*/
  double d = CalcDist(points[p1], points[p2]);
//  cache[make_pair(p1, p2)] = d;
  return d;
}

pair<double, pair<int, int> > Path::GetBest2OptMove(const Path& blocking_path, int p1, int p2,
    const vector<vector<int> >& closest,
    const vector<pair<double, double> >& points) const {
//  if (next_[p1] != p2 && next_[p2] != p1);

  if (next_[p2] == p1) swap(p1, p2);

  // Now we have next_[p1] == p2

  double best_improv = 0;
  int best_p3 = -1;

  for (int i = 0; i < closest[p1].size(); i++) {
    int p3 = closest[p1][i];
    if (p3 == p2) continue;
    int p4 = next_[p3];
    if (p4 == -1) continue;
    if (p4 == p1) continue;
    if (p4 == p2) continue;

    // Now we exchange p1-p2, p3-p4 for p1-p3 p2-p4
    if (blocking_path.HasEdge(p1, p3)) continue;
    if (blocking_path.HasEdge(p2, p4)) continue;
    double cost_prev = GetDistance(p1, p2, points) + GetDistance(p3, p4, points);
    double cost_new = GetDistance(p1, p3, points) + GetDistance(p2, p4, points);
    double improv = cost_prev - cost_new;
    if (best_p3 == -1 || improv > best_improv) {
      best_improv = improv;
      best_p3 = p3;
    }
  }

  return make_pair(best_improv, make_pair(p1, best_p3));
}

void Path::Make2OptMove(pair<int, int> p1p3) {
  int p1 = p1p3.first;
  int p2 = next_[p1];
  int p3 = p1p3.second;
  int p4 = next_[p3];
  if (ip_[p3] > ip_[p1]) {
    Reverse(p2, p3);
  } else {
    Reverse(p4, p1);
  }
}

pair<double, vector<int> > Path::GetBest3OptMove(const Path& blocking_path, int p1, int p2,
    const vector<vector<int> >& closest,
    const vector<pair<double, double> >& points) const {
//  if (next_[p1] != p2 && next_[p2] != p1);

  if (next_[p2] == p1) swap(p1, p2);

  double best_improv = 0;
  vector<int> best_path;

  // swap edges
  // p1->p2, p3->p4, p5->p6:
  // p1->p3, p2->p5, p4->p6
  for (int i = 0; i < closest[p1].size(); i++) {
    int p3 = closest[p1][i];
    if (p3 == p2) continue;
    int p4 = next_[p3];
    if (p4 == -1) continue;
    if (p4 == p1) continue;
    for (int j = 0; j < closest[p2].size(); j++) {
      int p5 = closest[p2][j];
      
      if (p5 == p1 || p5 == p3 || p5 == p4) continue;
      int p6 = next_[p5];
      if (p6 == -1 || p6 == p1 || p6 == p2 || p6 == p3 || p6 == p4) continue;
      int i1 = ip_[p1];
      int i3 = ip_[p3];
      int i5 = ip_[p5];
      if (i1 < i5 && i5 < i3) continue;
      if (i3 < i1 && i1 < i5) continue;
      if (i5 < i3 && i3 < i1) continue;
      if (blocking_path.HasEdge(p1, p3)) continue;
      if (blocking_path.HasEdge(p2, p5)) continue;
      if (blocking_path.HasEdge(p4, p6)) continue;
      double cost_prev = GetDistance(p1, p2, points) + GetDistance(p3, p4, points) + 
                         GetDistance(p5, p6, points);
      double cost_new = GetDistance(p1, p3, points) + GetDistance(p2, p5, points) + 
                        GetDistance(p4, p6, points);

      double improv = cost_prev - cost_new;
      if (best_path.empty() || improv > best_improv) {
        vector<int> pp;
        pp.push_back(p1); pp.push_back(p2); pp.push_back(p3);
        pp.push_back(p4);; pp.push_back(p5); pp.push_back(p6);
        best_improv = improv;
        best_path = pp;
      }
    }
  }
  return make_pair(best_improv, best_path);
}

void Path::Make3OptMove(vector<int> pp) {
  int i1 = ip_[pp[0]];
  int i3 = ip_[pp[2]];
  int i5 = ip_[pp[4]];

  if (i1 < i3 && i3 < i5) {
    Reverse(pp[1], pp[2]);
    Reverse(pp[3], pp[4]);
  } else if (i3 < i5 && i5 < i1) {
    Reverse(pp[3], pp[4]);
    Reverse(pp[4], pp[0]);
  } else {
    Reverse(pp[1], pp[2]);    
    Reverse(pp[5], pp[1]);
  }
}

KOptMove Path::GetKOptMove(const Path& blocking_path, int p1, int p2,
    const vector<vector<int> >& closest,
    const vector<pair<double, double> >& points, int Klim, double allowed_badness) {
  if (next_[p2] == p1) swap(p1, p2);
 
  KOptMove move;
  move.out.push_back(make_pair(p1, p2));
  move.improv = GetDistance(p1, p2, points);
  unordered_map<int, int> dist;
  GetVertexesInDistance(p1, Klim+2, closest, dist);

  GetKOptMoveExpandY(blocking_path, p1, p2, closest, points, Klim-1, move, allowed_badness,
      dist);
  if (!move.valid)
    GetKOptMoveExpandY(blocking_path, p2, p1, closest, points, Klim-1, move, allowed_badness,
        dist);
  return move;
}

void Path::GetKOptMoveExpandY(const Path& blocking_path, int pb, int p,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness, unordered_map<int, int>& dist) {
  // Klim -> edges_left: 1 -> 3; 2 -> 5; ...
  int edges_left = Klim*2 + 1;
  for (int i = 0; i < closest[p].size(); i++) {
    int p2 = closest[p][i];
/*    if (!dist.empty()) {
      if (dist.count(p2) == 0 || dist[p2] > edges_left) continue;
    }*/
    if (IsPairIn(make_pair(p, p2), move.in)) continue;
    if (IsPairIn(make_pair(p, p2), move.out)) continue;
//    pair<int, int> kicked_zal = move.kicked;
    bool blocked = blocking_path.HasEdge(p, p2);
    bool kick_add = false;
    if (blocked) {
      if (move.kicked.size() < kMaxKicks) {
        move.kicked.push_back(make_pair(p, p2));
        kick_add = true;
        blocked = false;
      }
    }
    if (!blocked) {
      double improv_zal = move.improv;
      move.in.push_back(make_pair(p, p2));
      move.improv -= GetDistance(p, p2, points);
      if (move.improv > 0)
        GetKOptMoveExpandX(blocking_path, pb, p2, closest, points,
                           Klim-1, move, allowed_badness,
                           dist);
      if (move.valid) return;
      move.in.pop_back();
      move.improv = improv_zal;
    }
    if (kick_add)
      move.kicked.pop_back();
  }
}

void Path::GetKOptMoveExpandX(const Path& blocking_path, int pb, int p,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness, unordered_map<int, int>& dist) {
  for (int i = 0; i < 2; i++) {
    int p2;
    if (i == 0) p2 = prev_[p];
    else p2 = next_[p];
    if (p2 == -1) continue;
    if (IsPairIn(make_pair(p, p2), move.in)) continue;
    if (IsPairIn(make_pair(p, p2), move.out)) continue;
    move.out.push_back(make_pair(p, p2));
    double improv_zal = move.improv;
    move.improv += GetDistance(p, p2, points);

    // first try to close path
//    bool can_be_closed = false;
    if (p2 != pb && !IsPairIn(make_pair(p2, pb), move.in) && !IsPairIn(make_pair(p2, pb), move.out)) {
      bool blocked = blocking_path.HasEdge(p2, pb);
      bool kick_add = false;
      if (blocked) {
        if (move.kicked.size() < kMaxKicks) {
          move.kicked.push_back(make_pair(p2, pb));
          kick_add = true; blocked = false;
        }
      }
      if (!blocked) {
        move.in.push_back(make_pair(p2, pb));
        double improv_zal2 = move.improv;
        move.improv -= GetDistance(p2, pb, points);
  //      can_be_closed = IsValidMove(move);
        if (move.improv > 0 && IsValidMove(move)) {
          if (move.kicked.empty()) {
            move.valid = true;
            return;
          } else {
            allowed_ = move.out;
            banned_ = move.in;
            vector<pair<int, int> > kicked = move.kicked;
            pair<double, vector<pair<int, int> > > kick = blocking_path.GetBestPatch(
                *this, kicked, closest, points, -move.improv*allowed_badness);
            if (!kick.second.empty() && move.improv * allowed_badness > -kick.first + 0.01) {
              move.valid = true;
              move.improv += kick.first;
              move.kick_in = kick.second;
              move.kick_out = kicked;
              ClearAllowedBanned();
              return;
            }
            ClearAllowedBanned();
          }
        } else if (move.improv > 1 && Klim >= 2) {
          PatchNonSeqMove(blocking_path, closest, points, Klim, move, allowed_badness);
          if (move.valid) return;
        }
        move.in.pop_back();
        move.improv = improv_zal2;
      }
      if (kick_add)
        move.kicked.pop_back();
    }
    if (Klim != 0) {
      GetKOptMoveExpandY(blocking_path, pb, p2, closest, points, Klim, move,
                         allowed_badness, dist);
      if (move.valid) return;
    }
    move.out.pop_back();
    move.improv = improv_zal;
  }
}

void Path::PatchNonSeqMove(const Path& blocking_path, const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness) {
  vector<pair<vector<int>, int> > cycles = GetCycles(move);
  for (int i = 1; i < cycles.size(); i++) {
    if (cycles[0].first.size() != 2 && cycles[i].first.size() == 2) {
      swap(cycles[0], cycles[i]);
    } else if (cycles[0].first.size() == 2 && cycles[i].first.size() == 2 &&
        cycles[0].second == 2 && cycles[i].second != 2) {
      swap(cycles[0], cycles[i]);
    } else if (cycles[0].first.size() == 2 && cycles[i].first.size() == 2 &&
        cycles[0].second > cycles[i].second && cycles[i].second != 2) {
      swap(cycles[0], cycles[i]);
    }
  }
  if (cycles[0].first.size() == 2 && cycles[0].second != 2 && cycles.size() <= Klim) {
    int ip1 = ip_[cycles[0].first[0]];
    int ip2 = ip_[cycles[0].first[1]];
    if (ip1 > ip2) swap(ip1, ip2);

//    if (ip2 - ip1 > 1000) printf("cycle %d\n", ip2 - ip1);
    for (int i = ip1; i < ip2 && i - ip1 < 6; i++) {
      int p1 = p_[i];
      int p2 = p_[i+1];
      move.out.push_back(make_pair(p1, p2));
      double improv_zal = move.improv;
      move.improv += GetDistance(p1, p2, points);
      unordered_map<int, int> dist;
      GetKOptMoveExpandY(blocking_path, p1, p2, closest, points, Klim, move, allowed_badness,
          dist);
      if (!move.valid)
        GetKOptMoveExpandY(blocking_path, p1, p2, closest, points, Klim, move, allowed_badness,
            dist);
      if (move.valid) {
        move.non_seq = true;
        return;
      }
      move.out.pop_back();
      move.improv = improv_zal;
    }
  }
/*  vector<int> poses;
  for (int i = 0; i < cycles[0].first.size(); i++) {
    poses.push_back(ip_[cycles[0].first[i]]);
  }
  sort(poses.begin(), poses.end());
  int ip1 = poses[0];
  int ip2 = poses[1];
  if (ip2 - ip1 == 1) { return; }

//    if (ip2 - ip1 > 1000) printf("cycle %d\n", ip2 - ip1);
  for (int i = ip1; i < ip2 && i - ip1 < 5; i++) {
//    printf("ii %d\n", i);
    int p1 = p_[i];
    int p2 = p_[i+1];
    if (!IsPairIn(make_pair(p1, p2), move.out)) {
      move.out.push_back(make_pair(p1, p2));
      double improv_zal = move.improv;
      move.improv += GetDistance(p1, p2, points);
      GetKOptMoveExpandY(blocking_path, p1, p2, closest, points, Klim, move, allowed_badness);
      if (!move.valid)
        GetKOptMoveExpandY(blocking_path, p1, p2, closest, points, Klim, move, allowed_badness);
      if (move.valid) {
        move.non_seq = true;
        return;
      }
      move.out.pop_back();
      move.improv = improv_zal;
    }
  }*/
//  printf("patch end\n");
}

vector<pair<vector<int>, int> > Path::GetCycles(const KOptMove& move) const {
  vector<int> outposes;
  unordered_map<int, vector<int> > edges;
  unordered_map<int, vector<int> > edge_lens;
  unordered_map<int, bool> visited;
  for (int i = 0; i < move.out.size(); i++) {
    outposes.push_back(ip_[move.out[i].first]);
    outposes.push_back(ip_[move.out[i].second]);
  }
  outposes.push_back(0);
  outposes.push_back(p_.size()-1);
  sort(outposes.begin(), outposes.end());

  for (int i = 0; i < outposes.size(); i++)
    visited[p_[outposes[i]]] = false;

  for (int i = 0; i < outposes.size(); i+=2) {
    if (outposes[i] != outposes[i+1]) {
      edges[p_[outposes[i]]].push_back(p_[outposes[i+1]]);
      edge_lens[p_[outposes[i]]].push_back(outposes[i+1] - outposes[i]);
      edges[p_[outposes[i+1]]].push_back(p_[outposes[i]]);
      edge_lens[p_[outposes[i+1]]].push_back(outposes[i+1] - outposes[i]);
    }
  }
  for (int i = 0; i < move.in.size(); i++) {
    edges[move.in[i].first].push_back(move.in[i].second);
    edge_lens[move.in[i].first].push_back(1);
    edges[move.in[i].second].push_back(move.in[i].first);
    edge_lens[move.in[i].second].push_back(1);
  }
  int cur = p_[0];
  while (cur != p_.back()) {
    visited[cur] = true;
    int nn = -1;
    for (int i = 0; i < edges[cur].size(); i++) {
      int np = edges[cur][i];
      if (visited[np]) continue;
      nn = np;
      break;
    }
    cur = nn;
  }
  visited[p_.back()] = true;
  vector<pair<vector<int>, int> > ret;
  for (auto it = visited.begin(); it != visited.end(); ++it) {
    if (it->second) continue;
    vector<int> cycle;
    int start = it->first;
    if (edges[start][0] == edges[start][1]) {
      cycle.push_back(start);
      cycle.push_back(edges[start][0]);
      int len = edge_lens[start][0] + edge_lens[start][1];
      ret.push_back(make_pair(cycle, len));
      visited[start] = true;
      visited[edges[start][0]] = true;
      continue;
    }
    cycle.push_back(start);
    int cur = edges[start][0];
    int len = edge_lens[start][0];
    int prev = start;
    while (cur != start) {
      visited[cur] = true;
      cycle.push_back(cur);
      for (int i = 0; i < edges[cur].size(); i++) {
        if (visited[edges[cur][i]] == false && edges[cur][i] != prev) {
          len += edge_lens[cur][i];
          prev = cur;
          cur = edges[cur][i];
          break;
        }
      }
    }
    visited[start] = true;
    ret.push_back(make_pair(cycle, len));
  }
  return ret;
}

bool Path::IsValidMove(const KOptMove& move) const {
  vector<int> outposes;
  for (int i = 0; i < move.out.size(); i++) {
    outposes.push_back(ip_[move.out[i].first]);
    outposes.push_back(ip_[move.out[i].second]);
  }
  outposes.push_back(0);
  outposes.push_back(p_.size()-1);
  sort(outposes.begin(), outposes.end());

  unordered_map<int, int> edges_a;
  unordered_map<int, vector<int> > edges_b;
  for (int i = 0; i < outposes.size(); i+=2) {
    if (outposes[i] != outposes[i+1]) {
      edges_a[p_[outposes[i]]] = p_[outposes[i+1]];
      edges_a[p_[outposes[i+1]]] = p_[outposes[i]];
    }
  }
  for (int i = 0; i < move.in.size(); i++) {
    edges_b[move.in[i].first].push_back(move.in[i].second);
    edges_b[move.in[i].second].push_back(move.in[i].first);
  }
  int cur = p_[0];
  int prev = -1;
  int bcount = 0;
  while (cur != p_.back()) {
    int nn = -1;
    if (edges_a.count(cur)) {
      if (prev != edges_a[cur])
        nn = edges_a[cur];
    }
    if (nn == -1) {
      for (int i = 0; i < edges_b[cur].size(); i++) {
        if (edges_b[cur][i] != prev) {
          nn = edges_b[cur][i]; 
          bcount++;
          break;
        }
      }
    }
    if (nn == -1) return false;
    prev = cur;
    cur = nn;
  }
  return bcount == move.in.size();
}

void Path::ExecuteMove(const vector<pair<int, int> >& in, const vector<pair<int, int> >& out) {
/*  printf("move start\n");
  for (int i = 0; i < out.size(); i++)
    printf("out %d %d %d %d\n", out[i].first, out[i].second,
        ip_[out[i].first], ip_[out[i].second]);
  for (int i = 0; i < in.size(); i++)
    printf("in %d %d %d %d\n", in[i].first, in[i].second,
        ip_[in[i].first], ip_[in[i].second]);*/
  vector<int> outposes;
  for (int i = 0; i < out.size(); i++) {
    outposes.push_back(ip_[out[i].first]);
    outposes.push_back(ip_[out[i].second]);
  }
  outposes.push_back(0);
  outposes.push_back(p_.size()-1);
  sort(outposes.begin(), outposes.end());

  unordered_map<int, int> edges_a;
  unordered_map<int, vector<int> > edges_b;
  for (int i = 0; i < outposes.size(); i+=2) {
    if (outposes[i] != outposes[i+1]) {
      edges_a[p_[outposes[i]]] = p_[outposes[i+1]];
      edges_a[p_[outposes[i+1]]] = p_[outposes[i]];
    }
  }
  for (int i = 0; i < in.size(); i++) {
    edges_b[in[i].first].push_back(in[i].second);
    edges_b[in[i].second].push_back(in[i].first);
  }
  int cur = p_[0];
  int prev = -1;
  int end = p_.back();
  while (cur != end) {
    int nn = -1;
    bool b_move = false;
    if (edges_a.count(cur)) {
      if (prev != edges_a[cur])
        nn = edges_a[cur];
    }
    if (nn == -1) {
      for (int i = 0; i < edges_b[cur].size(); i++) {
        if (edges_b[cur][i] != prev) {
          nn = edges_b[cur][i]; 
          b_move = true;
          break;
        }
      }
    }
    int x2 = next_[cur];
    if (b_move) {
      if (edges_a.count(nn) == 0 || ip_[edges_a[nn]] <= ip_[nn]) {
        Reverse(x2, nn);
      } else {
        int nn2 = edges_a[nn];
        Reverse(x2, nn2);
        Reverse(nn2, nn);
      }
    }
    prev = cur;
    cur = nn;
  }
  //printf("move end\n");
}

void Path::ExecuteMove(const KOptMove& move) {
  ExecuteMove(move.in, move.out);
}

pair<double, vector<pair<int, int> > > Path::GetBestPatch(const Path& blocking_path,
    vector<pair<int, int> >& kicked, const vector<vector<int> >& closest,
    const vector<pair<double, double> >&points, double lower_bound) const {
  if (kicked.size() == 1) {
    double best_improv = -1000000000;
    vector<pair<int, int> > best_patch;
    vector<pair<int, int> > best_kicked;
    pair<double, pair<int, int> > kick2 = GetBest2OptMove(
        blocking_path, kicked[0].first, kicked[0].second, closest, points);
    if (kick2.second.second != -1) {
      best_improv = kick2.first;
      best_patch.clear();
      best_patch.push_back(kick2.second);
      best_patch.push_back(make_pair(next_[kick2.second.first], next_[kick2.second.second]));
      best_kicked = kicked;
      best_kicked.push_back(make_pair(kick2.second.second, next_[kick2.second.second]));
    }
    if (best_improv > lower_bound) {
      kicked = best_kicked;
      return make_pair(best_improv, best_patch);
    }

    pair<double, vector<int> > kick3 = GetBest3OptMove(
        blocking_path, kicked[0].first, kicked[0].second, closest, points);
    if (!kick3.second.empty()) {
      if (kick3.first > best_improv) {
        best_improv = kick3.first;
        best_patch.clear();
        best_patch.push_back(make_pair(kick3.second[0], kick3.second[2]));
        best_patch.push_back(make_pair(kick3.second[1], kick3.second[4]));
        best_patch.push_back(make_pair(kick3.second[3], kick3.second[5]));
        best_kicked = kicked;
        best_kicked.push_back(make_pair(kick3.second[2], kick3.second[3]));
        best_kicked.push_back(make_pair(kick3.second[4], kick3.second[5]));
      }
    }
    if (!best_patch.empty()) {
      kicked = best_kicked;
      return make_pair(best_improv, best_patch);
    }
  }
  if (kicked.size() == 2) {
    int p1 = kicked[0].first;
    int p2 = kicked[0].second;
    if (next_[p2] == p1) swap(p1, p2);
    int p3 = kicked[1].first;
    int p4 = kicked[1].second;
    if (next_[p4] == p3) swap(p3, p4);
    if (!blocking_path.HasEdge(p1, p3) && !blocking_path.HasEdge(p2, p4)) {
      double improv = GetDistance(p1, p2, points) + GetDistance(p3, p4, points) -
          GetDistance(p1, p3, points) - GetDistance(p2, p4, points);

      vector<pair<int, int> > patch;
      patch.push_back(make_pair(p1, p3));
      patch.push_back(make_pair(p2, p4));
      return make_pair(improv, patch);
    }
  }
  if (kicked.size() == 3) {
    vector<int> perm;
    for (int i = 0; i < 3; i++) perm.push_back(i);
    
    double best_improv = -1000000000;
    vector<pair<int, int> > best_patch;
    vector<pair<int, int> > patch;

    do {
      int p1 = kicked[perm[0]].first;
      int p2 = kicked[perm[0]].second;
      int p3 = kicked[perm[1]].first;
      int p4 = kicked[perm[1]].second;
      int p5 = kicked[perm[2]].first;
      int p6 = kicked[perm[2]].second;

      if (next_[p2] == p1) swap(p1, p2);
      if (next_[p4] == p3) swap(p3, p4);
      if (next_[p6] == p5) swap(p5, p6);
      int i1 = ip_[p1];
      int i3 = ip_[p3];
      int i5 = ip_[p5];
      if (i1 < i5 && i5 < i3) continue;
      if (i3 < i1 && i1 < i5) continue;
      if (i5 < i3 && i3 < i1) continue;
      if (blocking_path.HasEdge(p1, p3)) continue;
      if (blocking_path.HasEdge(p2, p5)) continue;
      if (blocking_path.HasEdge(p4, p6)) continue;
      double cost_prev = GetDistance(p1, p2, points) + GetDistance(p3, p4, points) + 
                         GetDistance(p5, p6, points);
      double cost_new = GetDistance(p1, p3, points) + GetDistance(p2, p5, points) + 
                        GetDistance(p4, p6, points);
      double improv = cost_prev - cost_new;
      if (best_patch.empty() || improv > best_improv) {
        patch.clear();
        patch.push_back(make_pair(p1, p3));
        patch.push_back(make_pair(p2, p5));
        patch.push_back(make_pair(p4, p6));
        best_improv = improv;
        best_patch = patch;
      }
    } while (next_permutation(perm.begin(), perm.end()));
    if (!best_patch.empty()) {
      return make_pair(best_improv, best_patch);
    }
  } else {
    vector<int> ips;
    double start_improv = 0;
    for (int i = 0; i < kicked.size(); i++) {
      ips.push_back(ip_[kicked[i].first]);
      ips.push_back(ip_[kicked[i].second]);
      start_improv += GetDistance(kicked[i].first, kicked[i].second, points);
    }
    sort(ips.begin(), ips.end());
    int start = ips[0];
    int end = ips.back();
    vector<pair<int, int> > segments;
    for (int i = 1; i + 1 < ips.size(); i+=2) {
      segments.push_back(make_pair(ips[i], ips[i+1]));
    }
    vector<int> perm;
    for (int i = 0; i < segments.size(); i++)
      perm.push_back(i);
    double best_improv = -1000000000;
    vector<pair<int, int> > best_patch;
    vector<pair<int, int> > patch;
    do {
      for (int o = 0; o < (1 << segments.size()); o++) {
        int last = start;
        bool good = true;
        double improv = start_improv;
        patch.clear();
        for (int j = 0; j < segments.size(); j++) {
          int e1, e2;
          if ((o & (1 << j)) == 0) {
            e1 = segments[perm[j]].first;
            e2 = segments[perm[j]].second;
          } else {
            e2 = segments[perm[j]].first;
            e1 = segments[perm[j]].second;
          }
          if (blocking_path.HasEdge(p_[last], p_[e1])) {
            good = false;
            break;
          }
          improv -= GetDistance(p_[last], p_[e1], points);
          if (improv < lower_bound) {
            good = false;
            break;
          }
          patch.push_back(make_pair(p_[last], p_[e1]));
          last = e2;
        }
        if (blocking_path.HasEdge(p_[last], p_[end])) good = false;
        if (good) {
          improv -= GetDistance(p_[last], p_[end], points);
          patch.push_back(make_pair(p_[last], p_[end]));
          if (improv > best_improv) {
            best_improv = improv;
            best_patch = patch;
          }
        }
      }
    } while(next_permutation(perm.begin(), perm.end()));
    if (!best_patch.empty()) {
      return make_pair(best_improv, best_patch);
    }
  }
  return make_pair(-1000000000, vector<pair<int, int> >());
}

void Path::GetVertexesInDistance(int start, int dist_lim, const vector<vector<int> >& closest,
      unordered_map<int, int>& dist) const {
  queue<int> fr;
  fr.push(start);
  dist[start] = 0;
  while (!fr.empty()) {
    int x = fr.front(); fr.pop();
    int dx = dist[x];
    if (dx == dist_lim) continue;
    for (int i = 0; i < closest[x].size(); i++) {
      int n = closest[x][i];
      if (dist.count(n) == 0) {
        dist[n] = dx + 1;
        fr.push(n);
      }
    }
    int n = next_[x];
    if (dist.count(n) == 0) {
      dist[n] = dx + 1;
      fr.push(n);
    }
    n = prev_[x];
    if (dist.count(n) == 0) {
      dist[n] = dx + 1;
      fr.push(n);
    }
  }
}

void PrintPaths(const Path& p1, const Path& p2, string filename) {
  FILE* f = fopen(filename.c_str(), "w");
  for (int i = 0; i < p1.p_.size(); i++)
    fprintf(f, "%d,%d\n", p1.p_[i], p2.p_[i]);
  fclose(f);
}
