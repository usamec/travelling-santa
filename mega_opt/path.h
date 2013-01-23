#ifndef PATH_H__
#define PATH_H__

#include "utility.h"
#include "persistent_value.h"
#include <vector>
#include <tr1/unordered_map>

using namespace std;

struct Move2OptWith3Opt {
  double improv;
  pair<int, int> opt2;
  vector<int> opt3;

  Move2OptWith3Opt(double i, pair<int, int> o2, vector<int> o3) :
    improv(i), opt2(o2), opt3(o3) {}
};

struct Move3OptWith3Opt {
  double improv;
  vector<int> opt3;
  vector<int> kick3;

  Move3OptWith3Opt(double i, vector<int> o3, vector<int> k3) :
    improv(i), opt3(o3), kick3(k3) {}
};

struct KOptMove {
  double improv;
  bool valid;
  vector<pair<int, int> > in;
  vector<pair<int, int> > out;
  vector<pair<int, int> > kicked;
  vector<pair<int, int> > kick_in, kick_out;
  bool non_seq;

  KOptMove() : improv(-1), valid(false), non_seq(false) {}
};

class Path {
 public:
  vector<int> next_, prev_;
  vector<int> p_; // Path
  vector<int> ip_; // Inverse mapping (id -> path index)
  vector<pair<int, int> > allowed_, banned_;


  Path(vector<int> p) : next_(p.size(), -1), prev_(p.size(), -1), p_(p), ip_(p.size()),
      allowed_(3, make_pair(-1,-1)), banned_(3, make_pair(-1,-1)) {
    for (int i = 1; i < p.size(); i++) {
      prev_[p[i]] = p[i-1];
      next_[p[i-1]] = p[i];
      ip_[p[i]] = i;
    }
    ip_[p[0]] = 0;
  }

  void ClearAllowedBanned() {
    allowed_.resize(3, make_pair(-1, -1));
    banned_.resize(3, make_pair(-1, -1));
  }

  vector<int> OutputPathToVector() const {
    return p_;
  }

  bool HasEdge(int a, int b) const {
    for (int i = 0; i < allowed_.size(); i++) {
      if (allowed_[i] == make_pair(a, b) || allowed_[i] == make_pair(b, a)) return false;
    }
    for (int i = 0; i < banned_.size(); i++) {
      if (banned_[i] == make_pair(a, b) || banned_[i] == make_pair(b, a)) return true;
    }
    if (next_[a] == b) return true;
    if (next_[b] == a) return true;
    return false;
  }

  void Reverse(int p1, int p2) {
    int ip1 = ip_[p1];
    int ip2 = ip_[p2]+1;
    reverse(p_.begin() + ip1, p_.begin() + ip2);
    for (int i = ip1; i < ip2; i++) {
      ip_[p_[i]] = i;
      if (i > 0) {
        prev_[p_[i]] = p_[i-1];
        next_[p_[i-1]] = p_[i];
      }
      if (i + 1 < p_.size()) {
        prev_[p_[i+1]] = p_[i];
        next_[p_[i]] = p_[i+1];
      }
    }
    next_[p_.back()] = -1;
  }

  double GetLength(const vector<pair<double, double> >& points) const {
    double ret = 0;
    for (int i = 1; i < p_.size(); i++)
      ret += CalcDist(points[p_[i]], points[p_[i-1]]);
    return ret;
  }



  // Return best 2-opt move starting with kicking edge p1-p2
  // return pair<improvemt, <p1, p3> >
  pair<double, pair<int, int> > GetBest2OptMove(const Path& blocking_path, int p1, int p2,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points) const;

  void Make2OptMove(pair<int, int> p1p3);

  pair<double, vector<int> > GetBest3OptMove(const Path& blocking_path, int p1, int p2,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >&points) const;

  void Make3OptMove(vector<int> pp);

  KOptMove GetKOptMove(const Path& blocking_path, int p1, int p2,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, double allowed_badness);

  void GetKOptMoveExpandY(const Path& blocking_path, int pb, int p,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness, unordered_map<int, int>& dist);

  void GetKOptMoveExpandX(const Path& blocking_path, int pb, int p,
      const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness, unordered_map<int, int>& dist);

  pair<double, vector<pair<int, int> > > GetBestPatch(const Path& blocking_path,
      vector<pair<int, int> >& kicked, const vector<vector<int> >& closest,
      const vector<pair<double, double> >&points, double lower_bound) const;

  void PatchNonSeqMove(const Path& blocking_path, const vector<vector<int> >& closest,
      const vector<pair<double, double> >& points, int Klim, KOptMove& move,
      double allowed_badness);

  void ExecuteMove(const KOptMove& move);
  void ExecuteMove(const vector<pair<int, int> >& in, const vector<pair<int, int> >& out);

  bool IsValidMove(const KOptMove& move) const;
  vector<pair<vector<int>, int> > GetCycles(const KOptMove& move) const;

  void GetVertexesInDistance(int start, int dist_lim, const vector<vector<int> >& closest,
      unordered_map<int, int>& dist) const;
};

void PrintPaths(const Path& p1, const Path& p2, string filename);

#endif
