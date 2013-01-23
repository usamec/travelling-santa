#ifndef ALT_CYCLE_H__
#define ALT_CYCLE_H__

#include "utility.h"

class Graph4 {
 public:
  static vector<int> NormalizePath(vector<int>& p);

  vector<unordered_set<int> > next_;
  vector<unordered_set<int> > good_next_;
  vector<vector<int> >& closest_;
  vector<pair<double, double> >& points_;

  Graph4(vector<int>& path1, vector<int>& path2, vector<vector<int> >& edge4, 
      vector<vector<int> >& closest, vector<pair<double, double> >& points) : 
      next_(path1.size()), good_next_(path1.size()),
      closest_(closest), points_(points) {
    for (int i = 1; i < path1.size(); i++) {
      next_[path1[i]].insert(path1[i-1]);
      next_[path1[i-1]].insert(path1[i]);
    }
    for (int i = 1; i < path2.size(); i++) {
      next_[path2[i]].insert(path2[i-1]);
      next_[path2[i-1]].insert(path2[i]);
    }
    for (int i = 0; i < edge4.size(); i++) {
      for (int j = 0; j < edge4[i].size(); j++) {
        int b = edge4[i][j];
        good_next_[i].insert(b);
        good_next_[b].insert(i);
      }
    }
  }

  vector<int> FindAlternatingCycleDFS(int start, int dfs_limit,
      unordered_set<vector<int> >& banned_cycles) const;

  vector<int> FindAlternatingCycleBellman(int start, int limit,
      unordered_set<vector<int> >& banned_cycles) const;

  bool FindAlternatingCycleDFSExpand(
      vector<int>& cycle, unordered_set<pair<int, int> >& cycle_edge_set, 
      double cur_cost, int dfs_limit,
      unordered_set<vector<int> >& banned_cycles) const;

  vector<int> GetNearPoints(vector<int>& start, int dist_lim) const;
  void Alternate(vector<int>& cycle);
};

#endif
