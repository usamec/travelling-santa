int swap_length = 100;
double kill_coef = 2;
double kill_limit = 2000;

bool KillEdge(
    vector<int>&path, unordered_set<pair<int, int> >& used1,
    vector<int>&vpos, double& d, unordered_set<pair<int, int> >& used2,
    int p, int p2,
    unordered_map<int, pair<double, double> >& points,
    double allowed_loss) {

  if (rand() % 5 == 0) return false;
  int pos = min(vpos[p], vpos[p2]);
  if (pos < 3 || pos > path.size() - 5) return false;
  int best_e = -1;
  double best_loss = 0;
  for (int e = pos+2; e < path.size() - 3 && e < pos + 20; e++) {
    // swap pos-pos+1, e-e+1
    // to pos-e, pos+1-e+1
    if (used2.count(make_pair(path[pos], path[e]))) continue;
    if (used2.count(make_pair(path[pos+1], path[e+1]))) continue;
    double dp = CalcDist(points[path[pos]], points[path[pos+1]]) + 
                CalcDist(points[path[e]], points[path[e+1]]);
    double dn = CalcDist(points[path[pos]], points[path[e]]) + 
                CalcDist(points[path[pos+1]], points[path[e+1]]);
    if (dn - dp > allowed_loss) continue;
    if (best_e == -1 || dn - dp < best_loss) {
      best_e = e;
      best_loss = dn - dp;
    }
  }
  if (best_e != -1) {
    int e = best_e;
    used1.erase(make_pair(path[pos], path[pos+1]));
    used1.erase(make_pair(path[pos+1], path[pos]));
    used1.erase(make_pair(path[e], path[e+1]));
    used1.erase(make_pair(path[e+1], path[e]));
    used1.insert(make_pair(path[e], path[pos]));
    used1.insert(make_pair(path[pos], path[e]));
    used1.insert(make_pair(path[e+1], path[pos+1]));
    used1.insert(make_pair(path[pos+1], path[e+1]));
//    d -= dp; d += dn;
    d += best_loss;
    reverse(path.begin() + pos + 1, path.begin() + e + 1);
    for (int j = pos; j < e+2; j++)
      vpos[path[j]] = j;
    return true;
  }
  return false;
}

double UpdateRange(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points, int b, int e) {
  // edges path[b] -> path[b+1], path[e] -> path[e+1]
  // change to path[b] -> path[e], path[b+1] -> path[e+1]
  // sequence path[b+1] .. path[e] gets reversed

  if (used_edges.count(make_pair(path[b+1], path[e+1])))
    return d;
  if (used_edges.count(make_pair(path[b], path[e])))
    return d;

  double dp = CalcDist(points[path[b+1]], points[path[b]]) +
    CalcDist(points[path[e]], points[path[e+1]]);
  double dn = CalcDist(points[path[b]], points[path[e]]) +
    CalcDist(points[path[b+1]], points[path[e+1]]);
  if (dn < dp) {
    d -= dp;
    d += dn;
    used_edges.erase(make_pair(path[b+1], path[b]));
    used_edges.erase(make_pair(path[b], path[b+1]));
    used_edges.erase(make_pair(path[e+1], path[e]));
    used_edges.erase(make_pair(path[e], path[e+1]));
    used_edges.insert(make_pair(path[b], path[e]));
    used_edges.insert(make_pair(path[e], path[b]));
    used_edges.insert(make_pair(path[e+1], path[b+1]));
    used_edges.insert(make_pair(path[b+1], path[e+1]));
    reverse(path.begin() + b + 1, path.begin() + e + 1);
  }
  return d;
}

double OptPath(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  vector<double> ddist;
  for (int b = 1; b < path.size(); b++) {
    ddist.push_back(CalcDist(points[path[b-1]], points[path[b]]));
  }
  sort(ddist.rbegin(), ddist.rend());
  double threshold = ddist[50];
  printf("threshold %lf\n", threshold);
  PrintSeriesStats(ddist, true);
  for (int b = 0; b < path.size() - 1; b++) {
    if (b % 20000 == 0) {
      printf("b %d %lf\n", b, d);
    }
    int max_try = b + swap_length;
    int start = b+2;
    if (CalcDist(points[path[b]], points[path[b+1]]) > threshold) {
      printf("long edge %d %lf\n", b, CalcDist(points[path[b+1]], points[path[b]]));
      start = 0;
      max_try = path.size() - 1;
    }
    for (int e = start; e < max_try && e < path.size() - 1; e++) {
      if (e == b) continue;
      d = UpdateRange(path, d, used_edges, points, min(b, e), max(b, e));
    }
  }
  return d;
}

double sector_div = 557;

int GetPointSector(pair<double, double>& point) {
  int ms = 39;
  int xs = point.first / sector_div;
  if (xs < 0) xs = 0;
  if (xs > ms) xs = ms;
  int ys = point.second / sector_div;
  if (ys < 0) ys = 0;
  if (ys >= ms) ys = ms;
  return xs*(ms+1) + ys;
}

double LocalOptPath(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < path.size() - 1; i++) {
    next[path[i]] = path[i+1];
    prev[path[i+1]] = path[i];
  }
  int first = path[0];

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  for (int i = 1; i < path.size() - 2; i++) {
    if (i % 20000 == 0) {
      printf("%d: %lf\n", i, d);
    }
    int pr = prev[path[i]];
    int nn = next[path[i]];
    int p = path[i];

    if (used_edges.count(make_pair(pr, nn))) continue;

    int best_v = -1;
    double best_improv = 0;
    int s = GetPointSector(points[path[i]]);
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      if (p2 == p) continue;
      if (next[p2] == p) continue;
    //  if (next[p2] == -1) continue;
      // deleted edges:
      // prev[p] -> p, p -> next[p], p2 -> next[p2]
      // new edges:
      // prev[p] -> next[p], p2 -> p, p -> next[p2]
      if (used_edges.count(make_pair(p2, p))) continue;
      if (next[p2] != -1) 
        if (used_edges.count(make_pair(p, next[p2]))) continue;
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn]);
      if (next[p2] != -1)
        dp += CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p]);
      if (next[p2] != -1)
        dn += CalcDist(points[p], points[next[p2]]);
      if (dn < dp && (best_v == -1 || best_improv < dp - dn)) {
        best_v = p2;
        best_improv = dp - dn;
      }
      /*else if (dn < dp*1.1 && rand() % 500 == 0) {
        best_v = p2;
        best_improv = dp - dn;
      }
      else if (dn < dp*1.3 && rand() % 50000 == 0) {
        best_v = p2;
        best_improv = dp - dn;
      }*/
    }
    if (best_v != -1) {
      int p2 = best_v;
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn]);
      if (next[p2] != -1)
        dp += CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p]);
      if (next[p2] != -1) 
        dn += CalcDist(points[p], points[next[p2]]);
      d -= dp; d += dn;
      used_edges.erase(make_pair(prev[p], p));
      used_edges.erase(make_pair(p, prev[p]));
      used_edges.erase(make_pair(next[p], p));
      used_edges.erase(make_pair(p, next[p]));
      used_edges.insert(make_pair(prev[p], next[p]));
      used_edges.insert(make_pair(next[p], prev[p]));
      next[prev[p]] = next[p];
      prev[next[p]] = prev[p];
      used_edges.erase(make_pair(p2, next[p2]));
      used_edges.erase(make_pair(next[p2], p2));
      used_edges.insert(make_pair(p2, p));
      used_edges.insert(make_pair(p, p2));
      if (next[p2] != -1) {
        used_edges.insert(make_pair(p, next[p2]));
        used_edges.insert(make_pair(next[p2], p));
      }

      next[p] = next[p2];
      if (next[p2] != -1)
        prev[next[p2]] = p;
      next[p2] = p;
      prev[p] = p2;
    }
  }
  path.clear();
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  return d;
}

double LocalOptPathR(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points, int r) {
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < path.size() - 1; i++) {
    next[path[i]] = path[i+1];
    prev[path[i+1]] = path[i];
  }
  int first = path[0];

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  for (int i = 1; i < path.size() - 2; i++) {
    if (i % 20000 == 0) {
      printf("%d: %lf\n", i, d);
    }
    int pr = prev[path[i]];
    int p = path[i];
    int l = path[i];
    unordered_set<int> banned;
    banned.insert(l);
    bool bad = false;
    for (int j = 0; j < r-1; j++) {
      if (l == -1) {bad = true; break;}
      l = next[l];
      banned.insert(l);
    }
    if (bad) continue;
    int nn = next[l];
    if (nn == -1) continue;

    if (used_edges.count(make_pair(pr, nn))) continue;

    int best_v = -1;
    double best_improv = 0;
    int s = GetPointSector(points[path[i]]);
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      if (p2 == p) continue;
      if (next[p2] == p) continue;
      if (p2 == l) continue;
      if (next[p2] == l) continue;
      if (next[p2] == -1) continue;
      if (next[l] == p2) continue;
      if (banned.count(p2)) continue;
      // deleted edges:
      // prev[p] -> p, l -> next[p], p2 -> next[p2]
      // new edges:
      // prev[p] -> next[p], p2 -> p, l -> next[p2]
      if (used_edges.count(make_pair(p2, p))) continue;
      if (used_edges.count(make_pair(l, next[p2]))) continue;
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[l], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[l], points[next[p2]]);
      if (dn < dp && (best_v == -1 || best_improv < dp - dn)) {
        best_v = p2;
        best_improv = dp - dn;
      }
      /*else if (dn < dp*1.1 && rand() % 500 == 0) {
        best_v = p2;
        best_improv = dp - dn;
      }
      else if (dn < dp*1.3 && rand() % 50000 == 0) {
        best_v = p2;
        best_improv = dp - dn;
      }*/
    }
    if (best_v != -1) {
      int p2 = best_v;
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[l], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[l], points[next[p2]]);
      d -= dp; d += dn;
      used_edges.erase(make_pair(pr, p));
      used_edges.erase(make_pair(p, pr));
      used_edges.erase(make_pair(nn, l));
      used_edges.erase(make_pair(l, nn));
      used_edges.insert(make_pair(pr, nn));
      used_edges.insert(make_pair(nn, pr));
      next[pr] = nn;
      prev[nn] = pr;
      used_edges.erase(make_pair(p2, next[p2]));
      used_edges.erase(make_pair(next[p2], p2));
      used_edges.insert(make_pair(p2, p));
      used_edges.insert(make_pair(p, p2));
      used_edges.insert(make_pair(l, next[p2]));
      used_edges.insert(make_pair(next[p2], l));

      next[l] = next[p2];
      prev[next[p2]] = l;
      next[p2] = p;
      prev[p] = p2;
    }
  }
  path.clear();
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  return d;
}

double KillCrossingsBegin(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  vector<int> vpos(path.size());
  for (int i = 0; i < path.size(); i++)
    vpos[path[i]] = i;

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }

  //special case for begin
  {
    int s = GetPointSector(points[path[0]]);
    int b = path[0];
    for (int i = 0; i < sectors[s].size(); i++) {
      int p = sectors[s][i];
      if (p == b || p == path[1]) continue;
      if (vpos[p] >= path.size() -1) continue;
      int np = path[vpos[p]+1];
      // edge swap:
      // p->np to b->np
      // p is new beginning
      double dp = CalcDist(points[p], points[np]);
      double dn = CalcDist(points[b], points[np]);
      if (used_edges.count(make_pair(b, np))) continue;
      if (dn < dp) {
        d -= dp; d += dn;
        used_edges.erase(make_pair(p, np));
        used_edges.erase(make_pair(np, p));

        used_edges.insert(make_pair(b, np));
        used_edges.insert(make_pair(np, b));

        int vnp = vpos[np];
        reverse(path.begin(), path.begin() + vnp);
        for (int k = 0; k < vnp; k++)
          vpos[path[k]] = k;
        break;
      }
    }
  }
  return d;
}

double KillCrossings(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  printf("Kill crossings begin\n");
  vector<int> vpos(path.size());
  for (int i = 0; i < path.size(); i++)
    vpos[path[i]] = i;

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  printf("sectors filled\n");

  vector<double> ups;
  for (int s = 0; s < sectors.size(); s++) {
    if (s % 100 == 0)
      printf("%d %lf %d\n", s, d, sectors[s].size());
    for (int i = 0; i < sectors[s].size(); i++) {
      for (int j = i+2; j < sectors[s].size(); j++) {
        int p = sectors[s][i];
        int p2 = sectors[s][j];
        if (vpos[p] == path.size() - 1 || vpos[p2] == path.size() - 1) continue;
        int np = path[vpos[p]+1];
        int np2 = path[vpos[p2]+1];
        // swap from p-np, p2-np2
        // to p-p2, np-np2
        double dp = CalcDist(points[p], points[np]) + CalcDist(points[p2], points[np2]);
        double dn = CalcDist(points[p], points[p2]) + CalcDist(points[np], points[np2]);
        if (dn < dp) {
          ups.push_back(dp - dn);
        }
        if (used_edges.count(make_pair(p, p2))) continue;
        if (used_edges.count(make_pair(np, np2))) continue;
        if (dn < dp) {
          d -= dp; d += dn;
          used_edges.erase(make_pair(p, np));
          used_edges.erase(make_pair(np, p));
          used_edges.erase(make_pair(p2, np2));
          used_edges.erase(make_pair(np2, p2));

          used_edges.insert(make_pair(p2, p));
          used_edges.insert(make_pair(p, p2));
          used_edges.insert(make_pair(np2, np));
          used_edges.insert(make_pair(np, np2));

          int vnp = vpos[np];
          int vnp2 = vpos[np2];
          if (vpos[np] < vpos[np2]) {
            reverse(path.begin() + vnp, path.begin() + vnp2);
            for (int k = vnp; k < vnp2; k++)
              vpos[path[k]] = k;
          } else {
            reverse(path.begin() + vnp2, path.begin() + vnp);
            for (int k = vnp2; k < vnp; k++)
              vpos[path[k]] = k;
          }
        }
      }
    }
  }
  PrintSeriesStats(ups);
  return d;
}


void KillCrossingsHarder(
    vector<int>& path, double& d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    vector<int>& path2, double& d2) {
  vector<int> vpos(path.size()), vpos2(path.size());
  for (int i = 0; i < path.size(); i++) {
    vpos[path[i]] = i;
    vpos2[path2[i]] = i;
  }
  double limit = max(kill_limit, (d - d2)/2);

  unordered_set<pair<int, int> > used1, used2;
  for (int i = 1; i < path.size(); i++) {
    used1.insert(make_pair(path[i-1], path[i]));
    used1.insert(make_pair(path[i], path[i-1]));
    used2.insert(make_pair(path2[i-1], path2[i]));
    used2.insert(make_pair(path2[i], path2[i-1]));
  }

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }

  vector<double> ups;
  for (int s = 0; s < sectors.size(); s++) {
//    if (s % 500 == 0)
//      printf("%d %lf %lf %d\n", s, d, d2, sectors[s].size());
    for (int i = 0; i < sectors[s].size(); i++) {
      for (int j = i+2; j < sectors[s].size(); j++) {
        int p = sectors[s][i];
        int p2 = sectors[s][j];
        if (vpos[p] == path.size() - 1 || vpos[p2] == path.size() - 1) continue;
        int np = path[vpos[p]+1];
        int np2 = path[vpos[p2]+1];
        // swap from p-np, p2-np2
        // to p-p2, np-np2
        double dp = CalcDist(points[p], points[np]) + CalcDist(points[p2], points[np2]);
        double dn = CalcDist(points[p], points[p2]) + CalcDist(points[np], points[np2]);
        if (dp < dn) {
          continue;
        }
        if (used2.count(make_pair(p, p2)) && used2.count(make_pair(np, np2))) continue;
        double d2b = d2;
        if (used2.count(make_pair(p, p2))) {
          if(!KillEdge(path2, used2, vpos2, d2, used1, p, p2, points, min(limit, (dp - dn) *
                  kill_coef))) continue;
        }
        if (used2.count(make_pair(np, np2))) {
          if(!KillEdge(path2, used2, vpos2, d2, used1, np, np2, points, min(limit, (dp - dn) *
                  kill_coef))) continue;
        }
        limit -= d2 - d2b;
        if (dn < dp) {
          d -= dp; d += dn;
          used1.erase(make_pair(p, np));
          used1.erase(make_pair(np, p));
          used1.erase(make_pair(p2, np2));
          used1.erase(make_pair(np2, p2));

          used1.insert(make_pair(p2, p));
          used1.insert(make_pair(p, p2));
          used1.insert(make_pair(np2, np));
          used1.insert(make_pair(np, np2));

          int vnp = vpos[np];
          int vnp2 = vpos[np2];
          if (vpos[np] < vpos[np2]) {
            reverse(path.begin() + vnp, path.begin() + vnp2);
            for (int k = vnp; k < vnp2; k++)
              vpos[path[k]] = k;
          } else {
            reverse(path.begin() + vnp2, path.begin() + vnp);
            for (int k = vnp2; k < vnp; k++)
              vpos[path[k]] = k;
          }
        }
      }
    }
  }

  used_edges.clear();
  for (int i = 1; i < path.size(); i++) {
    used_edges.insert(make_pair(path[i-1], path[i]));
    used_edges.insert(make_pair(path[i], path[i-1]));
    used_edges.insert(make_pair(path2[i-1], path2[i]));
    used_edges.insert(make_pair(path2[i], path2[i-1]));
  }
}

void LocalOptPathRHarder(
    vector<int>& path, double& d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    vector<int>& path2, double& d2, int r) {
  vector<int> vpos(path.size()), vpos2(path.size());
  for (int i = 0; i < path.size(); i++) {
    vpos[path[i]] = i;
    vpos2[path2[i]] = i;
  }
  double limit = max(kill_limit, (d - d2)/2);

  unordered_set<pair<int, int> > used1, used2;
  for (int i = 1; i < path.size(); i++) {
    used1.insert(make_pair(path[i-1], path[i]));
    used1.insert(make_pair(path[i], path[i-1]));
    used2.insert(make_pair(path2[i-1], path2[i]));
    used2.insert(make_pair(path2[i], path2[i-1]));
  }
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < path.size() - 1; i++) {
    next[path[i]] = path[i+1];
    prev[path[i+1]] = path[i];
  }
  int first = path[0];

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  for (int i = 1; i < path.size() - 2; i++) {
/*    if (i % 20000 == 0) {
      printf("%d: %lf %lf\n", i, d, d2);
    }*/
    int pr = prev[path[i]];
    int p = path[i];
    int l = path[i];
    unordered_set<int> banned;
    banned.insert(l);
    bool bad = false;
    for (int j = 0; j < r-1; j++) {
      if (l == -1) {bad = true; break;}
      l = next[l];
      banned.insert(l);
    }
    if (bad) continue;
    int nn = next[l];
    if (nn == -1) continue;

//    if (used1.count(make_pair(pr, nn))) continue;

    int best_v = -1;
    double best_improv = 0;
    int s = GetPointSector(points[path[i]]);
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      if (p2 == p) continue;
      if (next[p2] == p) continue;
      if (next[l] == p2) continue;
      if (next[p2] == -1) continue;
      if (next[next[p2]] == p) continue;
      if (banned.count(p2)) continue;
      // deleted edges:
      // prev[p] -> p, p -> next[p], p2 -> next[p2]
      // new edges:
      // prev[p] -> next[p], p2 -> p, p -> next[p2]
      int bc = 0;
      if (used2.count(make_pair(pr, nn))) bc++;
      if (used2.count(make_pair(p2, p))) bc++;
      if (used2.count(make_pair(l, next[p2]))) bc++;
      if (bc > 1) continue;

      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[l], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[l], points[next[p2]]);
      if (dn < dp && (best_v == -1 || best_improv < dp - dn)) {
        best_v = p2;
        best_improv = dp - dn;
      }
    }
    if (best_v != -1) {
      int p2 = best_v;
      bool bad = false;
      used1.insert(make_pair(pr, nn));
      used1.insert(make_pair(pr, nn));
      used1.insert(make_pair(p2, p));
      used1.insert(make_pair(p, p2));
      used1.insert(make_pair(l, next[p2]));
      used1.insert(make_pair(next[p2], l));
      double d2b = d2;
      if (used2.count(make_pair(pr, nn))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, pr, nn, points, min(limit, best_improv * kill_coef)))
          bad = true;         
      } 
      else if (used2.count(make_pair(p2, p))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, p2, p, points, min(limit, best_improv *
                kill_coef)))
          bad = true;       
      }
      else if (used2.count(make_pair(l, next[p2]))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, l, next[p2], points, min(limit, best_improv *
                kill_coef)))
          bad = true;       
      }
      limit -= d2 - d2b;
      if (bad) {
        used1.erase(make_pair(pr, nn));
        used1.erase(make_pair(pr, nn));
        used1.erase(make_pair(p2, p));
        used1.erase(make_pair(p, p2));
        used1.erase(make_pair(l, next[p2]));
        used1.erase(make_pair(next[p2], l));
        continue;
      }
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[l], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[l], points[next[p2]]);
      d -= dp; d += dn;
      used_edges.erase(make_pair(pr, p));
      used_edges.erase(make_pair(p, pr));
      used_edges.erase(make_pair(nn, l));
      used_edges.erase(make_pair(l, nn));
      next[pr] = nn;
      prev[nn] = pr;
      used1.erase(make_pair(p2, next[p2]));
      used1.erase(make_pair(next[p2], p2));

      next[l] = next[p2];
      prev[next[p2]] = l;
      next[p2] = p;
      prev[p] = p2;
    }
  }
  path.clear();
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  used_edges.clear();
  for (int i = 1; i < path.size(); i++) {
    used_edges.insert(make_pair(path[i-1], path[i]));
    used_edges.insert(make_pair(path[i], path[i-1]));
    used_edges.insert(make_pair(path2[i-1], path2[i]));
    used_edges.insert(make_pair(path2[i], path2[i-1]));
  }
}

void LocalOptPathHarder(
    vector<int>& path, double& d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    vector<int>& path2, double& d2) {
  vector<int> vpos(path.size()), vpos2(path.size());
  for (int i = 0; i < path.size(); i++) {
    vpos[path[i]] = i;
    vpos2[path2[i]] = i;
  }

  unordered_set<pair<int, int> > used1, used2;
  for (int i = 1; i < path.size(); i++) {
    used1.insert(make_pair(path[i-1], path[i]));
    used1.insert(make_pair(path[i], path[i-1]));
    used2.insert(make_pair(path2[i-1], path2[i]));
    used2.insert(make_pair(path2[i], path2[i-1]));
  }
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < path.size() - 1; i++) {
    next[path[i]] = path[i+1];
    prev[path[i+1]] = path[i];
  }
  int first = path[0];

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  for (int i = 1; i < path.size() - 2; i++) {
/*    if (i % 20000 == 0) {
      printf("%d: %lf %lf\n", i, d, d2);
    }*/
    int pr = prev[path[i]];
    int nn = next[path[i]];
    int p = path[i];

//    if (used1.count(make_pair(pr, nn))) continue;

    int best_v = -1;
    double best_improv = 0;
    int s = GetPointSector(points[path[i]]);
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      if (p2 == p) continue;
      if (next[p2] == p) continue;
      if (next[p] == p2) continue;
      if (next[p2] == -1) continue;
      if (next[next[p2]] == p) continue;
      // deleted edges:
      // prev[p] -> p, p -> next[p], p2 -> next[p2]
      // new edges:
      // prev[p] -> next[p], p2 -> p, p -> next[p2]
      int bc = 0;
      if (used2.count(make_pair(pr, nn))) bc++;
      if (used2.count(make_pair(p2, p))) bc++;
      if (used2.count(make_pair(p, next[p2]))) bc++;
      if (bc > 1) continue;

      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[p], points[next[p2]]);
      if (dn < dp && (best_v == -1 || best_improv < dp - dn)) {
        best_v = p2;
        best_improv = dp - dn;
      }
    }
    if (best_v != -1) {
      int p2 = best_v;
      bool bad = false;
      used1.insert(make_pair(prev[p], next[p]));
      used1.insert(make_pair(next[p], prev[p]));
      used1.insert(make_pair(p2, p));
      used1.insert(make_pair(p, p2));
      used1.insert(make_pair(p, next[p2]));
      used1.insert(make_pair(next[p2], p));
      if (used2.count(make_pair(pr, nn))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, pr, nn, points, best_improv * kill_coef))
          bad = true;         
      } 
      else if (used2.count(make_pair(p2, p))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, p2, p, points, best_improv * kill_coef))
          bad = true;       
      }
      else if (used2.count(make_pair(p, next[p2]))) {
        if(!KillEdge(path2, used2, vpos2, d2, used1, p, next[p2], points, best_improv * kill_coef))
          bad = true;       
      }
      if (bad) {
        used1.erase(make_pair(prev[p], next[p]));
        used1.erase(make_pair(next[p], prev[p]));
        used1.erase(make_pair(p2, p));
        used1.erase(make_pair(p, p2));
        used1.erase(make_pair(p, next[p2]));
        used1.erase(make_pair(next[p2], p));
        continue;
      }
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn])
        + CalcDist(points[p2], points[next[p2]]);
      double dn = CalcDist(points[pr], points[nn]) + CalcDist(points[p2], points[p])
        + CalcDist(points[p], points[next[p2]]);
      d -= dp; d += dn;
      used1.erase(make_pair(prev[p], p));
      used1.erase(make_pair(p, prev[p]));
      used1.erase(make_pair(next[p], p));
      used1.erase(make_pair(p, next[p]));
      next[prev[p]] = next[p];
      prev[next[p]] = prev[p];
      used1.erase(make_pair(p2, next[p2]));
      used1.erase(make_pair(next[p2], p2));

      next[p] = next[p2];
      prev[next[p2]] = p;
      next[p2] = p;
      prev[p] = p2;
    }
  }
  path.clear();
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  used_edges.clear();
  for (int i = 1; i < path.size(); i++) {
    used_edges.insert(make_pair(path[i-1], path[i]));
    used_edges.insert(make_pair(path[i], path[i-1]));
    used_edges.insert(make_pair(path2[i-1], path2[i]));
    used_edges.insert(make_pair(path2[i], path2[i-1]));
  }
}

double LocalPerms(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    int perm_size, int jump_size) {
  int rem = rand()%jump_size;
  for (int i = rem; i + perm_size + 1 < path.size(); i+=jump_size) {
    if (i % 10000 == rem)
      printf("perm %d %lf\n", i, d);
    vector<int> perm;
    for (int j = 0; j < perm_size; j++) perm.push_back(j+1);
    vector<int> best_perm;
    double best_dist = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.erase(make_pair(path[i+j], path[i+j+1]));
      used_edges.erase(make_pair(path[i+j+1], path[i+j]));
    }
    double d0 = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      d0 += CalcDist(points[path[i+j]], points[path[i+j+1]]);
    }
    int ii = 0;
    do {
      ++ii;
      vector<int> pp;
      pp.push_back(path[i]);
      for (int j = 0; j < perm_size; j++) {
        pp.push_back(path[i+perm[j]]);
      }
      pp.push_back(path[i+perm_size+1]);
      bool bad = false;
      double d = 0;
      for (int j = 1; j < pp.size(); j++) {
        if (used_edges.count(make_pair(pp[j-1], pp[j]))) bad = true;
        d += CalcDist(points[pp[j-1]], points[pp[j]]);
      }
      if (!bad) {
        if (best_perm.empty() || d < best_dist) {
          best_dist = d;
          best_perm = perm;
        }
      }
//      if (i > 1 && d > d0 && d < d0*1.7) printf("pp %lf %lf\n", d, d/d0);
    } while (next_permutation(perm.begin(), perm.end()));

    d -= d0;
    d += best_dist;
    vector<int> p2;
    for (int j = 0; j < best_perm.size(); j++) {
      p2.push_back(path[i+best_perm[j]]);
    }
    for (int j = 0; j < best_perm.size(); j++)
      path[i+j+1] = p2[j];
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.insert(make_pair(path[i+j], path[i+j+1]));
      used_edges.insert(make_pair(path[i+j+1], path[i+j]));
    }
  }
  return d;
}

double SpanningTree(vector<pair<double, double> >& pp) {
  vector<int> col(pp.size());
  int n = pp.size();
  for (int i = 0; i < col.size(); i++)
    col[i] = i;
  vector<pair<double, pair<int, int> > > h;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < i; j++) {
      double d = CalcDist(pp[i], pp[j]);
      h.push_back(make_pair(d, make_pair(i, j)));
    }
  }
  sort(h.begin(), h.end());
  double sp = 0;
  for (int i = 0; i < h.size(); i++) {
    int v1 = h[i].second.first;
    int v2 = h[i].second.second;
    int c1 = col[v1];
    int c2 = col[v2];
    if (c1 == c2) continue;
    for (int j = 0; j < n; j++) {
      if (col[j] == c1)
        col[j] = c2;
    }
    sp += h[i].first;
  }
  return sp;
}

double LocalPermsB(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    int perm_size) {
  vector<pair<pair<double, double>, int> > best_imp;
  for (int i = 0; i + perm_size + 47 < path.size(); i+=perm_size+3) {
    double dd = 0;
    // i (i+1, ..., i + perm_size) i + perm_size + 1
    for (int j = 1; j < perm_size + 2; j++) {
      dd += CalcDist(points[path[i+j]], points[path[i+j-1]]);
    }
    vector<pair<double, double> > pp;
    for (int j = 0; j < perm_size + 2; j++)
      pp.push_back(points[path[i+j]]);
    double dist = dd;
    dd -= SpanningTree(pp);
    best_imp.push_back(make_pair(make_pair(dd, dist), i));
  }
  sort(best_imp.rbegin(), best_imp.rend());
  for (int k = 0; k < 20; k++) {
    int i = best_imp[k].second;
    vector<int> perm;
    for (int j = 0; j < perm_size; j++) perm.push_back(j+1);
    vector<int> best_perm;
    double best_dist = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.erase(make_pair(path[i+j], path[i+j+1]));
      used_edges.erase(make_pair(path[i+j+1], path[i+j]));
    }
    double d0 = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      d0 += CalcDist(points[path[i+j]], points[path[i+j+1]]);
    }
    int ii = 0;
    do {
      ++ii;
      vector<int> pp;
      pp.push_back(path[i]);
      for (int j = 0; j < perm_size; j++) {
        pp.push_back(path[i+perm[j]]);
      }
      pp.push_back(path[i+perm_size+1]);
      bool bad = false;
      double d = 0;
      for (int j = 1; j < pp.size(); j++) {
        if (used_edges.count(make_pair(pp[j-1], pp[j]))) bad = true;
        d += CalcDist(points[pp[j-1]], points[pp[j]]);
      }
      if (!bad) {
        if (best_perm.empty() || d < best_dist) {
          best_dist = d;
          best_perm = perm;
        }
      }
//      if (i > 1 && d > d0 && d < d0*1.7) printf("pp %lf %lf\n", d, d/d0);
    } while (next_permutation(perm.begin(), perm.end()));

    d -= d0;
    d += best_dist;
    vector<int> p2;
    for (int j = 0; j < best_perm.size(); j++) {
      p2.push_back(path[i+best_perm[j]]);
    }
    for (int j = 0; j < best_perm.size(); j++)
      path[i+j+1] = p2[j];
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.insert(make_pair(path[i+j], path[i+j+1]));
      used_edges.insert(make_pair(path[i+j+1], path[i+j]));
    }
  }
  return d;
}

void SwapCommonParts(vector<int>& p1, vector<int>& p2, double& d1, double& d2,
    unordered_map<int, pair<double, double> >& points) {
  int seq_size = 5;
  map<set<int>, int> seqs;
  for (int i = 1; i + seq_size - 1 < p1.size(); i++) {
    set<int> s;
    for (int j = 0; j < seq_size; j++) {
      s.insert(p1[i+j]);
    }
    seqs[s] = i;
  }
  for (int i = 1; i + seq_size - 1 < p1.size(); i++) {
    set<int> s;
    for (int j = 0; j < seq_size; j++) {
      s.insert(p2[i+j]);
    }
    if (seqs.count(s)) {
      int pp = seqs[s];
      set<int> s2;
      for (int j = 0; j < seq_size; j++) {
        s2.insert(p1[pp+j]);
      }
      if (s2 == s && rand() % 10 == 0) {
        for (int j = 0; j < seq_size; j++) {
          swap(p1[pp+j], p2[i+j]);
        }
      }
    }
  }
  d1 = 0, d2 = 0;
  for (int i = 1; i < p1.size(); i++) {
    d1 += CalcDist(points[p1[i-1]], points[p1[i]]);
    d2 += CalcDist(points[p2[i-1]], points[p2[i]]);
  }
}

double LocalPermsD(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    int perm_size) {
  vector<pair<pair<double, double>, int> > best_imp;
  for (int i = 0; i + perm_size + 47 < path.size(); i+=perm_size+3) {
    double dd = 0;
    // i (i+1, ..., i + perm_size) i + perm_size + 1
    vector<pair<double, double> > pp;
    for (int j = 0; j < perm_size + 2; j++)
      pp.push_back(points[path[i+j]]);
    for (int j = 1; j < perm_size + 2; j++) {
      dd += CalcDist(points[path[i+j]], points[path[i+j-1]]);
    }
    double dist = dd;
    dd -= SpanningTree(pp);
    best_imp.push_back(make_pair(make_pair(dd, dist), i));
  }
  sort(best_imp.rbegin(), best_imp.rend());
  double threshold = best_imp[best_imp.size()/4].first.first;
  for (int i = 0; i < 100; i++) {
    int start = best_imp[i].second;
//    int start = i;
    vector<pair<double, double> > ppx;
    for (int j = 0; j < perm_size + 2; j++)
      ppx.push_back(points[path[start+j]]);
    double d0 = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      d0 += CalcDist(points[path[start+j]], points[path[start+j+1]]);
    }
    double estimate = d0 - SpanningTree(ppx);
    if (estimate < threshold)
      continue;
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.erase(make_pair(path[start+j], path[start+j+1]));
      used_edges.erase(make_pair(path[start+j+1], path[start+j]));
    }
    vector<int> pp;
    for (int j = 0; j < perm_size; j++)
      pp.push_back(path[start+j+1]);
    int end = start + perm_size + 1;
    vector<vector<double> > best(1 << perm_size, vector<double>(perm_size, 1000000000));
    vector<vector<int> > from(1 << perm_size, vector<int>(perm_size, -1));
    for  (int i = 1; i < 1 << perm_size; i++) {
      for (int j = 0; j < perm_size; j++) {
        if (((1 << j) & i) == 0) continue;
        int pr = (i ^ (1 << j));
        if (pr == 0) {
          if (used_edges.count(make_pair(path[start], pp[j]))) continue;
          double dd = CalcDist(points[path[start]], points[pp[j]]);
          best[i][j] = min(best[i][j], dd);
        } else {
          for (int k = 0; k < perm_size; k++) {
            if (((1 << k) & pr) == 0) continue;
            if (used_edges.count(make_pair(pp[k], pp[j]))) continue;
            double dd = CalcDist(points[pp[k]], points[pp[j]]);
            if (best[pr][k] + dd < best[i][j]) {
              best[i][j] = best[pr][k] + dd;
              from[i][j] = k;
            }
          }
        }
      }
    }

    double bestl = 1000000000;
    int bestf = 0;
    int ll = (1 << perm_size) - 1;
    for (int j = 0; j < perm_size; j++) {
      if (used_edges.count(make_pair(pp[j], path[end]))) continue;
      double dd = CalcDist(points[pp[j]], points[path[end]]);
      if (bestl > best[ll][j] + dd) {
        bestl = best[ll][j] + dd;
        bestf = j;
      }
    }
    vector<int> pp2;
    int cur = ll;
    while (cur != 0) {
      pp2.push_back(pp[bestf]);
      int pr = cur ^ (1 << bestf);
      bestf = from[cur][bestf];
      cur = pr;
    }
    reverse(pp2.begin(), pp2.end());
    d -= d0;
    d += bestl;
    for (int j = 0; j < perm_size; j++)
      path[start+j+1] = pp2[j];
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.insert(make_pair(path[start+j], path[start+j+1]));
      used_edges.insert(make_pair(path[start+j+1], path[start+j]));
    }
    if (i % 5000 == 0)
      printf("pd %d %lf\n", i, d);
  }
  return d;
}

double LocalPermsC(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    int perm_size) {
  vector<pair<pair<double, double>, int> > best_imp;
  for (int i = 0; i + perm_size + 47 < path.size(); i+=perm_size+3) {
    double dd = 0;
    // i (i+1, ..., i + perm_size) i + perm_size + 1
    vector<pair<double, double> > pp;
    for (int j = 0; j < perm_size + 2; j++)
      pp.push_back(points[path[i+j]]);
    for (int j = 1; j < perm_size + 2; j++) {
      dd += CalcDist(points[path[i+j]], points[path[i+j-1]]);
    }
    double dist = dd;
    dd -= SpanningTree(pp);
    best_imp.push_back(make_pair(make_pair(dd, dist), i));
  }
  sort(best_imp.rbegin(), best_imp.rend());
  double threshold = best_imp[best_imp.size()/10].first.first;
//  for (int k = 0; k < 20; k++) {
  for (int i = 1; i + perm_size + 4 < path.size(); i++) {
    if (i % 10000 == 0)
      printf("pd %d %lf\n", i, d);
//    int start = best_imp[k].second;
    int start = i;
    vector<pair<double, double> > ppx;
    for (int j = 0; j < perm_size + 2; j++)
      ppx.push_back(points[path[start+j]]);
    double d0 = 0;
    for (int j = 0; j < perm_size + 1; j++) {
      d0 += CalcDist(points[path[start+j]], points[path[start+j+1]]);
    }
    double estimate = d0 - SpanningTree(ppx);
    if (estimate < threshold)
      continue;
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.erase(make_pair(path[start+j], path[start+j+1]));
      used_edges.erase(make_pair(path[start+j+1], path[start+j]));
    }
    vector<int> pp;
    for (int j = 0; j < perm_size; j++)
      pp.push_back(path[start+j+1]);
    int end = start + perm_size + 1;
    vector<vector<double> > best(1 << perm_size, vector<double>(perm_size, 1000000000));
    vector<vector<int> > from(1 << perm_size, vector<int>(perm_size, -1));
    for  (int i = 1; i < 1 << perm_size; i++) {
      for (int j = 0; j < perm_size; j++) {
        if (((1 << j) & i) == 0) continue;
        int pr = (i ^ (1 << j));
        if (pr == 0) {
          if (used_edges.count(make_pair(path[start], pp[j]))) continue;
          double dd = CalcDist(points[path[start]], points[pp[j]]);
          best[i][j] = min(best[i][j], dd);
        } else {
          for (int k = 0; k < perm_size; k++) {
            if (((1 << k) & pr) == 0) continue;
            if (used_edges.count(make_pair(pp[k], pp[j]))) continue;
            double dd = CalcDist(points[pp[k]], points[pp[j]]);
            if (best[pr][k] + dd < best[i][j]) {
              best[i][j] = best[pr][k] + dd;
              from[i][j] = k;
            }
          }
        }
      }
    }

    double bestl = 1000000000;
    int bestf = 0;
    int ll = (1 << perm_size) - 1;
    for (int j = 0; j < perm_size; j++) {
      if (used_edges.count(make_pair(pp[j], path[end]))) continue;
      double dd = CalcDist(points[pp[j]], points[path[end]]);
      if (bestl > best[ll][j] + dd) {
        bestl = best[ll][j] + dd;
        bestf = j;
      }
    }
    vector<int> pp2;
    int cur = ll;
    while (cur != 0) {
      pp2.push_back(pp[bestf]);
      int pr = cur ^ (1 << bestf);
      bestf = from[cur][bestf];
      cur = pr;
    }
    reverse(pp2.begin(), pp2.end());
    d -= d0;
    d += bestl;
    for (int j = 0; j < perm_size; j++)
      path[start+j+1] = pp2[j];
    for (int j = 0; j < perm_size + 1; j++) {
      used_edges.insert(make_pair(path[start+j], path[start+j+1]));
      used_edges.insert(make_pair(path[start+j+1], path[start+j]));
    }
  }
  return d;
}
double UglyStuff(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points,
    vector<int>& path2) {
  printf("ugly begin %lf\n", d);
  unordered_set<pair<int, int> > used1, used2;
  for (int i = 1; i < path.size(); i++) {
    used1.insert(make_pair(path[i-1], path[i]));
    used1.insert(make_pair(path[i], path[i-1]));
    used2.insert(make_pair(path2[i-1], path2[i]));
    used2.insert(make_pair(path2[i], path2[i-1]));
  }
  vector<int> pathx;
  vector<int> free;
  for (int i = 0; i < path.size(); i++) {
    if (i + 4 >= path.size()) {
      pathx.push_back(path[i]);
      continue;
    }
    if (used2.count(make_pair(path[i], path[i+4]))) {
      pathx.push_back(path[i]);
      continue;
    }
    if (rand() % 10 == 0) {
      pathx.push_back(path[i]);
      free.push_back(path[i+1]);
      free.push_back(path[i+2]);
      free.push_back(path[i+3]);
      pathx.push_back(path[i+4]);
      i += 4;
    } else {
      pathx.push_back(path[i]);
      continue;
    }
  }
  printf("free size %d pathx %d\n", free.size(), pathx.size());
  LocalPerms(pathx, d, used2, points, 3, 1);
  used2.clear();
  for (int i = 1; i < path.size(); i++) {
    used2.insert(make_pair(path2[i-1], path2[i]));
    used2.insert(make_pair(path2[i], path2[i-1]));
  }
  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < pathx.size(); i++) {
    int s = GetPointSector(points[pathx[i]]);
    sectors[s].push_back(pathx[i]);
  }
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < pathx.size() - 1; i++) {
    next[pathx[i]] = pathx[i+1];
    prev[pathx[i+1]] = pathx[i];
  }
  int first = pathx[0];
  for (int i = 0; i < free.size(); i++) {
    int p = free[i];
    int s = GetPointSector(points[p]);
    double best_d = 0;
    int best_v = -1;
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      int np2 = next[p2];
      if (np2 == -1) continue;
      if (used2.count(make_pair(p2, p))) continue;
      if (used2.count(make_pair(p, np2))) continue;
      double dd = CalcDist(points[p2], points[p]) + CalcDist(points[p], points[np2]);
      if (dd < best_d || best_v == -1) {
        best_d = dd; best_v = p2;
      }
    }
    if (best_v == -1) {
      for (int s = 0; s < sectors.size(); s++) {
        for (int j = 0; j < sectors[s].size(); j++) {
          int p2 = sectors[s][j];
          int np2 = next[p2];
          if (np2 == -1) continue;
          if (used2.count(make_pair(p2, p))) continue;
          if (used2.count(make_pair(p, np2))) continue;
          double dd = CalcDist(points[p2], points[p]) + CalcDist(points[p], points[np2]);
          if (dd < best_d || best_v == -1) {
            best_d = dd; best_v = p2;
          }
        }
      }
    }
    if (best_v == -1) {
      printf("fuuuuuuuuuck %d\n", i);
      return d;
    }
    int p2 = best_v;
    int np2 = next[p2];
    next[p] = np2;
    prev[np2] = p;
    prev[p] = p2;
    next[p2] = p;
  }
  path.clear();
  d = 0;
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  used_edges.clear();
  for (int i = 1; i < path.size(); i++) {
    used_edges.insert(make_pair(path[i-1], path[i]));
    used_edges.insert(make_pair(path[i], path[i-1]));
    used_edges.insert(make_pair(path2[i-1], path2[i]));
    used_edges.insert(make_pair(path2[i], path2[i-1]));
    d += CalcDist(points[path[i]], points[path[i-1]]);
  }
  printf("ugly end %lf\n", d);
  return d;
}

int TestPath(const vector<int>& p1, const vector<int>& p2,
             unordered_map<int, pair<double, double> >& points) {
  if (p1.size() != points.size()) {
    printf("too few points");
    return 1;
  }
  unordered_set<pair<int, int> > bad_edges;
  double dist1 = 0, dist2 = 0;
  int previd = p1[0];
  unordered_set<int> visited;
  visited.insert(previd);
  for (int i = 1; i < p1.size(); i++) {
    int id = p1[i];
    visited.insert(id);
    if (points.count(id) == 0 || points.count(previd) == 0) {
      printf("unknown point\n");
      return 1;
    }
    dist1 += CalcDist(points[id], points[previd]);
    bad_edges.insert(make_pair(id, previd));
    bad_edges.insert(make_pair(previd, id));
    previd = id;
  }
  if (visited.size() != points.size()) {
    printf("didn't visit all points on path 1\n");
    return 1;
  }
  printf("path 1 length: %lf\n", dist1);

  if (p2.size() != points.size()) {
    printf("too few points");
    return 1;
  }
  previd = p2[0];

  visited.clear();
  visited.insert(previd);
  for (int i = 1; i < points.size(); i++) {
    int id = p2[i];
    visited.insert(id);
    if (points.count(id) == 0 || points.count(previd) == 0) {
      printf("unknown point\n");
      return 1;
    }
    dist2 += CalcDist(points[id], points[previd]);
    if (bad_edges.count(make_pair(id, previd))) {
      printf("double edge %d %d\n", id, previd);
      return 1;
    }
    previd = id;
  }
  if (visited.size() != points.size()) {
    printf("didn't visit all points on path 2\n");
    return 1;
  }
  printf("path 2 length: %lf\n", dist2);
  printf("score: %lf\n", max(dist1, dist2));
  return 0;
}

void ShufflePoints(vector<int>& path, 
             unordered_map<int, pair<double, double> >& points) {
  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(i);
  }
  int s1 = rand()%sectors.size();
  int s2 = rand()%sectors.size();
  while(s1 == s2) s2 = rand()%sectors.size();
  int a = sectors[s1][rand()%sectors[s1].size()];
  int b = sectors[s2][rand()%sectors[s2].size()];
  int mi = min(a,b);
  int ma = max(a,b);
  reverse(path.begin() + mi, path.begin() + ma);
}

void ShufflePath(vector<int> &path, unordered_set<pair<int, int> >& used_edges) {
  for (int i = 0; i + 200 < path.size(); i+=500) {
    int b = i + rand()%70+10;
    int c = b + 5;
    int d = b + 10;
    // (b+1...c) (c+1...d)
    // cesty: b->b+1, c->c+1, d->d+1
    // b->c+1, d->b+1, c->d+1
    if (used_edges.count(make_pair(path[b], path[c+1]))) continue;
    if (used_edges.count(make_pair(path[d], path[b+1]))) continue;
    if (used_edges.count(make_pair(path[c], path[d+1]))) continue;
    used_edges.erase(make_pair(path[b], path[b+1]));
    used_edges.erase(make_pair(path[b+1], path[b]));
    used_edges.erase(make_pair(path[c], path[c+1]));
    used_edges.erase(make_pair(path[c+1], path[c]));
    used_edges.erase(make_pair(path[d], path[d+1]));
    used_edges.erase(make_pair(path[d+1], path[d]));

    used_edges.insert(make_pair(path[b], path[c+1]));
    used_edges.insert(make_pair(path[c+1], path[b]));
    used_edges.insert(make_pair(path[d], path[b+1]));
    used_edges.insert(make_pair(path[b+1], path[d]));
    used_edges.insert(make_pair(path[c], path[d+1]));
    used_edges.insert(make_pair(path[d+1], path[c]));

    reverse(path.begin() + b + 1, path.begin() + d + 1);
    reverse(path.begin() + b + 1, path.begin() + c + 1);
    reverse(path.begin() + c + 1, path.begin() + d + 1);
  }
}

void LocalOptPaths(vector<int>&p1, double&d1, vector<int>& p2, double&d2,
    unordered_set<pair<int, int> >&used_edges, unordered_map<int, pair<double, double> >& points) {
  unordered_set<pair<int, int> > used1, used2;
  for (int i = 1; i < p1.size(); i++) {
    used1.insert(make_pair(p1[i-1], p1[i]));
    used1.insert(make_pair(p1[i], p1[i-1]));
    used2.insert(make_pair(p2[i-1], p2[i]));
    used2.insert(make_pair(p2[i], p2[i-1]));
  }

  vector<int> next1(p1.size(), -1); 
  vector<int> prev1(p1.size(), -1); 
  vector<int> next2(p1.size(), -1); 
  vector<int> prev2(p1.size(), -1); 
  for (int i = 0; i < p1.size() - 1; i++) {
    next1[p1[i]] = p1[i+1];
    prev1[p1[i+1]] = p1[i];
  }
  for (int i = 0; i < p1.size() - 1; i++) {
    next2[p2[i]] = p2[i+1];
    prev2[p2[i+1]] = p2[i];
  }
  int first1 = p1[0];
  int first2 = p2[0];
  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < p1.size(); i++) {
    int s = GetPointSector(points[p1[i]]);
    sectors[s].push_back(p1[i]);
  }

  for (int i = 0; i < p1.size(); i++) {
    if (i % 100 == 0) printf("%d\n", i);
    if (i == first1) continue;
    if (i == first2) continue;
    if (next2[i] == -1) continue;
    if (next1[i] == -1) continue;

    int pr1 = prev1[i];
    int nn1 = next1[i];
    int pr2 = prev2[i];
    int nn2 = next2[i];
    if (used1.count(make_pair(pr2, nn2))) continue;
    if (used2.count(make_pair(pr1, nn1))) continue;

    int s = GetPointSector(points[i]);
    bool done = false;
    for (int j = 0; j < sectors[s].size() && !done; j++) {
      int x1 = sectors[s][j];
      if (x1 == i) continue;
      if (next1[x1] == i) continue;
      if (next1[x1] == -1) continue;
      int nx1 = next1[x1];
      // pr1->i, i -> nn1, x1->nx1 to
      // pr1->nn1, x1->i, i->nx1
      if (used2.count(make_pair(x1, i))) continue;
      if (used2.count(make_pair(i, nx1))) continue;
      double nd1 = d1 - CalcDist(points[pr1], points[i]) - CalcDist(points[i], points[nn1])
          - CalcDist(points[x1], points[nx1]) + CalcDist(points[pr1], points[nn1])
          + CalcDist(points[x1], points[i]) + CalcDist(points[i], points[nx1]);

      for (int k = 0; k < sectors[s].size() && !done; k++) {
        if (k == j) continue;
        int x2 = sectors[s][j];
        if (x2 == i) continue;
        int nx2 = next2[x2];
        if (nx2 == -1) continue;
        if (nx2 == x1) continue;
        if (nx1 == x2) continue;
        // pr2->i, i -> nn2, x2->nx2 to
        // pr2->nn1, x2->i, i->nx2
        if (used1.count(make_pair(x2, i))) continue;
        if (used1.count(make_pair(i, nx2))) continue;
        double nd2 = d2 - CalcDist(points[pr2], points[i]) - CalcDist(points[i], points[nn2])
          - CalcDist(points[x2], points[nx2]) + CalcDist(points[pr2], points[nn2])
          + CalcDist(points[x2], points[i]) + CalcDist(points[i], points[nx2]);
        if (max(nd1, nd2) < max(d1, d2)) {
          d1 = nd1; d2 = nd2; done = true;
          printf("d1 %lf d2 %lf\n", d1, d2);
          used1.erase(make_pair(pr1, i));    used1.erase(make_pair(i, pr1));
          used1.erase(make_pair(nn1, i));    used1.erase(make_pair(i, nn1));
          used1.erase(make_pair(x1, nx1));   used1.erase(make_pair(nx1, x1));
          used1.insert(make_pair(pr1, nn1)); used1.insert(make_pair(nn1, pr1));
          used1.insert(make_pair(x1, i));    used1.insert(make_pair(i, x1));
          used1.insert(make_pair(nx1, i));   used1.insert(make_pair(i, nx1));

          used2.erase(make_pair(pr2, i));    used2.erase(make_pair(i, pr2));
          used2.erase(make_pair(nn2, i));    used2.erase(make_pair(i, nn2));
          used2.erase(make_pair(x2, nx2));   used2.erase(make_pair(nx2, x2));
          used2.insert(make_pair(pr2, nn2)); used2.insert(make_pair(nn2, pr2));
          used2.insert(make_pair(x2, i));    used2.insert(make_pair(i, x2));
          used2.insert(make_pair(nx2, i));   used2.insert(make_pair(i, nx2));

          next1[pr1] = nn1; prev1[nn1] = pr1;
          next1[x1] = i; next1[i] = nx1; prev1[i] = x1; prev1[nx1] = i;

          next2[pr2] = nn2; prev2[nn2] = pr2;
          next2[x2] = i; next2[i] = nx2; prev2[i] = x2; prev2[nx2] = i;
        }
      }
    }
  }
  p1.clear();
  int cur = first1;
  while (cur != -1) {
    p1.push_back(cur);
    cur = next1[cur];
  }
  
  p2.clear();
  cur = first2;
  while (cur != -1) {
    p2.push_back(cur);
    cur = next2[cur];
  }

  used_edges.clear();
  for (int i = 1; i < p1.size(); i++) {
    used_edges.insert(make_pair(p1[i-1], p1[i]));
    used_edges.insert(make_pair(p1[i], p1[i-1]));
    used_edges.insert(make_pair(p2[i-1], p2[i]));
    used_edges.insert(make_pair(p2[i], p2[i-1]));
  }
}

double LocalPermsX(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  int perm_size = 5;
  int run_size = 5;
  int opt_size = perm_size * run_size;
  for (int i = 1; i + perm_size*run_size + 3 < path.size(); i+=1) {
    if (i % 10000 == 0)
      printf("perm %d %lf\n", i, d);
    vector<int> perm;
    for (int j = 0; j < perm_size; j++) perm.push_back(j);
    vector<int> best_perm;
    double best_dist = 0;
    for (int j = 0; j < opt_size + 1; j++) {
      used_edges.erase(make_pair(path[i+j], path[i+j+1]));
      used_edges.erase(make_pair(path[i+j+1], path[i+j]));
    }
    double d0 = 0;
    for (int j = 0; j < opt_size + 1; j++) {
      d0 += CalcDist(points[path[i+j]], points[path[i+j+1]]);
    }
    vector<vector<int> > px(perm_size);
    for (int j = 0; j < perm_size; j++) {
      for (int k = 0; k < run_size; k++)
        px[j].push_back(path[i+1+k+j*run_size]);
    }
    do {
      vector<int> pp;
      pp.push_back(path[i]);
      for (int j = 0; j < perm_size; j++) {
        for (int k = 0; k < run_size; k++) 
          pp.push_back(px[perm[j]][k]);
      }
      pp.push_back(path[i+opt_size+1]);
      bool bad = false;
      double d = 0;
      for (int j = 1; j < pp.size(); j++) {
        if (used_edges.count(make_pair(pp[j-1], pp[j]))) bad = true;
        d += CalcDist(points[pp[j-1]], points[pp[j]]);
      }
      if (!bad) {
        if (best_perm.empty() || d < best_dist) {
          best_dist = d;
          best_perm = perm;
        }
      }
    } while (next_permutation(perm.begin(), perm.end()));

    d -= d0;
    d += best_dist;
    vector<int> p2;
    for (int j = 0; j < best_perm.size(); j++) {
      for (int k = 0; k < run_size; k++) 
        p2.push_back(px[best_perm[j]][k]);
    }
    for (int j = 0; j < p2.size(); j++)
      path[i+j+1] = p2[j];
    for (int j = 0; j < opt_size + 1; j++) {
      used_edges.insert(make_pair(path[i+j], path[i+j+1]));
      used_edges.insert(make_pair(path[i+j+1], path[i+j]));
    }
  }
  return d;
}

double LocalOptPathY(
    vector<int>& path, double d, unordered_set<pair<int, int> >& used_edges,
    unordered_map<int, pair<double, double> >& points) {
  vector<int> next(path.size(), -1);
  vector<int> prev(path.size(), -1);
  for (int i = 0; i < path.size() - 1; i++) {
    next[path[i]] = path[i+1];
    prev[path[i+1]] = path[i];
  }
  int first = path[0];

  vector<vector<int> > sectors(40*40);
  for (int i = 1; i < path.size(); i++) {
    int s = GetPointSector(points[path[i]]);
    sectors[s].push_back(path[i]);
  }
  for (int i = 1; i < path.size() - 2; i++) {
    if (i % 20000 == 0) {
      printf("%d: %lf\n", i, d);
    }
    int pr = prev[path[i]];
    int nn = next[path[i]];
    int p = path[i];

//    if (used_edges.count(make_pair(pr, nn))) continue;

    int best_v = -1;
    double best_improv = 0;
    int s = GetPointSector(points[path[i]]);
    for (int j = 0; j < sectors[s].size(); j++) {
      int p2 = sectors[s][j];
      if (p2 == p) continue;
      if (next[p2] == p) continue;
      if (next[p2] == -1) continue;
      if (prev[p2] == -1) continue;
      if (prev[p2] == p) continue;
      int pr2 = prev[p2];
      int nn2 = next[p2];

      // deleted edges:
      // prev[p] -> p, p -> next[p], prev[p2] -> p2, p2 -> next[p2]
      // new edges:
      // prev[p] -> p2 -> next[p], prev[p2] -> p -> next[p2]
      if (used_edges.count(make_pair(pr, p2))) continue;
      if (used_edges.count(make_pair(p2, nn))) continue;
      if (used_edges.count(make_pair(pr2, p))) continue;
      if (used_edges.count(make_pair(p, nn2))) continue;
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn]);
      dp += CalcDist(points[pr2], points[p2]) + CalcDist(points[p2], points[nn2]);
      double dn = CalcDist(points[pr], points[p2]) + CalcDist(points[p2], points[nn]);
      dn += CalcDist(points[pr2], points[p]) + CalcDist(points[p], points[nn2]);
      if (dn < dp && (best_v == -1 || best_improv < dp - dn)) {
        best_v = p2;
        best_improv = dp - dn;
      }
    }
    if (best_v != -1) {
      int p2 = best_v;
      int pr2 = prev[p2], nn2 = next[p2];
      double dp = CalcDist(points[pr], points[p]) + CalcDist(points[p], points[nn]);
      dp += CalcDist(points[pr2], points[p2]) + CalcDist(points[p2], points[nn2]);
      double dn = CalcDist(points[pr], points[p2]) + CalcDist(points[p2], points[nn]);
      dn += CalcDist(points[pr2], points[p]) + CalcDist(points[p], points[nn2]);
      d -= dp; d += dn;
      used_edges.erase(make_pair(prev[p], p));
      used_edges.erase(make_pair(p, prev[p]));
      used_edges.erase(make_pair(next[p], p));
      used_edges.erase(make_pair(p, next[p]));
      used_edges.erase(make_pair(prev[p2], p2));
      used_edges.erase(make_pair(p2, prev[p2]));
      used_edges.erase(make_pair(next[p2], p2));
      used_edges.erase(make_pair(p2, next[p2]));

      used_edges.insert(make_pair(prev[p], p2));
      used_edges.insert(make_pair(p2, prev[p]));
      used_edges.insert(make_pair(next[p], p2));
      used_edges.insert(make_pair(p2, next[p]));
      used_edges.insert(make_pair(prev[p2], p));
      used_edges.insert(make_pair(p, prev[p2]));
      used_edges.insert(make_pair(next[p2], p));
      used_edges.insert(make_pair(p, next[p2]));

      next[pr] = p2; next[p2] = nn;
      prev[p2] = pr; prev[nn] = p2;

      next[pr2] = p; next[p] = nn2;
      prev[p] = pr2; prev[nn2] = p;
    }
  }
  path.clear();
  int cur = first;
  while (cur != -1) {
    path.push_back(cur);
    cur = next[cur];
  }
  return d;
}

