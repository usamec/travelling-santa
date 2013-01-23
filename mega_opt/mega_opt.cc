#include "utility.h"
#include "path.h"

#include <boost/thread/mutex.hpp>
#include <boost/thread.hpp>

vector<vector<int> > closest;
vector<pair<double, double> > points(150000);

void TestValidMove() {
  vector<int> p;
  for (int i = 0; i < 20; i++)
    p.push_back(i);
  Path pp(p);
  KOptMove move1;
  move1.out.push_back(make_pair(3, 4));
  move1.out.push_back(make_pair(7, 8));
  move1.in.push_back(make_pair(3, 7));
  move1.in.push_back(make_pair(4, 8));
  printf("TRUE %d\n", pp.IsValidMove(move1));
  move1.in[0] = make_pair(3, 8);
  move1.in[1] = make_pair(4, 7);
  printf("FALSE %d\n", pp.IsValidMove(move1));

  KOptMove move2;
  move2.out.push_back(make_pair(0, 1));
  move2.out.push_back(make_pair(7, 8));
  move2.in.push_back(make_pair(0, 7));
  move2.in.push_back(make_pair(1, 8));
  printf("TRUE %d\n", pp.IsValidMove(move2));

  KOptMove move3;
  move3.out.push_back(make_pair(0, 1));
  move3.out.push_back(make_pair(18, 19));
  move3.in.push_back(make_pair(0, 18));
  move3.in.push_back(make_pair(1, 19));
  printf("TRUE %d\n", pp.IsValidMove(move3));

  KOptMove move4;
  move4.out.push_back(make_pair(0, 1));
  move4.out.push_back(make_pair(1, 2));
  move4.out.push_back(make_pair(18, 19));
  move4.in.push_back(make_pair(0, 2));
  move4.in.push_back(make_pair(1, 19));
  move4.in.push_back(make_pair(1, 18));
  printf("TRUE %d\n", pp.IsValidMove(move4));

  KOptMove move5;
  move5.out.push_back(make_pair(4, 5));
  move5.out.push_back(make_pair(7, 8));
  move5.out.push_back(make_pair(10, 11));
  move5.in.push_back(make_pair(4, 11));
  move5.in.push_back(make_pair(5, 8));
  move5.in.push_back(make_pair(7, 10));
  printf("FALSE %d\n", pp.IsValidMove(move5));

  KOptMove move6;
  move6.out.push_back(make_pair(4, 5));
  move6.out.push_back(make_pair(6, 7));
  move6.out.push_back(make_pair(15, 16));
  move6.in.push_back(make_pair(5, 15));
  move6.in.push_back(make_pair(6, 16));
  move6.in.push_back(make_pair(4, 7));
  printf("FALSE %d\n", pp.IsValidMove(move6));
}

void Unstuck(Path& p1, Path& p2, vector<pair<double, double> >& points) {
  static double unstuck_coef = 1.001;
  double d1 = p1.GetLength(points), d2 = p2.GetLength(points);
  printf("unstuck begin %lf %lf %lf\n", d1, d2, (d1+d2)/2);
  int suc = 0;
  int unstucks = 1000 + rand()%600;
  int k = 0;
  for (; suc < unstucks && k < 150000; k++) {
    if (k % 5000 == 0) printf("k %d %d/%d\n", k, suc, unstucks);
    if (d1 > d2) {
      int i = 4 +rand() % (p1.p_.size()-10);
      KOptMove res = p1.GetKOptMove(p2, p1.p_[i], p1.p_[i+1], closest, points, 5, unstuck_coef);
      if (res.valid) {
        p1.ExecuteMove(res);
        for (int j = 0; j < res.in.size(); j++)
          d1 += CalcDist(points[res.in[j].first], points[res.in[j].second]);
        for (int j = 0; j < res.out.size(); j++)
          d1 -= CalcDist(points[res.out[j].first], points[res.out[j].second]);
        suc++;
        if (!res.kick_in.empty()) {
          p2.ExecuteMove(res.kick_in, res.kick_out);
          for (int j = 0; j < res.kick_in.size(); j++)
            d2 += CalcDist(points[res.kick_in[j].first], points[res.kick_in[j].second]);
          for (int j = 0; j < res.kick_out.size(); j++)
            d2 -= CalcDist(points[res.kick_out[j].first], points[res.kick_out[j].second]);
        }
      }
    } else {
      int i = 4 +rand() % (p1.p_.size()-10);
      KOptMove res = p2.GetKOptMove(p1, p2.p_[i], p2.p_[i+1], closest, points, 5, unstuck_coef);
      if (res.valid) {
        p2.ExecuteMove(res);
        for (int j = 0; j < res.in.size(); j++)
          d2 += CalcDist(points[res.in[j].first], points[res.in[j].second]);
        for (int j = 0; j < res.out.size(); j++)
          d2 -= CalcDist(points[res.out[j].first], points[res.out[j].second]);
        suc++;
        if (!res.kick_in.empty()) {
          p1.ExecuteMove(res.kick_in, res.kick_out);;
          for (int j = 0; j < res.kick_in.size(); j++)
            d1 += CalcDist(points[res.kick_in[j].first], points[res.kick_in[j].second]);
          for (int j = 0; j < res.kick_out.size(); j++)
            d1 -= CalcDist(points[res.kick_out[j].first], points[res.kick_out[j].second]);
        
        }
      }
    }
  }
  if (suc + 10< unstucks) {
    unstuck_coef += 0.1;
  }
  d1 = p1.GetLength(points), d2 = p2.GetLength(points);
  printf("unstuck end %lf %lf %lf %d in %d\n", d1, d2, (d1+d2)/2, unstucks, k);
}

boost::mutex mutex;
double best_max;
double best_av;
extern vector<double> good_non_seq;
extern int visits;

void Optimize(int Klim) {
  while (true) {
    mutex.lock();
    FILE *fpath = fopen("sub.avopt", "r");
    vector<int> vp1, vp2;
    for (int i = 0; i < points.size(); i++) {
      int a, b; fscanf(fpath, "%d,%d", &a, &b);
      vp1.push_back(a); vp2.push_back(b);
    }
    fclose(fpath);
    mutex.unlock();

    Path p1(vp1);
    Path p2(vp2);
    Unstuck(p1, p2, points);

    double d1 = p1.GetLength(points), d2 = p2.GetLength(points);

    printf("loaded\n");

    visits = 0;
    while (true) {
      double d1 = p1.GetLength(points), d2 = p2.GetLength(points);
      mutex.lock();
      if (max(d1, d2) < best_max) {
        best_max = max(d1, d2);
        PrintPaths(p1, p2, "sub.maxopt");
      }
      if ((d1 + d2)/2 < best_av) {
        best_av = (d1 + d2)/2;
        PrintPaths(p1, p2, "sub.avopt");
      }
      printf("start %lf %lf %lf %lf %lf\n", d1, d2, (d1+d2)/2, best_max, best_av);
      FILE *flog = fopen("opt.log", "a");
      fprintf(flog, "%d %lf %lf\n", time(NULL)-1355852120, max(d1, d2), (d1+d2)/2);
      fclose(flog);
      mutex.unlock();
      vector<int> oldp1 = p1.p_;
      vector<int> oldp2 = p2.p_;
      for (int i = 10; i < 150000 - 10; i++) {
//        printf("i %d\n", i);
        if (i % 5000 == 0) {
          double d1 = p1.GetLength(points), d2 = p2.GetLength(points);
          mutex.lock();
          if (max(d1, d2) < best_max) {
            best_max = max(d1, d2);
            PrintPaths(p1, p2, "sub.maxopt");
          }
          if ((d1 + d2)/2 < best_av) {
            best_av = (d1 + d2)/2;
            PrintPaths(p1, p2, "sub.avopt");
          }
          printf("%d %lf %lf %lf %lf %lf %d\n", i, d1, d2, (d1+d2)/2, best_max, best_av, visits);
          mutex.unlock();
//          PrintSeriesStats(good_non_seq, true);
        }
        if (d1 > d2) {    
          KOptMove res = p1.GetKOptMove(p2, p1.p_[i], p1.p_[i+1], closest, points, Klim, 1);
          if (res.valid) {
            p1.ExecuteMove(res);
//            printf("done 1\n");
            for (int j = 0; j < res.in.size(); j++)
              d1 += CalcDist(points[res.in[j].first], points[res.in[j].second]);
            for (int j = 0; j < res.out.size(); j++)
              d1 -= CalcDist(points[res.out[j].first], points[res.out[j].second]);
            if (!res.kick_in.empty()) {
//              printf("kick\n");
              p2.ExecuteMove(res.kick_in, res.kick_out);
              for (int j = 0; j < res.kick_in.size(); j++)
                d2 += CalcDist(points[res.kick_in[j].first], points[res.kick_in[j].second]);
              for (int j = 0; j < res.kick_out.size(); j++)
                d2 -= CalcDist(points[res.kick_out[j].first], points[res.kick_out[j].second]);
            }
/*            printf("done2\n");
            if (res.non_seq)
              break;*/
          } 
        } else {
          KOptMove res = p2.GetKOptMove(p1, p2.p_[i], p2.p_[i+1], closest, points, Klim, 1);
          if (res.valid) {
//            printf("base move\n");
            p2.ExecuteMove(res);
            for (int j = 0; j < res.in.size(); j++)
              d2 += CalcDist(points[res.in[j].first], points[res.in[j].second]);
            for (int j = 0; j < res.out.size(); j++)
              d2 -= CalcDist(points[res.out[j].first], points[res.out[j].second]);
            if (!res.kick_in.empty()) {
//              printf("kick move\n");
              p1.ExecuteMove(res.kick_in, res.kick_out);;
              for (int j = 0; j < res.kick_in.size(); j++)
                d1 += CalcDist(points[res.kick_in[j].first], points[res.kick_in[j].second]);
              for (int j = 0; j < res.kick_out.size(); j++)
                d1 -= CalcDist(points[res.kick_out[j].first], points[res.kick_out[j].second]);
            }
/*            if (res.non_seq)
              break;*/
          } 
        }
      }
      d1 = p1.GetLength(points), d2 = p2.GetLength(points);
      mutex.lock();
      if (max(d1, d2) < best_max) {
        best_max = max(d1, d2);
        PrintPaths(p1, p2, "sub.maxopt");
      }
      if ((d1 + d2)/2 < best_av) {
        best_av = (d1 + d2)/2;
        PrintPaths(p1, p2, "sub.avopt");
      }
      mutex.unlock();
      /*if (oldp1 == p1.p_ && oldp2 == p2.p_) {
        break;
      }*/
      break;
    }
  }
}

int main(int argc, char**argv) {
  srand(time(NULL));
  //TestValidMove();
  //return 0;

  FILE *fpoints = fopen(argv[2], "r");
  FILE *fclosest = fopen(argv[3], "r");

  int id, x, y;
  while (fscanf(fpoints, "%d %d %d", &id, &x, &y)>0) {
    points[id] = make_pair(x, y);
  }

  for (int i = 0; i < points.size(); i++) {
    vector<int> x;
    for (int j = 0; j < 100; j++) {
      int a; fscanf(fclosest, "%d", &a);
      if (x.size() < 7)
        x.push_back(a);
    }
    closest.push_back(x);
  }
  FILE *fpath = fopen("sub.avopt", "r");
  vector<int> vp1, vp2;
  for (int i = 0; i < points.size(); i++) {
    int a, b; fscanf(fpath, "%d,%d", &a, &b);
    vp1.push_back(a); vp2.push_back(b);
  }
  fclose(fpath);

  Path p1(vp1);
  Path p2(vp2);

  double d1 = p1.GetLength(points), d2 = p2.GetLength(points);
  best_max = max(d1, d2);
  best_av = (d1 + d2)/2;

  boost::thread th1(Optimize, 5);
/*  boost::thread th2(Optimize, 5);
  boost::thread th3(Optimize, 5);*/
/*  boost::thread th3(Optimize);
  boost::thread th4(Optimize);
  boost::thread th5(Optimize);
  boost::thread th6(Optimize);*/

  th1.join();
}
