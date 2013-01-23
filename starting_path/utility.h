#ifndef UTILITY_H__
#define UTILITY_H__

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <cstdio>
#include <tr1/unordered_set>
#include <deque>
#include <cmath>

namespace std {
  namespace tr1 {

    template<typename a, typename b>
      struct hash< std::pair<a, b> > {
        private:
          const hash<a> ah;
          const hash<b> bh;
        public:
          hash() : ah(), bh() {}
          size_t operator()(const std::pair<a, b> &p) const {
            return ah(p.first) ^ bh(p.second);
          }
      };

  }}

using namespace std;
using namespace std::tr1;

set<string> stop_words;
typedef pair<double, double> pdd;

string UnescapeShit(const string& x) {
  string ret = "";
  char last = 'x';
  for (int i = 0; i < x.length(); i++) {
    if (last != '\\' && x[i] != '\\')
      ret += x[i];
    if (x[i] == '\\')
      ret += ' ';
    last = x[i];
  }
  return ret;
}

string StripTags(const string& x) {
  bool intag = false;
  string ret = "";
  for (int i = 0; i < x.length(); i++) {
    if (x[i] == '<') intag = true;
    else if (x[i] == '>') intag = false;
    else if (intag == false) ret += x[i];
  }
  return ret;
}

void PrepareStopWords() {
  if (!stop_words.empty()) return;
  FILE *f = fopen("stop_words", "r");
  char buf[1000];
  while (fscanf(f, "%s", buf)>0) {
    stop_words.insert(buf);
  }
  fclose(f);
}

template<class T>
void VectorToSet(const vector<T>&a, set<T>* b) {
  for (int i = 0; i < a.size(); i++)
    b->insert(a[i]);
}

template<class T, class T2>
void SetToVector(const set<T>&a, vector<T2>* b) {
  for (auto it = a.begin(); it != a.end(); ++it)
    b->push_back(make_pair(*it, 0));
}

template<class TIt, class TType>
void SortBySecond(TIt start, TIt end, vector<TType>* out, double limit = 0) {
  out->clear();
  vector<pair<double, TType> > temp;
  for (TIt it = start; it != end; ++it) {
    if (it->second < limit) continue;
    temp.push_back(make_pair(it->second, it->first));
  }
  sort(temp.rbegin(), temp.rend());
  for (int i = 0; i < temp.size(); i++)
    out->push_back(temp[i].second);
}

template<class T>
double CalcAP(const vector<T>& rec_items, const set<T>& good_items) {
  if (good_items.size() == 0) return 0;
  double good = 0;
  double ret = 0;
  for (int i = 0; i < rec_items.size(); i++) {
    if (good_items.count(rec_items[i])) {
      good+=1;
      ret += (good / (i+1));
    }
  }
  return ret / good_items.size();
}

template<class T>
int CountInSet(const vector<T>& rec_items, const set<T>& good_items) {
  int ret = 0;
  for (int i = 0; i < rec_items.size(); i++) {
    if (good_items.count(rec_items[i])) ret++;
  }
  return ret;
}

inline void SplitStringToVector(const string&s, char sep, vector<string>* out) {
  string cur = "";
  for (int i = 0; i < s.length(); ++i) {
    if (s[i] == sep && !cur.empty()) {
      out->push_back(cur);
      cur = "";
    } else {
      if (s[i] != '\n' && s[i] != '\r')
        cur += s[i];
    }
  }
  if (!cur.empty()) {
    out->push_back(cur);
  }
}

inline void SplitStringToVectorE(const string&s, char sep, vector<string>* out) {
  string cur = "";
  bool ss = true;
  for (int i = 0; i < s.length(); ++i) {
    if (s[i] == sep) {
      ss = true;
      out->push_back(cur);
      cur = "";
    } else {
      ss = false;
      if (s[i] != '\n' && s[i] != '\r')
        cur += s[i];
    }
  }
  out->push_back(cur);
}

inline double StringToDouble(const string& s) {
  return atof(s.c_str());
}

inline int StringToInt(const string& s) {
  return atoi(s.c_str());
}

inline string RemoveNonAlpha(string s) {
  string ret;
  for (int i = 0; i < s.length(); i++) {
    if (!isascii(s[i])) {
      ret += ' ';
      continue;
    }
    if (s[i] == ' ' || isalpha(s[i])) {
      ret += tolower(s[i]);
    } else {
      ret += ' ';
    }
  }
  return ret;
}

inline string RemoveNonAlphaNum(string s) {
  string ret;
  for (int i = 0; i < s.length(); i++) {
    if (!isascii(s[i])) continue;
    if (s[i] == ' ' || isalnum(s[i]) || s[i] == '+' || s[i] == '-') {
      ret += tolower(s[i]);
    }
  }
  return ret;
}

inline string Trim(string s) {
  string ret;
  int b = s.length();
  int e = 0;
  for (int i = 0; i < s.length(); i++) {
    if (s[i] != ' ' && b == s.length()) {
      b = i;
    }
    if (s[i] != ' ') e = i+1;
  }
  for (int i = b; i < e; i++) {
    ret += s[i];
  }
  return ret;
}

inline string FixSpaces(string s) {
  char last = 'x';
  string ret;
  for (int i = 0; i < s.length(); i++) {
    if (last != ' ' || s[i] != ' ')
      ret += s[i];
    last = s[i];
  }
  return ret;
}

inline string RemoveStopWords(string s) {
  vector<string> words;
  SplitStringToVector(s, ' ', &words);
  string ret;
  for (int i = 0; i < words.size(); i++) {
//    if (words[i].length() < 2) continue;
    if (stop_words.count(words[i])) continue;
    if (!ret.empty()) ret += ' ';
    ret += words[i];
  }
  return ret;
}

inline string FixJobName(string job) {
  PrepareStopWords();
  job = RemoveStopWords(FixSpaces(Trim(RemoveNonAlpha(job))));
  return job;
}

inline string FixJobReq(string job) {
  PrepareStopWords();
  job = RemoveStopWords(FixSpaces(Trim(RemoveNonAlphaNum(job))));
  return job;
}

template<class T>
double CalcSimilarity(const vector<T>&a, const vector<T>& b) {
  unordered_set<T> uni;
  int inter_size = 0;
  for (int i = 0; i < a.size(); i++)
    uni.insert(a[i]);
  for (int i = 0; i < b.size(); i++) {
    if (uni.count(b[i])) inter_size++;
    uni.insert(b[i]);
  }
  return 1.0*inter_size/uni.size();
}

template<class T>
double CalcSimilarityDice(const vector<T>&a, const vector<T>& b, double beta) {
  unordered_set<T> uni;
  int inter_size = 0;
  for (int i = 0; i < a.size(); i++)
    uni.insert(a[i]);
  for (int i = 0; i < b.size(); i++) {
    if (uni.count(b[i])) inter_size++;
    uni.insert(b[i]);
  }
  return (1 + beta)*inter_size/(a.size() + b.size()*beta);
}

template<class T>
double CalcSimilarityAToB(const vector<T>&a, const vector<T>& b) {
  unordered_set<T> bs;
  int inter_size = 0;
  for (int i = 0; i < b.size(); i++)
    bs.insert(b[i]);
  for (int i = 0; i < a.size(); i++) {
    if (bs.count(a[i])) inter_size++;
  }
  return 1.0*inter_size/(a.size());
}

template<class T>
double CalcSimilarityBToA(const vector<T>&a, const vector<T>& b) {
  unordered_set<T> bs;
  int inter_size = 0;
  for (int i = 0; i < a.size(); i++)
    bs.insert(a[i]);
  for (int i = 0; i < b.size(); i++) {
    if (bs.count(b[i])) inter_size++;
  }
  return 1.0*inter_size/(b.size());
}
double CalcJaccardsSimilarity(const set<string>&a, const set<string>& b) {
  set<string> uni;
  int inter_size = 0;
  for (set<string>::iterator it = a.begin(); it != a.end(); ++it)
    uni.insert(*it);
  for (set<string>::iterator it = b.begin(); it != b.end(); ++it) {
    if (uni.count(*it)) inter_size++;
    uni.insert(*it);
  }
  if (uni.size() == 0) return 0;
  return 1.0*inter_size/uni.size();
}

inline pair<double, double> CalcMeanAndVar(const vector<float>& x) {
  double av = 0;
  for (int i = 0; i < x.size(); i++) {
    av += x[i];
  }
  if (x.size() > 0)
    av /= x.size();
  double var = 0;
  for (int i = 0; i < x.size(); i++) 
    var += (x[i] - av)*(x[i] - av);
  if (x.size() > 0)
    var /= x.size();
  else
    var = 1;
  return make_pair(av, var);
}

inline pair<double, double> CalcMeanAndVar(const vector<double>& x) {
  double av = 0;
  for (int i = 0; i < x.size(); i++) {
    av += x[i];
  }
  if (x.size() > 0)
    av /= x.size();
  double var = 0;
  for (int i = 0; i < x.size(); i++) 
    var += (x[i] - av)*(x[i] - av);
  if (x.size() > 0)
    var /= x.size();
  else
    var = 1;
  return make_pair(av, var);
}

inline pair<double, double> CalcMeanAndVar(const deque<double>& x) {
  double av = 0;
  for (int i = 0; i < x.size(); i++) {
    av += x[i];
  }
  av /= x.size();
  double var = 0;
  for (int i = 0; i < x.size(); i++) 
    var += (x[i] - av)*(x[i] - av);
  var /= x.size();
  return make_pair(av, var);
}

inline void PrintSeriesStats(vector<double> series, bool detailed = false) {
  pdd mv = CalcMeanAndVar(series);
  double m3 = 0, m4 = 0;
  for (int i = 0; i < series.size(); i++) {
    double x = series[i];
    m3 += x*x*x;
    m4 += x*x*x*x;
  }
  m3 /= series.size(); m4 /= series.size();
  m4 /= mv.second;
  m4 /= mv.second;
  sort(series.begin(), series.end());
  int l = series.size() - 1;
  if (series.size() == 0)
    series.push_back(-1);
  int up = 0, down = 0;
  for (int i = 0; i < series.size(); i++) {
    if (series[i] > 0) up++;
    else down++;
  }
  if (detailed)
    fprintf(stderr, "%6d %8.4lf+-%6.3lf %5.2lf%% %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n",
        l+1,
        mv.first, sqrt(mv.second), 100.0*up/(up+down),
        series[0], series[l/10], series[l/4], series[l/2],
        series[l*3/4], series[l*9/10],
        series[l], m4);
  else
    fprintf(stderr, "%6d %8.5lf+-%6.4lf %5.2lf%% %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf\n", l+1,
        mv.first, sqrt(mv.second), 100.0*up/(up+down),
        series[0], series[l/10], series[l/2], series[l*9/10],
        series[l]);
}

inline void PrintSeriesStats2(vector<double> series, bool detailed = false) {
  pdd mv = CalcMeanAndVar(series);
  sort(series.begin(), series.end());
  int l = series.size() - 1;
  if (series.size() == 0)
    series.push_back(-1);
  int up = 0, down = 0;
  for (int i = 0; i < series.size(); i++) {
    if (series[i] > 0) up++;
    else down++;
  }
  if (detailed)
    fprintf(stderr, "%6d %8.4lf %5.2lf%% %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf %8.4lf\n",
        l+1,
        mv.first, 100.0*up/(up+down),
        series[0], series[l/10], series[l/4], series[l/2],
        series[l*3/4], series[l*9/10],
        series[l]);
  else
    fprintf(stderr, "%6d %8.5lf %5.2lf%% %8.5lf %8.5lf %8.5lf %8.5lf %8.5lf\n", l+1,
        mv.first, 100.0*up/(up+down),
        series[0], series[l/10], series[l/2], series[l*9/10],
        series[l]);
}

double torad(double x) {
  return x/180.0*M_PI;
}

double CalcEarthDist(pair<double, double> a, pair<double, double> b) {
  double lat1 = torad(a.first);
  double lat2 = torad(b.first);
  double dlat = torad(a.first-b.first);
  double dlong = torad(a.second-b.second);
  double R = 6371*1000;
  double q = sin(dlat/2)*sin(dlat/2) + cos(lat1)*cos(lat2)*sin(dlong/2)*sin(dlong/2);
  double c = atan2(sqrt(q), sqrt(1-q));
  return c*R;
}

double CalcLen(pair<double, double> a) {
  return sqrt(a.first*a.first + a.second*a.second);
}

double CalcDist(pair<double, double> a, pair<double, double> b) {
  return CalcLen(make_pair(a.first - b.first, a.second - b.second));
}

pair<double, double> NormalizePair(pair<double, double> x) {
  double l = CalcLen(x);
  return make_pair(x.first / l, x.second / l);
}

void StringToBigrams(const string& in, vector<string>* out) {
  vector<string> words;
  SplitStringToVector(in, ' ', &words);
  if (words.size() == 1) {
    out->push_back(words[0]);
  } else {
    for (int i = 1; i < words.size(); i++)
      out->push_back(words[i-1] + " " + words[i]);
  }
}

void RemoveFromVector(vector<int>& x, int pos) {
  swap(x[pos], x[x.size()-1]);
  x.pop_back();
}

#endif
