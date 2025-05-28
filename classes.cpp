#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <deque>
#include <fstream>
#include <iostream>
#include <math.h>
#include <random>
#include <tuple>
#include <vector>

#include "Eigen/Dense"

using namespace std;
#define PI 3.14159265

class Point {
public:
  vector<float> coordinates;
  vector<int> belongs_to;
  vector<int> not_in;

  Point(int dimension) {
    for (int i = 0; i < dimension; i++) {
      coordinates.push_back(static_cast<float>(rand()) /
                            static_cast<float>(RAND_MAX));
    }
  }

  Point(vector<float> coords) { coordinates = coords; }
};

struct {
  bool operator()(tuple<int, float> a, tuple<int, float> b) const {
    return get<1>(a) < get<1>(b);
  }
} floatWeightOrder;

class Edge {
public:
  int points[2];
  int weight;
  int next_ite_weight;

  Edge(int index1, int index2) {
    points[0] = index1;
    points[1] = index2;
    weight = 0;
    next_ite_weight = 0;
  }

  void increase() { weight++; }
  void decrease() { weight--; }
  void resetWeight() { weight = 0; }
};

class Set {
public:
  vector<bool> points;
  int weight;
  vector<int> points_indices;
  vector<int> complement_indices;

  Set(vector<bool> indices) {
    points = indices;
    weight = 0;
  }
  Set(int size) {
    for (int i = 0; i < size; i++) {
      points.push_back(0);
    }
    weight = 0;
  }

  void increase() { weight++; }
  void decrease() { weight--; }
  void resetWeight() { weight = 0; }
  void buildAdjacency() {
    points_indices.clear();
    complement_indices.clear();
    for (size_t i = 0; i < points.size(); i++) {
      if (points.at(i)) {
        points_indices.push_back(i);
      } else {
        complement_indices.push_back(i);
      }
    }
  }
};

class SetSystem {
public:
  vector<Point> points;
  vector<Set> sets;

  SetSystem(int n, int d, int m) {
    for (int i = 0; i < n; i++) {
      points.push_back(Point(d));
    }
    for (int i = 0; i < n; i++) {
      sets.push_back(Set(n));
    }
  }

  SetSystem(vector<Point> p, vector<Set> s) {
    points = p;
    sets = s;
  }

  SetSystem() {
    vector<Point> p;
    vector<Set> s;
    points = p;
    sets = s;
  }

  SetSystem(string filename) {
    vector<Point> p;
    vector<Set> s;

    ifstream ss(filename);
    string line;

    bool test = false;

    if (ss.is_open()) {
      while (getline(ss, line)) {
        if (!(line == "sets")) {
          if (!test) {
            float a = stof(line.substr(0, line.find(",")));
            line.erase(0, line.find(",") + 1);
            line.pop_back();
            float b = stof(line);
            p.push_back(Point({a, b}));
          } else {
            vector<bool> temp;
            for (size_t i = 0; i < line.size() / 2; i++) {
              temp.push_back(line.at(2 * i) - '0');
            }
            s.push_back(Set(temp));
          }
        } else {
          test = true;
        }
      }
    }
    points = p;
    sets = s;
  }

  void buildAdjacency(bool pts) {
    for (Point &p : points) {
      p.belongs_to.clear();
    }
#pragma omp parallel for
    for (size_t j = 0; j < sets.size(); j++) {
      sets.at(j).buildAdjacency();
      if (pts) {
        for (size_t i = 0; i < sets.at(j).points.size(); i++) {
          if (sets.at(j).points.at(i)) {
            points.at(i).belongs_to.push_back(j);
          } else {
            points.at(i).not_in.push_back(j);
          }
        }
      }
    }
  }

  void resetWeight() {
    for (Set &s : sets) {
      s.resetWeight();
    }
  }

  SetSystem project(vector<bool> projection) {
    assert(projection.size() == points.size());
    vector<Point> p;
    vector<Set> s;
    for (int i = 0; i < points.size(); i++) {
      if (projection.at(i)) {
        p.push_back(points.at(i));
      }
    }
    for (int j = 0; j < sets.size(); j++) {
      vector<bool> temp;
      for (int i = 0; i < points.size(); i++) {
        if (projection.at(i)) {
          temp.push_back(sets.at(j).points.at(i));
        }
      }
      s.push_back(temp);
    }
    return SetSystem(p, s);
  }
};

Set RandomSet(int size) {
  vector<bool> points;
  for (int i = 0; i < size; i++) {
    points.push_back(rand() > 0.5 * RAND_MAX);
  }
  return Set(points);
}

SetSystem Grid(int n, int d) {
  vector<Point> p;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    p.push_back(Point(d));
  }

  vector<bool> temp;
  for (int dim = 0; dim < d; dim++) {
    for (int i = 0; i < pow(n, 1.0 / d); i++) {
      temp.clear();
      for (const Point &pt : p) {
        if (pt.coordinates[dim] > i / (pow(n, 1.0 / d))) {
          temp.push_back(1);
        } else {
          temp.push_back(0);
        }
      }
      s.push_back(Set(temp));
    }
  }

  random_device rd;
  mt19937 g(rd());
  for (size_t i = 0; i < 10 * s.size(); i++) {
    shuffle(s.begin(), s.end(), g);
  }

  return SetSystem(p, s);
}

SetSystem GridWithCenters(int n, int d, int c) {
  vector<Point> p;
  vector<Set> s;
  vector<Point> centers;

  for (int i = 0; i < c; i++) {
    centers.push_back(Point(d));
  }

  random_device rd2{};
  mt19937 geng{rd2()};

  cout << "Centers are:" << endl;
  for (int j = 0; j < c; j++) {
    p.push_back(centers.at(j));
    for (int k = 0; k < d; k++) {
      cout << centers.at(j).coordinates.at(k) << " ";
    }
    cout << endl;
    for (int i = 0; i < n / c; i++) {
      vector<float> temp;
      for (int k = 0; k < d; k++) {
        normal_distribution<float> dist(centers.at(j).coordinates.at(k), 0.1);
        temp.push_back(dist(geng));
      }
      p.push_back(Point(temp));
    }
  }

  vector<bool> temp;
  for (int dim = 0; dim < d; dim++) {
    for (int i = 0; i < pow(n, 1.0 / d); i++) {
      temp.clear();
      for (const Point &pt : p) {
        if (pt.coordinates[dim] > i / (pow(n, 1.0 / d))) {
          temp.push_back(1);
        } else {
          temp.push_back(0);
        }
      }
      s.push_back(Set(temp));
    }
  }

  random_device rd;
  mt19937 g(rd());
  for (size_t i = 0; i < 10 * s.size(); i++) {
    shuffle(s.begin(), s.end(), g);
  }

  return SetSystem(p, s);
}

float dot(vector<float> p, vector<float> q) // dot product of two vectors
{
  assert(p.size() == q.size());
  float sum = 0.0f;
  for (size_t i = 0; i < p.size(); i++)
    sum += p[i] * q[i];
  return sum;
}

SetSystem RandomHyperplanes(int n, size_t d, int m) {
  vector<Point> p;
  vector<Set> s;
  vector<int> sample;
  for (int i = 0; i < n; i++) {
    p.push_back(Point(d));
  }
  for (int z = 0; z < m; z++) {
    sample.clear();
    while (sample.size() < d) {
      bool push = true;
      int temp = rand() % n;
      for (const int &c : sample) {
        if (c == temp) {
          push = false;
        }
      }
      if (push) {
        sample.push_back(rand() % n);
      }
    }
    Eigen::MatrixXd A(d, d);
    for (size_t i = 0; i < d; i++)
      for (size_t j = 0; j < d; j++)
        A(i, j) = p.at(sample.at(i)).coordinates.at(j);
    Eigen::VectorXd b(d);
    for (size_t i = 0; i < d; i++)
      b(i) = 1.0f;
    // VectorXd ans = A.colPivHouseholderQr().solve(b);
    Eigen::VectorXd ans = A.partialPivLu().solve(b);
    vector<float> N;
    N.reserve(d);
    for (size_t i = 0; i < d; i++)
      N.push_back(ans(i));
    Set set = Set(n);
    for (int i = 0; i < n; i++) {
      set.points.at(i) = (dot(p.at(i).coordinates, N) > 1.0) ? 1 : 0;
    }
    for (const int &index : sample) {
      set.points.at(index) = 1;
    }
    s.push_back(set);
  }

  return SetSystem(p, s);
}

SetSystem GridGraph(int n, int d) {
  vector<Point> pts;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      pts.push_back(Point(
          vector<float>{static_cast<float>(i) / n, static_cast<float>(j) / n}));
    }
  }

  const int pts_len = pts.size();

  vector<Set> sets;
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      Set temp = Set(vector<bool>(pts_len, 0));
      temp.points.at(i * n + j) = 1;
      for (int d1 = 0; d1 <= d; d1++) {
        for (int d2 = 0; d2 <= d1; d2++) {
          int d3 = d1 - d2;
          if (j >= d2) {
            if (i >= d3) {
              temp.points.at(i * n + j - d2 - n * d3) = 1;
            }
            if (i < n - d3) {
              temp.points.at(i * n + j - d2 + n * d3) = 1;
            }
          }
          if (j < n - d2) {
            if (i >= d3) {
              temp.points.at(i * n + j + d2 - n * d3) = 1;
            }
            if (i < n - d3) {
              temp.points.at(i * n + j + d2 + n * d3) = 1;
            }
          }
        }
      }
      sets.push_back(temp);
    }
  }

  return SetSystem(pts, sets);
}

SetSystem LinearGrid(int n, int d) {
  vector<Point> p;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    p.push_back(Point(d));
  }

  vector<bool> temp;
  for (int dim = 0; dim < d; dim++) {
    for (int i = 0; i < pow(n, 1.0 / d); i++) {
      temp.clear();
      for (const Point &pt : p) {
        if (pt.coordinates[dim] > i / (pow(n, 1.0 / d))) {
          temp.push_back(1);
        } else {
          temp.push_back(0);
        }
      }
      for (int j = 0; j < i; j++) {
        s.push_back(Set(temp));
      }
    }
  }

  random_device rd;
  mt19937 g(rd());
  for (size_t i = 0; i < 10 * s.size(); i++) {
    shuffle(s.begin(), s.end(), g);
  }

  return SetSystem(p, s);
}

SetSystem ExponentialGrid(int n, int d) {
  vector<Point> p;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    p.push_back(Point(d));
  }

  vector<bool> temp;
  for (int dim = 0; dim < d; dim++) {
    for (int i = 0; i < pow(n, 1.0 / d); i++) {
      temp.clear();
      for (const Point &pt : p) {
        if (pt.coordinates[dim] > i / (pow(n, 1.0 / d))) {
          temp.push_back(1);
        } else {
          temp.push_back(0);
        }
      }
      for (int j = 0; j < pow(2, i); j++) {
        s.push_back(Set(temp));
      }
    }
  }

  random_device rd;
  mt19937 g(rd());
  for (size_t i = 0; i < 10 * s.size(); i++) {
    shuffle(s.begin(), s.end(), g);
  }

  return SetSystem(p, s);
}

SetSystem DirectionalGrid(int n, int d) {
  vector<Point> p;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    p.push_back(Point(d));
  }

  vector<bool> temp;
  for (int dim = 0; dim < d; dim++) {
    for (int i = 0; i < pow(n, 1.0 / d); i++) {
      temp.clear();
      for (const Point &pt : p) {
        if (pt.coordinates[dim] > i / (pow(n, 1.0 / d))) {
          temp.push_back(1);
        } else {
          temp.push_back(0);
        }
      }
      if (dim == 0) {
        for (int j = 0; j < pow(n, 1.0 / d); j++) {
          s.push_back(Set(temp));
        }
      } else {
        s.push_back(Set(temp));
      }
    }
  }

  random_device rd;
  mt19937 g(rd());
  for (size_t i = 0; i < 10 * s.size(); i++) {
    shuffle(s.begin(), s.end(), g);
  }

  return SetSystem(p, s);
}

SetSystem Random(int n, int d, int m, float p) {
  vector<Point> pt;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    pt.push_back(Point(2));
  }
  for (int j = 0; j < m; j++) {
    vector<bool> temp;
    temp.clear();
    for (int i = 0; i < n; i++) {
      if (static_cast<float>(rand()) / RAND_MAX < p) {
        temp.push_back(1);
      } else {
        temp.push_back(0);
      }
    }
    s.push_back(temp);
  }
  return SetSystem(pt, s);
}

/*
Credits to Salmelu for its code that I just adapted to my data model:
https://github.com/Salmelu/ProjectivePlane
*/
SetSystem ProjectivePlane(int o) {
  vector<Point> p;
  vector<Set> s;

  p.push_back(Point({1., 0., 0.}));

  for (float i = 0; i < o; i++) {
    p.push_back(Point({i, 1., 0.}));
  }

  for (float i = 0; i < o; i++) {
    for (float j = 0; j < o; j++) {
      p.push_back(Point({i, j, 1.}));
    }
  }

  vector<bool> unique(o * o * o, true);
  unique.at(0) = false;

  for (int a = 0; a < o; a++) {
    for (int b = 0; b < o; b++) {
      for (int c = 0; c < o; c++) {
        if (!unique[(a * o * o) + (b * o) + c])
          continue;
        vector<bool> temp;
        temp.clear();
        for (size_t i = 0; i < p.size(); i++) {
          if ((static_cast<int>(a * p.at(i).coordinates.at(0) +
                                b * p.at(i).coordinates.at(1) +
                                c * p.at(i).coordinates.at(2)) %
               o) == 0) {
            temp.push_back(1);
          } else {
            temp.push_back(0);
          }
        }
        for (int i = 2; i < o; i++) {
          unique.at(((i * a) % o) * o * o + ((i * b) % o) * o + (i * c) % o) =
              false;
        }
        s.push_back(temp);
      }
    }
  }

  return SetSystem(p, s);
}

SetSystem ERGGraph(int n, int d, float p) {
  vector<Point> pt;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    pt.push_back(Point(2));
  }

  vector<vector<int>> edges;

  for (int i = 0; i < n - 1; i++) {
    vector<int> temp;
    temp.clear();
    for (int j = i + 1; j < n; j++) {
      if (rand() / static_cast<float>(RAND_MAX) < p) {
        temp.push_back(j);
      }
    }
    edges.push_back(temp);
  }
  edges.push_back({});
  for (int i = 1; i < n; i++) {
    for (int &j : edges.at(i)) {
      if (j > i) {
        edges.at(j).push_back(i);
      }
    }
  }

#pragma omp for
  for (int i = 0; i < n; i++) {
    vector<int> steps(n, -1);
    vector<bool> visited(n, false);
    deque<int> file = {i};
    steps.at(i) = 0;
    visited.at(i) = true;
    while (file.size() > 0) {
      int a = file.at(0);
      file.pop_front();
      for (int &x : edges.at(a)) {
        if (!visited.at(x)) {
          visited.at(x) = true;
          if (steps.at(x) == -1) {
            steps.at(x) = steps.at(a) + 1;
          } else {
            steps.at(x) = min(steps.at(x), steps.at(a) + 1);
          }
          if (steps.at(x) < d) {
            file.push_back(x);
          }
        }
      }
    }
    vector<bool> temp;
    for (int i = 0; i < n; i++) {
      if (steps.at(i) > -1 && steps.at(i) <= d) {
        temp.push_back(1);
      } else {
        temp.push_back(0);
      }
    }
    s.push_back(temp);
  }

  return SetSystem(pt, s);
}

/*
Credits to Laurent Vienno et al. for its code that I just adapted to my data
model: https://gitlab.inria.fr/viennot/graph-vcdim
*/
SetSystem PowerLaw(int n, int d, float beta, int seed) {
  vector<Point> pt;
  vector<Set> s;

  for (int i = 0; i < n; i++) {
    pt.push_back(Point(2));
  }

  vector<vector<int>> edges;

  vector<double> p_deg(n);
  double nf = n;
  for (int k = 1; k < n; ++k)
    p_deg[k] = nf / pow((float)k, beta);
  p_deg[0] = 0.;
  discrete_distribution<> power_law(
      p_deg.begin(),
      p_deg.end()); // normalization is done by std::discrete_distribution
  // generate degrees :
  mt19937 rand_gen(seed);
  vector<int> deg(n);
  int m = 0;
  for (int u = 0; u < n; ++u) {
    deg[u] = power_law(rand_gen);
    m += deg[u];
  }
  // generate edges :
  vector<int> dst, src;
  for (int u = 0; u < n; ++u) {
    for (int e = 0; e < deg[u]; ++e) {
      dst.push_back(u);
      src.push_back(u);
    }
  }
  uniform_int_distribution<> rnd_int(0, m - 1);
  for (int e = m - 1; e > 0; --e) {
    int f = rnd_int(rand_gen) % (e + 1);
    swap(dst[e], dst[f]);
  }
  for (int u = 0, e = 0; u < n; ++u) {
    for (int i = 0; i < deg[u]; ++i, ++e) {
      while (u == dst[e]) {
        int f = rnd_int(rand_gen) % m;
        while (u == dst[f] || src[f] == dst[e]) {
          f = rnd_int(rand_gen) % m;
        }
        swap(dst[e], dst[f]);
      }
    }
  }
  int current = 0;
  vector<int> temp;
  for (size_t e = 0; e < src.size(); e++) {
    if (src[e] != current) {
      edges.push_back(temp);
      current = src[e];
      temp.clear();
    }
    temp.push_back(dst[e]);
  }
  edges.push_back({});

  for (size_t e = 0; e < dst.size(); e++) {
    edges.at(dst[e]).push_back(src[e]);
  }

  // #pragma omp for
  for (int i = 0; i < n; i++) {
    vector<int> steps(n, -1);
    vector<bool> visited(n, false);
    deque<int> file = {i};
    steps.at(i) = 0;
    visited.at(i) = true;
    while (file.size() > 0) {
      int a = file.at(0);
      file.pop_front();
      for (int &x : edges.at(a)) {
        if (!visited.at(x)) {
          visited.at(x) = true;
          if (steps.at(x) == -1) {
            steps.at(x) = steps.at(a) + 1;
          } else {
            steps.at(x) = min(steps.at(x), steps.at(a) + 1);
          }
          if (steps.at(x) < d) {
            file.push_back(x);
          }
        }
      }
    }
    vector<bool> temp;
    for (int i = 0; i < n; i++) {
      if (steps.at(i) > -1 && steps.at(i) <= d) {
        temp.push_back(1);
      } else {
        temp.push_back(0);
      }
    }
    s.push_back(temp);
  }

  return SetSystem(pt, s);
}

float euclidian_distance(Point p1, Point p2) {
  return sqrt((p1.coordinates.at(0) - p2.coordinates.at(0)) *
                  (p1.coordinates.at(0) - p2.coordinates.at(0)) +
              (p1.coordinates.at(1) - p2.coordinates.at(1)) *
                  (p1.coordinates.at(1) - p2.coordinates.at(1)));
}

SetSystem ConcentricCircles(int n, int m) {
  vector<Point> pt;
  vector<Set> s;

  for (float j = 1.0; j < 5.0; j++) {
    for (int i = 0; i < n / 4; i++) {
      float theta = rand() / static_cast<float>(RAND_MAX) * 2 * PI;
      vector<float> temp;
      temp.push_back(0.5 + 0.1 * j * cos(theta));
      temp.push_back(0.5 + 0.1 * j * sin(theta));
      pt.push_back(Point(temp));
    }
  }

  for (int i = 0; i < m; i++) {
    vector<bool> temp;
    int c = rand() / static_cast<float>(RAND_MAX) * pt.size();
    int r = rand() / static_cast<float>(RAND_MAX) * pt.size();
    while (r == c) {
      r = rand() / static_cast<float>(RAND_MAX) * pt.size();
    }
    for (Point &p : pt) {
      if (euclidian_distance(p, pt.at(c)) <=
          euclidian_distance(pt.at(c), pt.at(r))) {
        temp.push_back(1);
      } else {
        temp.push_back(0);
      }
    }
    s.push_back(temp);
  }

  return SetSystem(pt, s);
}

class Result : public SetSystem {
public:
  vector<float> weights;
  vector<int> intersections;
  Result(int n, int d, int m) {
    for (int i = 0; i < n; i++) {
      points.push_back(Point(d));
    }
    for (int i = 0; i < n; i++) {
      sets.push_back(Set(n));
    }
  }

  Result(vector<Point> p, vector<Set> s) {
    points = p;
    sets = s;
  }

  Result() {
    vector<Point> p;
    vector<Set> s;
    points = p;
    sets = s;
  }

  Result(SetSystem ss) {
    points = ss.points;
    sets = ss.sets;
  }
};

class Coloring {
public:
  vector<float> colors;
  Coloring(Point p) {
    vector<float> c;
    for (int i = 0; i < p.coordinates.size(); i++) {
      c.push_back(p.coordinates.at(i));
    }
    colors = c;
  }
  Coloring(int n) {
    vector<float> c;
    for (int i = 0; i < n; i++) {
      c.push_back(0.);
    }
    colors = c;
  }
  Coloring() {
    vector<float> c;
    colors = c;
  }
};
