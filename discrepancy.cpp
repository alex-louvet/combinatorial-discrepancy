#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <math.h>
#include <ostream>
#include <sched.h>
#include <vector>

#include "classes.cpp"

using namespace std;

bool intersects(Edge e, Set s) {
  return s.points.at(e.points[0]) != s.points.at(e.points[1]);
}

int sumSet(Set s) {
  int res = 0;
  for (const bool &value : s.points) {
    res += value;
  }
  return res;
}

float n1(vector<float> v) {
  float res = 0;
  for (const float &value : v) {
    res += value;
  }
  return res;
}

float n2(vector<float> v) {
  float res = 0;
  for (const float &value : v) {
    res += value * value;
  }
  return sqrt(res);
}

int posNum(vector<unsigned long> weight) {
  int res = 0;
  for (const long unsigned int &value : weight) {
    if (value >= 0) {
      res++;
    }
  }
  return res;
}

bool valid_partition_size_list(vector<int> partition_size, int n) {
  int count = 0;
  for (int &x : partition_size) {
    if (x < 1) {
      return false;
    } else {
      count += x;
    }
  }
  if (count > n) {
    return false;
  }
  return true;
}

int partialSum(vector<int> v, int i) {
  int res = 0;
  for (int j = 0; j < i; j++) {
    res += v.at(j);
  }
  return res;
}

vector<float> vec_diff(vector<float> v1, vector<float> v2) {
  vector<float> res;
  for (int i = 0; i < v1.size(); i++) {
    res.push_back(v1.at(i) - v2.at(i));
  }
  return res;
}

vector<float> vec_sum(vector<float> v1, vector<float> v2) {
  vector<float> res;
  for (int i = 0; i < v1.size(); i++) {
    res.push_back(v1.at(i) + v2.at(i));
  }
  return res;
}

vector<float> vec_scal(vector<float> v1, float k) {
  vector<float> res;
  for (int i = 0; i < v1.size(); i++) {
    res.push_back(k * v1.at(i));
  }
  return res;
}

vector<float> vec_scalb(vector<bool> v1, float k) {
  vector<float> res;
  for (int i = 0; i < v1.size(); i++) {
    if (v1.at(i))
      res.push_back(k);
    else
      res.push_back(0);
  }
  return res;
}

float dots(vector<float> p, vector<bool> q) // dot product of two vectors
{
  assert(p.size() == q.size());
  float sum = 0.0f;
  for (size_t i = 0; i < p.size(); i++)
    if (q[i]) {
      sum += p[i];
    }
  return sum;
}

vector<vector<float>> convert_eigen(Eigen::MatrixXd M) {
  vector<vector<float>> res;
  for (int i = 0; i < M.cols(); i++) {
    vector<float> temp;
    for (int j = 0; j < M.rows(); j++) {
      temp.push_back(M(j, i));
    }
    res.push_back(temp);
  }
  return res;
}

Coloring lm_partial(SetSystem ss, vector<float> constraints, Point x0) {
  int n = ss.points.size();
  int m = ss.sets.size();
  std::random_device rd{};
  std::mt19937 gen{rd()};
  std::normal_distribution d{0.0, 1.0};
  auto random_int = [&d, &gen] { return d(gen); };
  Coloring res(x0);
  for (int t = 0; t < 64 * log(m) * log(m); t++) {
    int cst_count = 0;
    for (int i = 0; i < n; i++) {
      if (res.colors.at(i) > 1 - 1 / (8 * log(m)) or
          res.colors.at(i) < 1 / (8 * log(m)) - 1) {
        cst_count++;
      }
    }
    for (int j = 0; j < m; j++) {
      if (dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) >=
              constraints.at(j) - 1 / (8 * log(m)) or
          dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) <=
              1 / (8 * log(m)) - constraints.at(j)) {
        cst_count++;
      }
    }
    Eigen::MatrixXd A(max(1, cst_count), n);
    int A_counter = 0;
    for (int i = 0; i < n; i++) {
      if (res.colors.at(i) > 1 - 1 / (8 * log(m)) or
          res.colors.at(i) < 1 / (8 * log(m)) - 1) {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (k == i) ? 1.0f : 0.0f;
        }
        A_counter++;
      }
    }
    for (int j = 0; j < m; j++) {
      if (dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) >=
              constraints.at(j) - 1 / (8 * log(m)) or
          dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) <=
              1 / (8 * log(m)) - constraints.at(j)) {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (ss.sets.at(j).points.at(k)) ? 1.0f : 0.0f;
        }
        A_counter++;
      }
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    Eigen::MatrixXd A_null_space = lu.kernel();
    vector<vector<float>> A_null_space_vec = convert_eigen(A_null_space);
    for (int k = 0; k < A_null_space_vec.size(); k++) {
      vector<float> add = vec_scal(
          A_null_space_vec.at(k),
          random_int() / (8 * max(1.0f, n2(A_null_space_vec.at(k))) * log(m)));
      for (int i = 0; i < n; i++) {
        res.colors.at(i) += add.at(i);
      }
    }
  }
  return res;
}

Coloring lm(SetSystem ss, vector<float> constraints) {
  int n = ss.points.size();
  int m = ss.sets.size();
  Coloring res(n);
  vector<bool> colored(n, false);
  vector<int> color_map;
  for (int i = 0; i < n; i++) {
    color_map.push_back(i);
  }
  Point x0(vector<float>(n, 0.0));
  while (color_map.size() > 1) {
    vector<float> cst;
    for (int j = 0; j < m; j++)
      cst.push_back(2 * sqrt(sumSet(ss.sets.at(j))));
    Coloring temp = lm_partial(ss, cst, x0);
    vector<int> next_color_map;
    vector<float> next_x0;
    vector<bool> projection_vec;
    for (int i = 0; i < temp.colors.size(); i++) {
      if (!colored.at(color_map.at(i)) and
          (temp.colors.at(i) > 1 - 1 / (8 * log(m)) or
           temp.colors.at(i) < 1 / (8 * log(m)) - 1)) {
        res.colors.at(color_map.at(i)) = temp.colors.at(i);
        colored.at(color_map.at(i)) = true;
        projection_vec.push_back(false);
      } else {
        next_color_map.push_back(color_map.at(i));
        next_x0.push_back(temp.colors.at(i));
        projection_vec.push_back(true);
      }
    }
    color_map.assign(next_color_map.begin(), next_color_map.end());
    x0 = Point(next_x0);
    SetSystem tempp = ss.project(projection_vec);
    ss.points.assign(tempp.points.begin(), tempp.points.end());
    ss.sets.assign(tempp.sets.begin(), tempp.sets.end());
  }
  for (int i = 0; i < n; i++) {
    if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) >
        abs(res.colors.at(i))) {
      res.colors.at(i) = (res.colors.at(i) > 0) ? 1 : -1;
    } else {
      res.colors.at(i) = (res.colors.at(i) < 0) ? -1 : 1;
    }
  }
  return res;
}

bool weightcomp(tuple<float, int> a, tuple<float, int> b) {
  return get<0>(a) > get<0>(b);
}

Coloring lrr_partial(SetSystem ss, vector<float> constraints, Point x0) {
  assert(constraints.size() == ss.sets.size());
  int n = ss.points.size();
  int m = ss.sets.size();
  Coloring res(n);
  vector<double> weights(m, 0.);
  for (int j = 0; j < m; j++)
    weights.at(j) = exp(-1 * constraints.at(j) * constraints.at(j));
  int c = 0;
  while (true) {
    Eigen::MatrixXd A(max(n, 10), n);
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        A(i, j) = 0;
      }
    }
    int A_counter = 0;
    // U1
    for (int i = 0; i < n; i++) {
      if (res.colors.at(i) >= 0.99 or res.colors.at(i) <= -0.99) {
        res.colors.at(i) = (res.colors.at(i) > 0) ? 1 : -1;
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (k == i) ? 1.0f : 0.0f;
        }
        A_counter++;
      }
    }
    if (A_counter >= n / 2) {
      break;
    }
    // U2
    for (int k = 0; k < n; k++) {
      A(A_counter, k) = res.colors.at(k);
    }
    A_counter++;
    // U3
    vector<tuple<float, int>> tosort;
    for (int j = 0; j < m; j++) {
      tosort.push_back({weights.at(j), j});
    }
    sort(tosort.begin(), tosort.end(), weightcomp);
    for (int i = 0; i < min(n / 16, m); i++) {
      for (int k = 0; k < n; k++) {
        A(A_counter, k) =
            (ss.sets.at(get<1>(tosort.at(i))).points.at(k)) ? 1.0f : 0.0f;
      }
      A_counter++;
    }
    // U4
    for (int j = 0; j < m; j++) {
      if (constraints.at(j) <= 1 and constraints.at(j) > 0) {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (ss.sets.at(j).points.at(k)) ? 1.0f : 0.0f;
        }
        A_counter++;
      }
    }
    // U5
    vector<float> toadd(n, 0.);
    for (int j = 0; j < m; j++) {
      vector<float> temp = vec_sum(
          toadd,
          vec_scalb(ss.sets.at(j).points,
                    constraints.at(j) * weights.at(j) *
                        exp(-1 * (4 * constraints.at(j) * constraints.at(j)) /
                            (n * n))));
      toadd.assign(temp.begin(), temp.end());
    }
    for (int k = 0; k < n; k++)
      A(A_counter, k) = toadd.at(k);
    A_counter++;
    // U6
    Eigen::VectorXd b(n);
    for (int i = 0; i < n; i++)
      b(i) = ss.sets.at(0).points.at(i);
    Eigen::MatrixXd mat =
        (weights.at(0) * constraints.at(0) * constraints.at(0)) *
        (b * b.transpose());
    for (int j = 1; j < m; j++) {
      for (int i = 0; i < n; i++)
        b(i) = ss.sets.at(j).points.at(i);
      mat += (weights.at(j) * constraints.at(j) * constraints.at(j)) *
             (b * b.transpose());
    }
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver(mat);
    if (eigensolver.info() != Eigen::Success) {
      abort();
    }
    for (int k = n - n / 16 - 1; k < n; k++) {
      for (int i = 0; i < n; i++) {
        A(A_counter, i) = eigensolver.eigenvectors()(i, k);
      }
      A_counter++;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    Eigen::MatrixXd A_null_space = lu.kernel();
    vector<vector<float>> A_null_space_vec = convert_eigen(A_null_space);
    vector<float> add = A_null_space_vec.at(0);
    if (A_null_space_vec.size() == 1) {
      bool test = true;
      for (int i = 0; i < n; i++) {
        if (add.at(i) != 0)
          test = false;
        break;
      }
      if (test) {
        break;
      }
    }
    for (int k = 1; k < A_null_space_vec.size(); k++) {
      for (int i = 0; i < n; i++) {
        add.at(i) = add.at(i) + A_null_space_vec.at(k).at(i);
      }
    }
    float factor = 1;
    for (int i = 0; i < n; i++) {
      if (abs(res.colors.at(i)) < 1 &&
          1 / abs(add.at(i)) - res.colors.at(i) / add.at(i) < factor) {
        factor = 1 / abs(add.at(i)) - res.colors.at(i) / add.at(i);
      }
    }
    for (int i = 0; i < n; i++) {
      res.colors.at(i) += add.at(i) * factor;
    }

    for (int j = 0; j < m; j++) {
      weights.at(j) =
          weights.at(j) *
          exp(weights.at(j) / sqrt(n) *
              dots(vec_scal(add, factor), ss.sets.at(j).points)) *
          exp(-1 * (4 * constraints.at(j) * constraints.at(j)) / (n * n));
    }
  }
  return res;
}

Coloring lrr(SetSystem ss, vector<float> constraints) {
  int n = ss.points.size();
  int m = ss.sets.size();
  Coloring res(n);
  vector<bool> colored(n, false);
  vector<int> color_map;
  for (int i = 0; i < n; i++) {
    color_map.push_back(i);
  }
  Point x0(vector<float>(n, 0.0));
  while (color_map.size() > 2) {
    vector<float> cst;
    for (int j = 0; j < m; j++)
      cst.push_back(2 * sqrt(sumSet(ss.sets.at(j))));
    Coloring temp = lrr_partial(ss, cst, x0);
    vector<int> next_color_map;
    vector<float> next_x0;
    vector<bool> projection_vec;
    for (int i = 0; i < temp.colors.size(); i++) {
      if (!colored.at(color_map.at(i)) and
          (temp.colors.at(i) > 1 - 1 / (8 * log(m)) or
           temp.colors.at(i) < 1 / (8 * log(m)) - 1)) {
        res.colors.at(color_map.at(i)) = temp.colors.at(i);
        colored.at(color_map.at(i)) = true;
        projection_vec.push_back(false);
      } else {
        next_color_map.push_back(color_map.at(i));
        next_x0.push_back(temp.colors.at(i));
        projection_vec.push_back(true);
      }
    }
    if (color_map.size() == next_color_map.size()) {
      cout << "LRR exited at " << color_map.size() << " (might exit below 16)"
           << endl;
      for (int i = 0; i < temp.colors.size(); i++) {
        if (!colored.at(color_map.at(i))) {
          res.colors.at(color_map.at(i)) = (temp.colors.at(i) < 0) ? -1 : 1;
        }
      }
      break;
    }
    color_map.assign(next_color_map.begin(), next_color_map.end());
    x0 = Point(next_x0);
    SetSystem tempp = ss.project(projection_vec);
    ss.points.assign(tempp.points.begin(), tempp.points.end());
    ss.sets.assign(tempp.sets.begin(), tempp.sets.end());
  }
  for (int i = 0; i < n; i++) {
    if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) >
        abs(res.colors.at(i))) {
      res.colors.at(i) = (res.colors.at(i) > 0) ? 1 : -1;
    } else {
      res.colors.at(i) = (res.colors.at(i) < 0) ? -1 : 1;
    }
  }
  return res;
}
