#include <cmath>
#include <cstdlib>
#include <math.h>
#include <ostream>
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
    Eigen::MatrixXd A(n + m, n);
    int A_counter = 0;
    for (int i = 0; i < n; i++) {
      if (res.colors.at(i) > 1 - 1 / (8 * log(m)) or
          res.colors.at(i) < 1 / (8 * log(m)) - 1) {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (k == i) ? 1.0f : 0.0f;
        }
      } else {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = 0.0f;
        }
      }
      A_counter++;
    }
    for (int j = 0; j < m; j++) {
      if (dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) >=
              constraints.at(j) - 1 / (8 * log(m)) or
          dots(vec_diff(res.colors, x0.coordinates), ss.sets.at(j).points) <=
              1 / (8 * log(m)) - constraints.at(j)) {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = (ss.sets.at(j).points.at(k)) ? 1.0f : 0.0f;
        }
      } else {
        for (int k = 0; k < n; k++) {
          A(A_counter, k) = 0.0f;
        }
      }
      A_counter++;
    }
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    Eigen::MatrixXd A_null_space = lu.kernel();
    vector<vector<float>> A_null_space_vec = convert_eigen(A_null_space);
    for (int k = 0; k < A_null_space_vec.size(); k++) {
      vector<float> add =
          vec_scal(A_null_space_vec.at(k), random_int() / (8 * log(m)));
      for (int i = 0; i < n; i++)
        res.colors.at(i) += add.at(i);
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
  while (color_map.size() > 0) {
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
  return res;
}
