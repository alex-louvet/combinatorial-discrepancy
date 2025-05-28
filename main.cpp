#include <algorithm>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <pthread.h>
#include <stdlib.h>
#include <string>
#include <unistd.h>
#include <vector>

#include "discrepancy.cpp"

using namespace std;

// Write CSV file wit the data of a set system
void writeCSVFile(SetSystem res, string fileName) {
  ofstream MyFile(fileName);
  for (Point &p : res.points) {
    for (const float &f : p.coordinates) {
      MyFile << f << ",";
    }
    MyFile << "\n";
  }
  MyFile << "sets\n";
  for (Set &s : res.sets) {
    for (const bool &b : s.points) {
      MyFile << b << ",";
    }
    MyFile << "\n";
  }
}

vector<int> simple_tokenizer(string s) {
  vector<int> res;
  stringstream ss(s);
  string word;
  while (ss >> word) {
    res.push_back(stoi(word));
  }
  return res;
}

bool isPrime(int n) {
  for (int i = 2; i <= sqrt(n); i++) {
    if (n % i == 0) {
      return false;
    }
  }
  return true;
}

int main(int argc, char **argv) {
  int r = time(NULL);
  srand(r);
  int n = 100;
  size_t d = 2;
  int m = 100;
  int c;
  while ((c = getopt(argc, argv, "a:n:t:d:f:r:p:m:i:c:k:es")) != -1) {
    switch (c) {
    case 'n':
      n = stoi(optarg);
      break;
    case 'd':
      d = stoi(optarg);
      break;
    case 'm':
      m = stoi(optarg);
      break;
    case '?':
      if (optopt == 'a' || optopt == 'n' || optopt == 't' || optopt == 'd' ||
          optopt == 'f' || optopt == 'r' || optopt == 'p' || optopt == 'm' ||
          optopt == 'i' || optopt == 'c' || optopt == 'k' || optopt == 'g')
        fprintf(stderr, "Option -%c requires an argument.\n", optopt);
      else if (isprint(optopt))
        fprintf(stderr, "Unknown option `-%c'.\n", optopt);
      else
        fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
      return 1;
    default:
      abort();
    }
  }
  SetSystem ss = RandomHyperplanes(n, d, m);
  vector<float> cst;
  for (int j = 0; j < m; j++)
    cst.push_back(sqrt(sumSet(ss.sets.at(j))));
  Coloring res = lm(ss, cst);
  for (int i = 0; i < n; i++) {
    if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) >
        abs(res.colors.at(i))) {
      res.colors.at(i) = (res.colors.at(i) < 0) ? 1 : -1;
    } else {
      res.colors.at(i) = (res.colors.at(i) < 0) ? -1 : 1;
    }
  }
  vector<float> discrepancies;
  for (int j = 0; j < m; j++)
    discrepancies.push_back(abs(dots(res.colors, ss.sets.at(j).points)));
  cout << "discrepancy: "
       << *max_element(discrepancies.begin(), discrepancies.end()) << endl;
  return 0;
}
