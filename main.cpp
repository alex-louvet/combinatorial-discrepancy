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

int main() {
  int n = 16;
  int m = 16;
  size_t d = 2;
  Point x0(vector<float>(0, n));
  SetSystem ss = RandomHyperplanes(n, d, m);
  vector<float> cst;
  for (int j = 0; j < m; j++)
    cst.push_back(sqrt(sumSet(ss.sets.at(j))));
  Coloring res = lrr(ss, cst);
  Coloring reslm = lm(ss, cst);
  vector<float> discrepancies;
  for (int j = 0; j < m; j++)
    discrepancies.push_back(abs(dots(res.colors, ss.sets.at(j).points)));
  cout << "discrepancy lrr: "
       << *max_element(discrepancies.begin(), discrepancies.end()) << endl;
  vector<float> discrepancieslm;
  for (int j = 0; j < m; j++)
    discrepancieslm.push_back(abs(dots(reslm.colors, ss.sets.at(j).points)));
  cout << "discrepancy llm: "
       << *max_element(discrepancieslm.begin(), discrepancieslm.end()) << endl;
  return 0;
}

// int main(int argc, char **argv) {
//   int r = time(NULL);
//   srand(r);
//   int n = 1024;
//   int d = 2;
//   int m = d * sqrt(n);
//   float p = .1;
//   SetSystem ss;
//   string ss_type = "grid";
//   bool save = false;
//   int c;
//   int centers = 3;
//   vector<int> algoList;
//   string filename = "";
//   float constant = 2.0;
//   // Options are detailed on the main.cpp documentation file
//   while ((c = getopt(argc, argv, "a:n:d:f:r:p:m:i:c:s")) != -1) {
//     switch (c) {
//     case 'a':
//       algoList = simple_tokenizer(optarg);
//       break;
//     case 'n':
//       n = stoi(optarg);
//       break;
//     case 'd':
//       d = stoi(optarg);
//       break;
//     case 'm':
//       m = stoi(optarg);
//       break;
//     case 'f':
//       ss_type = optarg;
//       break;
//     case 's':
//       save = true;
//       break;
//     case 'r':
//       r = stoi(optarg);
//       srand(r);
//       break;
//     case 'p':
//       p = stof(optarg);
//       break;
//     case 'i':
//       filename = optarg;
//       break;
//     case 'g':
//       centers = stoi(optarg);
//       break;
//     case 'c':
//       constant = stof(optarg);
//       break;
//     case '?':
//       if (optopt == 'a' || optopt == 'n' || optopt == 't' || optopt == 'd' ||
//           optopt == 'f' || optopt == 'r' || optopt == 'p' || optopt == 'm' ||
//           optopt == 'i' || optopt == 'c' || optopt == 'k' || optopt == 'g')
//         fprintf(stderr, "Option -%c requires an argument.\n", optopt);
//       else if (isprint(optopt))
//         fprintf(stderr, "Unknown option `-%c'.\n", optopt);
//       else
//         fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
//       return 1;
//     default:
//       abort();
//     }
//   }
//
//   // Generates set system
//   if (ss_type == "grid") {
//     ss = Grid(n, d);
//   } else if (ss_type == "grid_centers") {
//     ss = GridWithCenters(n, d, centers);
//   } else if (ss_type == "random_hs") {
//     m = n * log(n);
//     ss = RandomHyperplanes(n, d, m);
//   } else if (ss_type == "grid_graph") {
//     ss = GridGraph(sqrt(n), d);
//   } else if (ss_type == "linear_grid") {
//     ss = LinearGrid(n, d);
//   } else if (ss_type == "exponential_grid") {
//     ss = ExponentialGrid(n, d);
//   } else if (ss_type == "directed_grid") {
//     ss = DirectionalGrid(n, d);
//   } else if (ss_type == "random") {
//     ss = Random(n, d, m, p);
//   } else if (ss_type == "projective_plane") {
//     if (!isPrime(n)) {
//       fprintf(stderr, "Projective plane needs n to be prime, given: '%d'.\n",
//               n);
//       return 1;
//     }
//     ss = ProjectivePlane(n);
//   } else if (ss_type == "ERGGraph") {
//     ss = ERGGraph(n, d, p);
//   } else if (ss_type == "power_law") {
//     ss = PowerLaw(n, d, p, r);
//   } else if (ss_type == "concentric_circles") {
//     ss = ConcentricCircles(n, m);
//   } else if (ss_type == "file") {
//     ss = SetSystem(filename);
//   } else {
//     fprintf(stderr, "Unknown set system, given: '%s'.\n", ss_type.c_str());
//     return 1;
//   }
//   m = ss.sets.size();
//   n = ss.points.size();
//   vector<float> cst;
//   for (int j = 0; j < m; j++)
//     cst.push_back(sqrt(sumSet(ss.sets.at(j))));
//   Coloring res = lm(ss, cst);
//   for (int i = 0; i < n; i++) {
//     if (static_cast<float>(rand()) / static_cast<float>(RAND_MAX) >
//         abs(res.colors.at(i))) {
//       res.colors.at(i) = (res.colors.at(i) < 0) ? 1 : -1;
//     } else {
//       res.colors.at(i) = (res.colors.at(i) < 0) ? -1 : 1;
//     }
//   }
//   vector<float> discrepancies;
//   for (int j = 0; j < m; j++)
//     discrepancies.push_back(abs(dots(res.colors, ss.sets.at(j).points)));
//   cout << "discrepancy: "
//        << *max_element(discrepancies.begin(), discrepancies.end()) << endl;
//   if (save) {
//     ofstream MyFile("results.csv", std::ios_base::app);
//     for (int i = 0; i < n; i++) {
//       MyFile << res.colors.at(i) << ";";
//     }
//     MyFile << *max_element(discrepancies.begin(), discrepancies.end()) <<
//     endl;
//   }
//   return 0;
// }
