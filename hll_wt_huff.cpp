#include <algorithm>
#include <bits/stdc++.h>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
#include <set>
#include <vector>
using namespace std;
using namespace sdsl;

void kmers(std::string seg, int k);
int b = 4;
// Const for b = 4
// stdl hash struct
hash<string> h;
// Total buckets
double m = pow(2, b);
// Const
double am = (0.7213 / (1 + 1.079 / m));
// Return an int that cointains only the first b bytes
unsigned long long first_bb(unsigned long long x) { return (x >> (64 - b)); }
// Return an int equal to x minus the first b bytes
unsigned long long resto(unsigned long long x) { return (x << b); }
int_vector<> v1(m, 0, 8 * b);
std::set<std::string> Real;

vector<int> unionHLL(vector<int> v, vector<int> v2) {
  vector<int> f(m, 0);
  for (auto i : f) {
    f[i] = max(v[i], v2[i]);
  }
  return f;
}

// Returns the cardinallity of the intersection
//  |A|+|B| - |A U B|
unsigned long long intersection(unsigned long long A, unsigned long long B,
                                unsigned long long U) {
  return (A + B - U);
}

unsigned long long jaccard(unsigned long long I, unsigned long long U) {
  return (I / U);
}

void procesar(string M) {
  // Apply the hash function to the string i
  unsigned long long x = h(M);
  // Index
  unsigned long long j = first_bb(x);
  // x - j bytes
  unsigned long long w = resto(x);
  // Save the maximum in the bucket j
  v1[j] = max((int)v1[j], __builtin_clz(w) + 1);
}

unsigned long long hll() {
  // Raw Estimate
  double E;

  double sum = 0.0;
  wt_huff_int<rrr_vector<64>> wt;
  construct_im(wt, v1);
  // v1.resize(0);

  for (int i = 0; i < m; ++i) {
    const double max = (double)wt[i];
    const double term = pow(2.0, -max);
    sum += term;
  }
  const double harmonic_mean = m * (1.0 / sum);
  const double estimate = am * m * harmonic_mean;
  assert(estimate >= 0.0);
  E = estimate;

  double Eaux;
  // Small range correction
  if (E <= (5.0 / 2.0) * m) {

    double V = 0;
    for (int i = 0; i < m; ++i)
      if (v1[i] == 0)
        V++;
    if (V != 0) {
      Eaux = m * log2f((double)m / V);
    } else
      Eaux = E;
  } else if (E <=
             (1.0 / 30.0) * pow(2, 32)) { // Intermediate range => no correction
    Eaux = E;
  } else {
    Eaux = E;
  }

  // v1.clear();
  //  Return the estimate cardinality
  return Eaux;
}

int main() {
  int j = 0;
  std::cout << "Ingrese nombre del archivo: " << std::endl;

  std::string filename, line;
  std::cin >> filename;
  int S = 0, x = 31;
  std::ifstream newFile(filename);
  std::vector<std::string> k1;
  if (newFile.is_open()) {
    while (newFile.peek() != EOF) {
      getline(newFile, line);
      if (line.size() >= x)
        x = line.size();
      if (line != "")
        kmers(line, 31);
      k1.clear();
      S++;
    }
    newFile.close();
  }
  cout << "respuesta: " << hll() << endl;
  return 0;
}

void kmers(std::string seg, int k) {
  int size = seg.size();
  int kmer_size = size - k + 1;
  for (int i = 0; i < kmer_size; ++i) {
    // Inserting the kmer directly in to the sketch
    // Edit function to ignore N and >'lines
    procesar(seg.substr(i, k));
  }
}
