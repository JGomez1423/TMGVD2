#include <iostream>
#include <vector>
#include <fstream>     
#include <string>
#include <math.h>
#include <algorithm>
#include <map>
#include <cstdint>
#include <set>
#include <climits>
#include "sdsl/vectors.hpp"
using namespace std;
using namespace sdsl;
uint64_t len = sizeof(long long);
int d = 8;
int w = 24;
int kcount=0;
//vector<vector<int> > C( d , vector<int>(w) ); /* Vector C to Countmin-CU */
//vector<int> murseed( d + 1 );  /* Vector of seed to murmushHash*/ 


void kmers(string seg, int k);
//void CountMin_Cu(unsigned long long x);

int b = 4;          // Const for b = 4
hash<string> h;       // stdl hash struct
double m = pow(2, b); // Total buckets


//Const
double am = (0.7213/(1 + 1.079/m));
// Return an int that cointains only the first b bytes
unsigned long long first_bb(unsigned long long x) { return (x >> (64 - b)); }
// Return an int equal to x minus the first b bytes
unsigned long long resto(unsigned long long x) { return (x << b); }

//vector<int> v1(m, 0);
sdsl::int_vector<> v1(m,0);
enc_vector<sdsl::coder::elias_delta, 128> ev;
set<string> Real;
map<string,int> Frecuencia;
void procesar(unsigned long long x) {
    // Index
    unsigned long long j = first_bb(x);
    // x - j bytes
    unsigned long long w = resto(x);
    // Save the maximum in the bucket j
    v1[j] = max((int)v1[j], __builtin_clz(w) + 1);//v1[j], );
}

unsigned long long hll() {
    // Raw Estimate
    double E;
    double sum = 0.0;

    for (int i = 0; i < m; ++i) {
        const double max = (double)v1[i];
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
    } else if (E <= (1.0 / 30.0) * pow(2, 32)) { // Intermediate range => no correction
        Eaux = E;
    }else{
        Eaux = E;
    }
    //v1.clear();
    // Return the estimate cardinality
    return Eaux;
}
unsigned long long hll2() {
    // Raw Estimate
    double E;
    double sum = 0.0;

    for (int i = 0; i < m; ++i) {
        const double max = (double)ev[i];
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
        if (ev[i] == 0)
            V++;
        if (V != 0) {
            Eaux = m * log2f((double)m / V);
        } else
            Eaux = E;
    } else if (E <= (1.0 / 30.0) * pow(2, 32)) { // Intermediate range => no correction
        Eaux = E;
    }else{
        Eaux = E;
    }
    //v1.clear();
    // Return the estimate cardinality
    return Eaux;
}





int main( int argc, char *argv[] ){


    srand(time(NULL));
    //Inicializedseed();
    int j = 0;
    string line;
    int S = 0, x = 31;
    ifstream newFile(argv[1]);    
    vector<string> k1;

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
    //for(map<string,int>::iterator it = Frecuencia.begin() ; it != Frecuencia.end(); ++it)
        //cout << it->first << " => " << it->second << '\n';
    
    ev=v1;
    cout << "respuesta: " << hll() << endl;
    cout << "respuesta con enc: " << hll2() << endl;
    cout<<v1<<endl;
    cout<<ev<<endl;
    cout<<size_in_bytes(v1)<<endl;
    cout<<size_in_bytes(ev)<<endl;

    return 0;
}


void kmers(string seg, int k) {
    int size = seg.size();
    int kmer_size = size - k + 1;
    string Km = "";
    for (int i = 0; i < kmer_size; ++i){
        // Inserting the kmer directly in to the sketch
        // Edit function to ignore N and >'lines
            // Apply the hash function to the string i
        Km = seg.substr(i,k);
        unsigned long long x = h(Km);
        procesar(x);
        Frecuencia[Km]++;
        kcount++;
        //CountMin_Cu(x);
    }
}
/*void CountMin_Cu(unsigned long long x){
    int minaux = INT_MAX;    
    for(int i = 0 ; i < d ; i++){
        unsigned long long  aux = MurmurHash2_x64_64A(&x, len ,murseed[i]) % w;
        if(C[i][aux] < minaux)
            minaux = C[i][aux];
    }
    for(int i = 0 ; i < d ; i++){
        unsigned long long  aux = MurmurHash2_x64_64A(&x, len ,murseed[i]) % w;
        if(C[i][aux] == minaux)
            C[i][aux]++;      
    }
}*/
