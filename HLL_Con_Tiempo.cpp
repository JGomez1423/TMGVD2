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
#include <sdsl/bit_vectors.hpp>
#include <sdsl/vectors.hpp>
#include <sdsl/wavelet_trees.hpp>
using namespace std;
using namespace sdsl;
uint64_t len = sizeof(long long);
int d = 8;
int w = 24;
int kcount=0;
//vector<vector<int> > C( d , vector<int>(w) ); /* Vector C to Countmin-CU */
//vector<int> murseed( d + 1 );  /* Vector of seed to murmushHash*/ 
int flag =0;

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


int_vector<> unionHLL_int(int_vector<> v, int_vector<> v2) {
  int_vector<> f(m, 0);
  for (int i=0; i<m;i++) {
    f[i] = max(v[i], v2[i]);
  }
  return f;
}
int_vector<> unionHLL_enc(enc_vector<> v, enc_vector<> v2) {
  int_vector<> f(m, 0);
  for (int i=0; i<m;i++) {
    f[i] = max(v[i], v2[i]);
  }
  return f;
}
int_vector<> unionHLL_wm(wm_int<rrr_vector<64>> v, wm_int<rrr_vector<64>> v2) {
  int_vector<> f(m, 0);
  for (int i=0; i<m;i++) {
    f[i] = max(v[i], v2[i]);
  }
  return f;
}
int_vector<> unionHLL_wt(wt_int<sd_vector<>> v, wt_int<sd_vector<>> v2) {
  int_vector<> f(m, 0);
  for (int i=0; i<m;i++) {
    f[i] = max(v[i], v2[i]);
  }
  return f;
}

int_vector<> unionHLL_whuff(wt_huff_int<rrr_vector<64>> v, wt_huff_int<rrr_vector<64>>  v2) {
  int_vector<> f(m, 0);
  for (int i=0; i<m;i++) {
    f[i] = max(v[i], v2[i]);
  } 
  return f;
}
//vector<int> v1(m, 0);
sdsl::int_vector<> v1(m,0);
sdsl::int_vector<> vunion(m,0);
enc_vector<sdsl::coder::elias_delta, 128> ev;
wm_int<rrr_vector<64>> wm;
wt_int<sd_vector<>> wt;
wt_huff_int<rrr_vector<64>> whuff;

enc_vector<sdsl::coder::elias_delta, 128> evunion;
wm_int<rrr_vector<64>> wmunion;
wt_int<sd_vector<>> wtunion;
wt_huff_int<rrr_vector<64>> whuffunion;
set<string> Real;
map<int,int> Frecuencia;
void procesar(unsigned long long x, int_vector<> &v) {
    // Index
    unsigned long long j = first_bb(x);
    // x - j bytes
    unsigned long long w = resto(x);
    // Save the maximum in the bucket j
    v[j] = max((int)v[j], __builtin_clz(w) + 1);//v1[j], );

}

double Entropia(){
    double H0_real = 0; 
    int M = 0;
    for(int i = 0; i < v1.size() ; i++){
        Frecuencia[v1[i]]++;
    }
    for(auto i : Frecuencia){
        M += i.second;
    }
    cout << "M: " << M << endl;
    for(auto i: Frecuencia){
        H0_real+= (double)(i.second/(double)M) * (double)log((double)(M/i.second));
    }
    return H0_real;
}
unsigned long long hll(int_vector<> v) {
    // Raw Estimate
    double E;
    double sum = 0.0;

    for (int i = 0; i < m; ++i) {
        const double max = (double)v[i];
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
        if (v[i] == 0)
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
unsigned long long hll_enc() {
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
unsigned long long hll_wm() {
    // Raw Estimate
    double E;
    double sum = 0.0;

    for (int i = 0; i < m; ++i) {
        const double max = (double)wm[i];
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
        if (wm[i] == 0)
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


unsigned long long hll_wt() {
    // Raw Estimate
    double E;
    double sum = 0.0;

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
        if (wt[i] == 0)
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
unsigned long long hll_whuff() {
    // Raw Estimate
    double E;
    double sum = 0.0;

    for (int i = 0; i < m; ++i) {
        const double max = (double)whuff[i];
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
        if (whuff[i] == 0)
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
    flag++;
    srand(time(NULL));
    //Inicializedseed()
    ifstream newFil(argv[2]);    

    if (newFil.is_open()) {
        while (newFil.peek() != EOF) {
            getline(newFil, line);
            if (line.size() >= x)
                x = line.size();
            if (line != "") 
                kmers(line, 31);
            k1.clear();
            S++;
        }
        newFile.close();
    }
    ev=v1;
    construct_im(wm, v1);
    construct_im(wt, v1);
    construct_im(whuff, v1);
    cout<<v1<<endl;
    cout<<vunion<<endl;

    cout<<unionHLL_int(v1,vunion)<<endl;
    cout << "respuesta real: "<<Frecuencia.size()<<endl;
    cout << "respuesta con int_vector: " << hll(v1) << endl;
    cout << "respuesta con int_vector: " << hll(vunion) << endl;
    cout << "respuesta con enc: " << hll_enc() << endl;
    cout << "respuesta con wm: "<<hll_wm()<<endl;
    cout << "respuesta con wt: "<<hll_wt()<<endl;
    cout << "respuesta con w_huff: "<<hll_whuff()<<endl;
    cout<<"Tamano del int_vector: "<<size_in_bytes(v1)<<endl;
    cout<<"Tamano del enc_vector: "<<size_in_bytes(ev)<<endl;
    cout<<"Tamano del wm_int: "<<size_in_bytes(wm)<<endl;
    cout<<"Tamano del wt_int: "<<size_in_bytes(wt)<<endl;
    cout<<"Tamano del wt_huff: "<<size_in_bytes(whuff)<<endl;
    cout<<"--------Analisis de Tiempo--------"<<endl;
    double tiempo_int=0,tiempo_enc=0, tiempo_wm=0, tiempo_wt =0, tiempo_w_huff=0;
    evunion=vunion;
    construct_im(wmunion, vunion);
    construct_im(wtunion, vunion);
    construct_im(whuffunion, vunion);

    double nano= 1000000000;
    clock_t start=clock();
    unionHLL_int(v1, vunion);
    tiempo_int=((double)clock()-start)/CLOCKS_PER_SEC;
    cout<<"Tiempo int_vector: "<<tiempo_int<<endl;

    start=clock();
    unionHLL_enc(ev, evunion);
    tiempo_enc=((double)clock()-start)/CLOCKS_PER_SEC;
    cout<<"Tiempo enc_vector: "<<tiempo_int<<endl;

    start=clock();
    unionHLL_wm(wm, wmunion);
    tiempo_wm=((double)clock()-start)/CLOCKS_PER_SEC;
    cout<<"Tiempo wm_int: "<<tiempo_int<<endl;

    start=clock();
    unionHLL_wt(wt, wtunion);
    tiempo_wt=((double)clock()-start)/CLOCKS_PER_SEC;
    cout<<"Tiempo wt_int: "<<tiempo_int<<endl;

    start=clock();
    unionHLL_whuff(whuff, whuffunion);
    tiempo_w_huff=((double)clock()-start)/CLOCKS_PER_SEC;
    cout<<"Tiempo wt_huff_int: "<<tiempo_int<<endl;

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
        if(flag==0){
            procesar(x,v1);
        }else procesar(x, vunion);
        //Frecuencia[Km]++;
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
