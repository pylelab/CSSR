/*
 *LinearPartition.cpp*
 The main code for LinearPartition: Linear-Time Approximation of 
                                    RNA Folding Partition Function 
                                    and Base Pairing Probabilities

 author: He Zhang
 created by: 03/2019
 modified by Chengxin Zhang 08/2021
*/

#include <fstream>
#include <iostream>
#include <sys/time.h>
#include <stack>
#include <tuple>
#include <cassert>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <map>
#include <stdio.h> 
#include <set>

#include "LinearPartition.h"
#include "Utils/utility.h"
#include "Utils/utility_v.h"

using namespace std;

/* LinearPartition.h inline function start */
inline pf_type Fast_LogExpPlusOne(pf_type x){
  
    // Bounds for tolerance of 7.05e-06: (0, 11.8625)
    // Approximating interval: (0, 0.661537) --> ((T(-0.0065591595)*x+T(0.1276442762))*x+T(0.4996554598))*x+T(0.6931542306);
    // Approximating interval: (0.661537, 1.63202) --> ((T(-0.0155157557)*x+T(0.1446775699))*x+T(0.4882939746))*x+T(0.6958092989);
    // Approximating interval: (1.63202, 2.49126) --> ((T(-0.0128909247)*x+T(0.1301028251))*x+T(0.5150398748))*x+T(0.6795585882);
    // Approximating interval: (2.49126, 3.37925) --> ((T(-0.0072142647)*x+T(0.0877540853))*x+T(0.6208708362))*x+T(0.5909675829);
    // Approximating interval: (3.37925, 4.42617) --> ((T(-0.0031455354)*x+T(0.0467229449))*x+T(0.7592532310))*x+T(0.4348794399);
    // Approximating interval: (4.42617, 5.78907) --> ((T(-0.0010110698)*x+T(0.0185943421))*x+T(0.8831730747))*x+T(0.2523695427);
    // Approximating interval: (5.78907, 7.81627) --> ((T(-0.0001962780)*x+T(0.0046084408))*x+T(0.9634431978))*x+T(0.0983148903);
    // Approximating interval: (7.81627, 11.8625) --> ((T(-0.0000113994)*x+T(0.0003734731))*x+T(0.9959107193))*x+T(0.0149855051);
    // 8 polynomials needed.
    
    assert(pf_type(0.0000000000) <= x && x <= pf_type(11.8624794162) && "Argument out-of-range.");
    if (x < pf_type(3.3792499610))
    {
        if (x < pf_type(1.6320158198))
        {
            if (x < pf_type(0.6615367791))
                return ((pf_type(-0.0065591595)*x+pf_type(0.1276442762))*x+pf_type(0.4996554598))*x+pf_type(0.6931542306);
            return ((pf_type(-0.0155157557)*x+pf_type(0.1446775699))*x+pf_type(0.4882939746))*x+pf_type(0.6958092989);
        }
        if (x < pf_type(2.4912588184))
            return ((pf_type(-0.0128909247)*x+pf_type(0.1301028251))*x+pf_type(0.5150398748))*x+pf_type(0.6795585882);
        return ((pf_type(-0.0072142647)*x+pf_type(0.0877540853))*x+pf_type(0.6208708362))*x+pf_type(0.5909675829);
    }
    if (x < pf_type(5.7890710412))
    {
        if (x < pf_type(4.4261691294))
            return ((pf_type(-0.0031455354)*x+pf_type(0.0467229449))*x+pf_type(0.7592532310))*x+pf_type(0.4348794399);
        return ((pf_type(-0.0010110698)*x+pf_type(0.0185943421))*x+pf_type(0.8831730747))*x+pf_type(0.2523695427);
    }
    if (x < pf_type(7.8162726752))
        return ((pf_type(-0.0001962780)*x+pf_type(0.0046084408))*x+pf_type(0.9634431978))*x+pf_type(0.0983148903);
    return ((pf_type(-0.0000113994)*x+pf_type(0.0003734731))*x+pf_type(0.9959107193))*x+pf_type(0.0149855051);
}

inline void Fast_LogPlusEquals (pf_type &x, pf_type y)
{
    if (x < y) std::swap (x, y);
    if (y > pf_type(NEG_INF/2) && x-y < pf_type(11.8624794162))
        x = Fast_LogExpPlusOne(x-y) + y;
}

inline pf_type Fast_Exp(pf_type x)
{
    // Bounds for tolerance of 4.96e-05: (-9.91152, 0)
    // Approximating interval: (-9.91152, -5.86228) --> ((T(0.0000803850)*x+T(0.0021627428))*x+T(0.0194708555))*x+T(0.0588080014);
    // Approximating interval: (-5.86228, -3.83966) --> ((T(0.0013889414)*x+T(0.0244676474))*x+T(0.1471290604))*x+T(0.3042757740);
    // Approximating interval: (-3.83966, -2.4915) --> ((T(0.0072335607)*x+T(0.0906002677))*x+T(0.3983111356))*x+T(0.6245959221);
    // Approximating interval: (-2.4915, -1.48054) --> ((T(0.0232410351)*x+T(0.2085645908))*x+T(0.6906367911))*x+T(0.8682322329);
    // Approximating interval: (-1.48054, -0.672505) --> ((T(0.0573782771)*x+T(0.3580258429))*x+T(0.9121133217))*x+T(0.9793091728);
    // Approximating interval: (-0.672505, -3.9145e-11) --> ((T(0.1199175927)*x+T(0.4815668234))*x+T(0.9975991939))*x+T(0.9999505077);
    // 6 polynomials needed.
    
    if (x < pf_type(-2.4915033807))
    {
        if (x < pf_type(-5.8622823336))
        {
            if (x < pf_type(-9.91152))
                return pf_type(0);
            return ((pf_type(0.0000803850)*x+pf_type(0.0021627428))*x+pf_type(0.0194708555))*x+pf_type(0.0588080014);
        }
        if (x < pf_type(-3.8396630909))
            return ((pf_type(0.0013889414)*x+pf_type(0.0244676474))*x+pf_type(0.1471290604))*x+pf_type(0.3042757740);
        return ((pf_type(0.0072335607)*x+pf_type(0.0906002677))*x+pf_type(0.3983111356))*x+pf_type(0.6245959221);
    }
    if (x < pf_type(-0.6725053211))
    {
        if (x < pf_type(-1.4805375919))
            return ((pf_type(0.0232410351)*x+pf_type(0.2085645908))*x+pf_type(0.6906367911))*x+pf_type(0.8682322329);
        return ((pf_type(0.0573782771)*x+pf_type(0.3580258429))*x+pf_type(0.9121133217))*x+pf_type(0.9793091728);
    }
    if (x < pf_type(0))
        return ((pf_type(0.1199175927)*x+pf_type(0.4815668234))*x+pf_type(0.9975991939))*x+pf_type(0.9999505077);
    return (x > pf_type(46.052) ? pf_type(1e20) : expf(x));
}
/* LinearPartition.h inline function end */

/* bpp.cpp
 The main code for base pair probability calculation.

 author: He Zhang
 created by: 04/2019
*/

void BeamCKYParser::output_to_file(string file_name, const char * type) {
    if(!file_name.empty()) {
        printf("Outputing base pairing probability matrix to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        int turn = no_sharp_turn?3:0;
        for (int i = 1; i <= seq_length; i++) {
            for (int j = i + turn + 1; j <= seq_length; j++) {
                pair<int, int> key = make_pair(i,j);
                auto got = Pij.find(key);
                if (got != Pij.end()){
                    fprintf(fptr, "%d %d %.4e\n", i, j, got->second);
                }
            }
        }
        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }

    return;
}

void BeamCKYParser::output_to_file_MEA_threshknot_bpseq(string file_name, const char * type, map<int,int>& pairs, string & seq) {

    int i,j;
    char nuc;
    if(!file_name.empty()) {
        printf("Outputing base pairs in bpseq format to %s...\n", file_name.c_str()); 
        FILE *fptr = fopen(file_name.c_str(), type); 
        if (fptr == NULL) { 
            printf("Could not open file!\n"); 
            return; 
        }

        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            fprintf(fptr, "%d %c %d\n", i, nuc, j);
        }

        fprintf(fptr, "\n");
        fclose(fptr); 
        printf("Done!\n"); 
    }
    else{
        for (int i = 1; i <= seq_length; i++) {
            if (pairs.find(i) != pairs.end()){
                j = pairs[i];
            }
            else{
                j = 0;
            }
            nuc = seq[i-1];
            printf("%d %c %d\n", i, nuc, j);
        }
        printf("\n");
    }

}

void BeamCKYParser::cal_PairProb(State& viterbi,
     vector<pair<float,pair<int,int> > >&RR_list)
{

    for(int j=0; j<seq_length; j++){
        for(auto &item : bestP[j]){
            int i = item.first;
            State state = item.second;
            
            pf_type temp_prob_inside = state.alpha + state.beta - viterbi.alpha;
            if (temp_prob_inside > pf_type(-9.91152)) {
                pf_type prob = Fast_Exp(temp_prob_inside);
                if(prob > pf_type(1.0)) prob = pf_type(1.0);
                if(prob < pf_type(bpp_cutoff)) continue;
                Pij[make_pair(i+1, j+1)] = prob;
                RR_list.push_back(make_pair(prob,make_pair(i,j)));
            }
        }
    }

    // -o mode: output to a single file with user specified name;
    // bpp matrices for different sequences are separated with empty lines
    if (!bpp_file.empty()){
        output_to_file(bpp_file, "a");
    } 

    // -prefix mode: output to multiple files with user specified prefix;
    else if (!bpp_file_index.empty()) {
        output_to_file(bpp_file_index, "w");
    }
    return;
}


string BeamCKYParser::back_trace(const int i, const int j, const vector<vector<int> >& back_pointer){

    if (i>j) return "";
    if (back_pointer[i][j] == -1){
        if (i == j) return ".";
        else return "." + back_trace(i+1,j, back_pointer);
    }else if (back_pointer[i][j] != 0){
        int k = back_pointer[i][j];
        assert(k + 1 > 0 && k + 1 <= seq_length);
        string temp;
        if (k == j) temp = "";
        else temp = back_trace(k+1,j, back_pointer);
        return "(" + back_trace(i+1,k-1, back_pointer) + ")" + temp;
    }
    assert(false);
    return "";
}

map<int, int> BeamCKYParser::get_pairs(string & structure){
    map<int, int> pairs;
    stack<int> s;
    int index = 1;
    int pre_index = 0;
    for (auto & elem : structure){
        if (elem == '(') s.push(index);
        else if(elem == ')'){
            pre_index = s.top();
            pairs[pre_index] = index;
            pairs[index] = pre_index;
            s.pop();
        }
        index++;
    }
    return pairs;
}

void BeamCKYParser::ThreshKnot(string & seq){
    
    map<int, pf_type> rowprob;
    vector<tuple<int, int, pf_type> > prob_list;

    map<int, int> pairs;
    set<int> visited;

    for(auto& pij : Pij){
        auto i = pij.first.first; //index starts from 1
        auto j = pij.first.second; 
        auto score = pij.second;

        if (score < threshknot_threshold) continue;

        prob_list.push_back(make_tuple(i,j,score));

        rowprob[i] = max(rowprob[i], score);
        rowprob[j] = max(rowprob[j], score);

    }

    for(auto& elem : prob_list){

        auto i = std::get<0>(elem);
        auto j = std::get<1>(elem);
        auto score =  std::get<2>(elem);

        if (score == rowprob[i] && score == rowprob[j]){

            if ((visited.find(i) != visited.end()) || (visited.find(j) != visited.end())) continue;
            visited.insert(i);
            visited.insert(j);

            pairs[i] = j;
            pairs[j] = i;
        }
    }

    if (is_verbose) fprintf(stderr, "%s\n", seq.c_str());
    output_to_file_MEA_threshknot_bpseq(threshknot_file_index, "w", pairs, seq);
}

void BeamCKYParser::PairProb_MEA(string & seq) {
    
    vector<vector<pf_type> > OPT;
    OPT.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) OPT[i].resize(seq_length);

    vector<vector<pf_type>> P;
    P.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) P[i].resize(seq_length);

    vector<vector<int> > back_pointer;
    back_pointer.resize(seq_length);

    for (int i = 0; i < seq_length; ++i) back_pointer[i].resize(seq_length);

    vector<vector<int>> paired;
    paired.resize(seq_length);

    vector<pf_type> Q;
    for (int i = 0; i < seq_length; ++i) Q.push_back(pf_type(1.0));

    for(auto& pij : Pij){
        auto i = pij.first.first-1;
        auto j = pij.first.second-1;
        auto score = pij.second;

        P[i][j] = score;

        paired[i].push_back(j);
        Q[i] -= score;
        Q[j] -= score;
    }

    for (int i = 0; i < seq_length; ++i) std::sort (paired[i].begin(), paired[i].end());
    for (int l = 0; l< seq_length; l++){
        for (int i = 0; i<seq_length - l; i++){
            int j = i + l;
            if (i == j){
                OPT[i][j] = Q[i];
                back_pointer[i][j] = -1;
                continue;
            }
            OPT[i][j] = OPT[i][i] + OPT[i+1][j];
            back_pointer[i][j] = -1;
            for (int k : paired[i]){
                if (k>j) break;
                pf_type temp_OPT_k1_j;
                if (k<j) temp_OPT_k1_j = OPT[k+1][j];
                else temp_OPT_k1_j = pf_type(0.);
                auto temp_score = 2 * gamma * P[i][k] + OPT[i+1][k-1] + temp_OPT_k1_j;
                if (OPT[i][j] < temp_score){
                    OPT[i][j] = temp_score;
                    back_pointer[i][j] = k;
                }
            }
        }
    }

    auto structure = back_trace(0,seq_length-1, back_pointer);

    if (!bpseq){
        if(!mea_file_index.empty()) {
            FILE *fptr = fopen(mea_file_index.c_str(), "w"); 
            if (fptr == NULL) { 
                printf("Could not open file!\n"); 
                return; 
            }
            fprintf(fptr, "%s\n", seq.c_str());
            fprintf(fptr, "%s\n\n", structure.c_str());
        }

        else{
            printf("%s\n", seq.c_str());
            printf("%s\n\n", structure.c_str());
        }
    }

    else{
        auto pairs = get_pairs(structure);
        output_to_file_MEA_threshknot_bpseq(mea_file_index, "w", pairs, seq);
    }
}


void BeamCKYParser::outside(vector<int> next_pair[]){
      
    struct timeval bpp_starttime, bpp_endtime;
    gettimeofday(&bpp_starttime, NULL);

    bestC[seq_length-1].beta = 0.0;

    // from right to left
    value_type newscore;
    for(int j = seq_length-1; j > 0; --j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
            Fast_LogPlusEquals(beamstepC.beta, (bestC[j+1].beta));
                    
#else
            newscore = score_external_unpaired(j+1, j+1);
            Fast_LogPlusEquals(beamstepC.beta, bestC[j+1].beta + newscore);
#endif
            }
        }
    
        // beam of M
        {
            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta);
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(state.beta, bestM[j+1][i].beta + newscore);
#endif
                }
            }
        }

        // beam of M2
        {
            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                {
                    for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int q = next_pair[nucp][j];
                        if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lpv
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta);
#else
                            newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                            Fast_LogPlusEquals(state.beta, bestMulti[q][p].beta + newscore);
#endif
                        }
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(state.beta, beamstepM[i].beta);
            }
        }

        // beam of P
        {  
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                if (i >0 && j<seq_length-1) {
#ifndef lpv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lpv
                                newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);
                                // SHAPE for Vienna only
                                if (use_shape)
                                {
                                    newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                }


                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            } else {
                                // single branch
#ifdef lpv
                                newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore/kT);
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed + 
                                        score_single_without_junctionB(p, q, i, j, nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(state.beta, bestP[q][p].beta + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lpv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore/kT);
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepM[i].beta + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = newscore/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = newscore;
#endif
                    pf_type m1_plus_P_alpha = state.alpha + m1_alpha;
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(state.beta, (beamstepM2[newi].beta + m_state.alpha + m1_alpha));
                        Fast_LogPlusEquals(m_state.beta, (beamstepM2[newi].beta + m1_plus_P_alpha));
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                                 nucj, nucj1, seq_length);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore/kT;

#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1, nucj, nucj1, seq_length);
                        pf_type external_paired_alpha_plus_beamstepC_beta = beamstepC.beta + newscore;
#endif
                        Fast_LogPlusEquals(bestC[k].beta, state.alpha + external_paired_alpha_plus_beamstepC_beta);
                        Fast_LogPlusEquals(state.beta, bestC[k].alpha + external_paired_alpha_plus_beamstepC_beta);
                    } else {
                        // value_type newscore;
#ifdef lpv
                        newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, (beamstepC.beta + newscore/kT));
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(state.beta, beamstepC.beta + newscore);
#endif
                    }
                }
            }
        }

        // beam of Multi
        {
            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lpv
                        Fast_LogPlusEquals(state.beta, (bestMulti[jnext][i].beta));
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(state.beta, bestMulti[jnext][i].beta + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore/kT);
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(state.beta, beamstepP[i].beta + newscore);
#endif
                }
            }
        }
    }  // end of for-loo j

    gettimeofday(&bpp_endtime, NULL);
    double bpp_elapsed_time = bpp_endtime.tv_sec - bpp_starttime.tv_sec + (bpp_endtime.tv_usec-bpp_starttime.tv_usec)/1000000.0;

    if(is_verbose) fprintf(stderr,"Base Pairing Probabilities Calculation Time: %.2f seconds.\n", bpp_elapsed_time);

    fflush(stdout);

    return;
}

/* bpp.cpp end */

#define SPECIAL_HP

unsigned long quickselect_partition(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper) {
    pf_type pivot = scores[upper].first;
    while (lower < upper) {
        while (scores[lower].first < pivot) ++lower;
        while (scores[upper].first > pivot) --upper;
        if (scores[lower].first == scores[upper].first) ++lower;
        else if (lower < upper) swap(scores[lower], scores[upper]);
    }
    return upper;
}

// in-place quick-select
pf_type quickselect(vector<pair<pf_type, int>>& scores, unsigned long lower, unsigned long upper, unsigned long k) {
    if ( lower == upper ) return scores[lower].first;
    unsigned long split = quickselect_partition(scores, lower, upper);
    unsigned long length = split - lower + 1;
    if (length == k) return scores[split].first;
    else if (k  < length) return quickselect(scores, lower, split-1, k);
    else return quickselect(scores, split+1, upper, k - length);
}


pf_type BeamCKYParser::beam_prune(std::unordered_map<int, State> &beamstep) {
    scores.clear();
    for (auto &item : beamstep) {
        int i = item.first;
        State &cand = item.second;
        int k = i - 1;
        pf_type newalpha = (k >= 0 ? bestC[k].alpha : pf_type(0.0)) + cand.alpha;
        scores.push_back(make_pair(newalpha, i));
    }
    if (scores.size() <= beam) return VALUE_MIN;
    pf_type threshold = quickselect(scores, 0, scores.size() - 1, scores.size() - beam);
    for (auto &p : scores) {
        if (p.first < threshold) beamstep.erase(p.second);
    }

    return threshold;
}


void BeamCKYParser::prepare(unsigned len) {
    seq_length = len;

    nucs = new int[seq_length];
    bestC = new State[seq_length];
    bestH = new unordered_map<int, State>[seq_length];
    bestP = new unordered_map<int, State>[seq_length];
    bestM = new unordered_map<int, State>[seq_length];
    bestM2 = new unordered_map<int, State>[seq_length];
    bestMulti = new unordered_map<int, State>[seq_length];
    
    scores.reserve(seq_length);
}

void BeamCKYParser::postprocess() {

    delete[] bestC;  
    delete[] bestH;  
    delete[] bestP;  
    delete[] bestM;  
    delete[] bestM2;  
    delete[] bestMulti;  

    delete[] nucs;  
}

void BeamCKYParser::parse(string& seq,
    vector<pair<float,pair<int,int> > >&RR_list)
{
      
    struct timeval parse_starttime, parse_endtime;

    gettimeofday(&parse_starttime, NULL);

    prepare(static_cast<unsigned>(seq.length()));

    for (int i = 0; i < seq_length; ++i)
        nucs[i] = GET_ACGU_NUM(seq[i]);

    vector<int> next_pair[NOTON];
    {
        for (int nuci = 0; nuci < NOTON; ++nuci) {
            // next_pair
            next_pair[nuci].resize(seq_length, -1);
            int next = -1;
            for (int j = seq_length-1; j >=0; --j) {
                next_pair[nuci][j] = next;
                if (_allowed_pairs[nuci][nucs[j]]) next = j;
            }
        }
    }

#ifdef SPECIAL_HP
#ifdef lpv
    v_init_tetra_hex_tri(seq, seq_length, if_tetraloops, if_hexaloops, if_triloops);
#endif
#endif

#ifdef lpv
        if(seq_length > 0) bestC[0].alpha = 0.0;
        if(seq_length > 1) bestC[1].alpha = 0.0;
#else
        if(seq_length > 0) Fast_LogPlusEquals(bestC[0].alpha, score_external_unpaired(0, 0));
        if(seq_length > 1) Fast_LogPlusEquals(bestC[1].alpha, score_external_unpaired(0, 1));
#endif

    value_type newscore;
    for(int j = 0; j < seq_length; ++j) {
        int nucj = nucs[j];
        int nucj1 = (j+1) < seq_length ? nucs[j+1] : -1;

        unordered_map<int, State>& beamstepH = bestH[j];
        unordered_map<int, State>& beamstepMulti = bestMulti[j];
        unordered_map<int, State>& beamstepP = bestP[j];
        unordered_map<int, State>& beamstepM2 = bestM2[j];
        unordered_map<int, State>& beamstepM = bestM[j];
        State& beamstepC = bestC[j];

        // beam of H
        {
            if (beam > 0 && beamstepH.size() > beam) beam_prune(beamstepH);

            {
                // for nucj put H(j, j_next) into H[j_next]
                int jnext = next_pair[nucj][j];
                if (no_sharp_turn) while (jnext - j < 4 && jnext != -1) jnext = next_pair[nucj][jnext];
                if (jnext != -1) {
                    int nucjnext = nucs[jnext];
                    int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;
#ifdef lpv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-j-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[j];
                        else if (jnext-j-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[j];
                        else if (jnext-j-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[j];
#endif
                        newscore = - v_score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore/kT);
#else
                        newscore = score_hairpin(j, jnext, nucj, nucj1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][j].alpha, newscore);
#endif
                }
            }

            {
                // for every state h in H[j]
                //   1. extend h(i, j) to h(i, jnext)
                //   2. generate p(i, j)
                for (auto &item : beamstepH) {
                    int i = item.first;
                    State &state = item.second;
                    int nuci = nucs[i];
                    int jnext = next_pair[nuci][j];

                    if (jnext != -1) {
                        int nuci1 = (i + 1) < seq_length ? nucs[i + 1] : -1;
                        int nucjnext = nucs[jnext];
                        int nucjnext_1 = (jnext - 1) > -1 ? nucs[jnext - 1] : -1;

                        // 1. extend h(i, j) to h(i, jnext)=
#ifdef lpv
                        int tetra_hex_tri = -1;
#ifdef SPECIAL_HP
                        if (jnext-i-1 == 4) // 6:tetra
                            tetra_hex_tri = if_tetraloops[i];
                        else if (jnext-i-1 == 6) // 8:hexa
                            tetra_hex_tri = if_hexaloops[i];
                        else if (jnext-i-1 == 3) // 5:tri
                            tetra_hex_tri = if_triloops[i];
#endif
                        newscore = - v_score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext, tetra_hex_tri);
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore/kT);
#else
                        newscore = score_hairpin(i, jnext, nuci, nuci1, nucjnext_1, nucjnext);
                        Fast_LogPlusEquals(bestH[jnext][i].alpha, newscore);
#endif
                    }

                    // 2. generate p(i, j)
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha);
                }
            }
        }
        if (j == 0) continue;

        // beam of Multi
        {
            if (beam > 0 && beamstepMulti.size() > beam) beam_prune(beamstepMulti);

            for(auto& item : beamstepMulti) {
                int i = item.first;
                State& state = item.second;

                int nuci = nucs[i];
                int nuci1 = nucs[i+1];
                int jnext = next_pair[nuci][j];

                // 1. extend (i, j) to (i, jnext)
                {
                    if (jnext != -1) {
#ifdef lpv
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha);
#else
                        newscore = score_multi_unpaired(j, jnext - 1);
                        Fast_LogPlusEquals(bestMulti[jnext][i].alpha, state.alpha + newscore);
#endif
                    }
                }

                // 2. generate P (i, j)
                {
#ifdef lpv
                    newscore = - v_score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore/kT);
#else
                    newscore = score_multi(i, j, nuci, nuci1, nucs[j-1], nucj, seq_length);
                    Fast_LogPlusEquals(beamstepP[i].alpha, state.alpha + newscore);
#endif
                }
            }
        }

        // beam of P
        {   
            if (beam > 0 && beamstepP.size() > beam) beam_prune(beamstepP);

            // for every state in P[j]
            //   1. generate new helix/bulge
            //   2. M = P
            //   3. M2 = M + P
            //   4. C = C + P
            for(auto& item : beamstepP) {
                int i = item.first;
                State& state = item.second;
                int nuci = nucs[i];
                int nuci_1 = (i-1>-1) ? nucs[i-1] : -1;

                // 1. generate new helix / single_branch
                // new state is of shape p..i..j..q
                if (i >0 && j<seq_length-1) {
#ifndef lpv
                    value_type precomputed = score_junction_B(j, i, nucj, nucj1, nuci_1, nuci);
#endif
                    for (int p = i - 1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                        int nucp = nucs[p];
                        int nucp1 = nucs[p + 1]; 
                        int q = next_pair[nucp][j];
                        while (q != -1 && ((i - p) + (q - j) - 2 <= SINGLE_MAX_LEN)) {
                            int nucq = nucs[q];
                            int nucq_1 = nucs[q - 1];

                            if (p == i - 1 && q == j + 1) {
                                // helix
#ifdef lpv
                                newscore = -v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                             nuci_1, nuci, nucj, nucj1);

                                // SHAPE for Vienna only
                                if (use_shape)
                                {
                                    newscore += -(pseudo_energy_stack[p] + pseudo_energy_stack[i] + pseudo_energy_stack[j] + pseudo_energy_stack[q]);
                                }


                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
#else
                                newscore = score_helix(nucp, nucp1, nucq_1, nucq);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);

#endif
                            } else {
                                // single branch
#ifdef lpv
                                newscore = - v_score_single(p,q,i,j, nucp, nucp1, nucq_1, nucq,
                                                   nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore/kT);
#else
                                newscore = score_junction_B(p, q, nucp, nucp1, nucq_1, nucq) +
                                        precomputed +
                                        score_single_without_junctionB(p, q, i, j,
                                                                       nuci_1, nuci, nucj, nucj1);
                                Fast_LogPlusEquals(bestP[q][p].alpha, state.alpha + newscore);
#endif
                            }
                            q = next_pair[nucp][q];
                        }
                    }
                }

                // 2. M = P
                if(i > 0 && j < seq_length-1){
#ifdef lpv
                        newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore/kT);
#else
                        newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha + newscore);
#endif
                }

                // 3. M2 = M + P
                int k = i - 1;
                if ( k > 0 && !bestM[k].empty()) {
#ifdef lpv
                    newscore = - v_score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = state.alpha + newscore/kT;
#else
                    newscore = score_M1(i, j, j, nuci_1, nuci, nucj, nucj1, seq_length);
                    pf_type m1_alpha = state.alpha + newscore;
#endif
                    for (auto &m : bestM[k]) {
                        int newi = m.first;
                        State& m_state = m.second;
                        Fast_LogPlusEquals(beamstepM2[newi].alpha, m_state.alpha + m1_alpha);
                    }
                }

                // 4. C = C + P
                {
                    int k = i - 1;
                    if (k >= 0) {
                        State& prefix_C = bestC[k];
                        int nuck = nuci_1;
                        int nuck1 = nuci;
#ifdef lpv
                        newscore = - v_score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore/kT);      
#else
                        newscore = score_external_paired(k+1, j, nuck, nuck1,
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, prefix_C.alpha + state.alpha + newscore);
#endif
                    } else {
#ifdef lpv
                        newscore = - v_score_external_paired(0, j, -1, nucs[0],
                                                                 nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore/kT);       
#else
                        newscore = score_external_paired(0, j, -1, nucs[0],
                                                             nucj, nucj1, seq_length);
                        Fast_LogPlusEquals(beamstepC.alpha, state.alpha + newscore);
#endif
                    }
                }
            }
        }

        // beam of M2
        {
            if (beam > 0 && beamstepM2.size() > beam) beam_prune(beamstepM2);

            for(auto& item : beamstepM2) {
                int i = item.first;
                State& state = item.second;

                // 1. multi-loop
                for (int p = i-1; p >= std::max(i - SINGLE_MAX_LEN, 0); --p) {
                    int nucp = nucs[p];
                    int q = next_pair[nucp][j];
                    if (q != -1 && ((i - p - 1) <= SINGLE_MAX_LEN)) {
#ifdef lpv
                    Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha);      

#else
                    newscore = score_multi_unpaired(p+1, i-1) +
                                    score_multi_unpaired(j+1, q-1);
                    Fast_LogPlusEquals(bestMulti[q][p].alpha, state.alpha + newscore);      
#endif
                    }
                }

                // 2. M = M2
                Fast_LogPlusEquals(beamstepM[i].alpha, state.alpha);  
            }
        }

        // beam of M
        {
            if (beam > 0 && beamstepM.size() > beam) beam_prune(beamstepM);

            for(auto& item : beamstepM) {
                int i = item.first;
                State& state = item.second;
                if (j < seq_length-1) {
#ifdef lpv
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha); 
#else
                    newscore = score_multi_unpaired(j + 1, j + 1);
                    Fast_LogPlusEquals(bestM[j+1][i].alpha, state.alpha + newscore); 
#endif
                }
            }
        }

        // beam of C
        {
            // C = C + U
            if (j < seq_length-1) {
#ifdef lpv
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha); 
                    
#else
                newscore = score_external_unpaired(j+1, j+1);
                Fast_LogPlusEquals(bestC[j+1].alpha, beamstepC.alpha + newscore); 
#endif
            }
        }
    }  // end of for-loo j

    State& viterbi = bestC[seq_length-1];

    gettimeofday(&parse_endtime, NULL);
    double parse_elapsed_time = parse_endtime.tv_sec - parse_starttime.tv_sec + (parse_endtime.tv_usec-parse_starttime.tv_usec)/1000000.0;

    // unsigned long nos_tot = nos_H + nos_P + nos_M2 + nos_Multi + nos_M + nos_C;

#ifdef lpv
    if (is_verbose) fprintf(stderr,"Free Energy of Ensemble: %.2f kcal/mol\n", -kT * viterbi.alpha / 100.0);
#else
    if (is_verbose) fprintf(stderr,"Log Partition Coefficient: %.5f\n", viterbi.alpha);
#endif

    if(is_verbose) fprintf(stderr,"Partition Function Calculation Time: %.2f seconds.\n", parse_elapsed_time);

    fflush(stdout);

    // lhuang
    if(pf_only && !forest_file.empty()) dump_forest(seq, true); // inside-only forest

    if(!pf_only){
        outside(next_pair);
    	if (!forest_file.empty())
    	  dump_forest(seq, false); // inside-outside forest
            cal_PairProb(viterbi,RR_list);

        if (mea_) PairProb_MEA(seq);

        if (threshknot_) ThreshKnot(seq);
    }
    postprocess();
    return;
}

void BeamCKYParser::print_states(FILE *fptr, unordered_map<int, State>& states, int j, string label, bool inside_only, double threshold) {    
    for (auto & item : states) {
        int i = item.first;
        State & state = item.second;
        if (inside_only) fprintf(fptr, "%s %d %d %.5lf\n", label.c_str(), i+1, j+1, state.alpha);
        else if (state.alpha + state.beta > threshold) // lhuang : alpha + beta - totalZ < ...
            fprintf(fptr, "%s %d %d %.5lf %.5lf\n", label.c_str(), i+1, j+1, state.alpha, state.beta);
    }
}

void BeamCKYParser::dump_forest(string seq, bool inside_only) {  
    printf("Dumping (%s) Forest to %s...\n", (inside_only ? "Inside-Only" : "Inside-Outside"), forest_file.c_str());
    FILE *fptr = fopen(forest_file.c_str(), "w");  // lhuang: should be fout >>
    fprintf(fptr, "%s\n", seq.c_str());
    int n = seq.length(), j;
    for (j = 0; j < n; j++) {
        if (inside_only) fprintf(fptr, "E %d %.5lf\n", j+1, bestC[j].alpha);
        else fprintf(fptr, "E %d %.5lf %.5lf\n", j+1, bestC[j].alpha, bestC[j].beta);
    }
    double threshold = bestC[n-1].alpha - 9.91152; // lhuang -9.xxx or ?
    for (j = 0; j < n; j++) 
        print_states(fptr, bestP[j], j, "P", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM[j], j, "M", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestM2[j], j, "M2", inside_only, threshold);
    for (j = 0; j < n; j++) 
        print_states(fptr, bestMulti[j], j, "Multi", inside_only, threshold);
}

BeamCKYParser::BeamCKYParser(int beam_size,
                             bool nosharpturn,
                             bool verbose,
                             string bppfile,
                             string bppfileindex,
                             bool pfonly,
                             float bppcutoff,
			                 string forestfile,
                             bool mea,
                             float MEA_gamma,
                             string MEA_file_index,
                             bool MEA_bpseq,
                             bool ThreshKnot,
                             float ThreshKnot_threshold,
                             string ThreshKnot_file_index,
                             string shape_file_path)
    : beam(beam_size), 
      no_sharp_turn(nosharpturn), 
      is_verbose(verbose),
      bpp_file(bppfile),
      bpp_file_index(bppfileindex),
      pf_only(pfonly),
      bpp_cutoff(bppcutoff),
      forest_file(forestfile), 
      mea_(mea),
      gamma(MEA_gamma),
      mea_file_index(MEA_file_index),
      bpseq(MEA_bpseq),
      threshknot_(ThreshKnot),
      threshknot_threshold(ThreshKnot_threshold),
      threshknot_file_index(ThreshKnot_file_index){
#ifdef lpv
        initialize();
#else
        initialize();
        initialize_cachesingle();
#endif

    if (shape_file_path != "" ){
        use_shape = true;
        int position;
        string data;

        double temp_after_mb_shape;

        ifstream in(shape_file_path);

        if (!in.good()){
            cout<<"Reading SHAPE file error!"<<endl;
            assert(false);
        }

        // actually, we can combine the SHAPE_data and the energy_stack together
        while (!(in >> position >> data).fail()) {
            // cout<<"position data "<< int(position)<<endl<<data<<endl;
            // assert(int(position) == SHAPE_data.size() + 1);
            // cout<<"data "<<data<<endl;
            if (isdigit(int(data[0])) == 0){
                SHAPE_data.push_back(double((-1.000000)));
            }

            else {
                SHAPE_data.push_back(stod(data));
            }
            

        }

        for (int i = 0; i<SHAPE_data.size(); i++){
            temp_after_mb_shape = SHAPE_data[i] < 0 ? 0. : (m * log(SHAPE_data[i] + 1) + b);

            pseudo_energy_stack.push_back((int)roundf(temp_after_mb_shape * 100.));

            assert(pseudo_energy_stack.size() == i + 1 );

            // cout<<"pseudo energy "<<i<<' '<<SHAPE_data[i]<<' '<<temp_after_mb_shape<<' '<<pseudo_energy_stack[i]<<' '<<pseudo_energy_stack.size()<<endl;

        }
    }



}

int LinearPartition_main(string seq, vector<pair<float,pair<int,int> > >&RR_list)
{

    struct timeval total_starttime, total_endtime;
    gettimeofday(&total_starttime, NULL);

    int beamsize = 100;
    bool sharpturn = false;
    bool is_verbose = false;
    string bpp_file;
    string bpp_prefix;
    bool pf_only = false;
    float bpp_cutoff = 0.0;
    string forest_file;

    float MEA_gamma = 3.0;
    bool mea = false;
    bool MEA_bpseq = false;
    string MEA_prefix;
    float ThreshKnot_threshold = 0.3;
    bool ThreshKnot = false;
    string ThresKnot_prefix;


    // SHAPE
    string shape_file_path = "";



    //if (argc > 1) {
        //beamsize = atoi(argv[1]);
        //sharpturn = atoi(argv[2]) == 1;
        //is_verbose = atoi(argv[3]) == 1;
        //bp_file = argv[4];
        //bpp_prefix = argv[5];
        //pf_only = atoi(argv[6]) == 1;
        //bpp_cutoff = atof(argv[7]);
        //forest_file = argv[8];
        //mea = atoi(argv[9]) == 1;
        //MEA_gamma = atof(argv[10]);
        //ThreshKnot = atoi(argv[11]) == 1;
        //ThreshKnot_threshold = atof(argv[12]);
        //ThresKnot_prefix = argv[13];
        //MEA_prefix = argv[14];
        //MEA_bpseq = atoi(argv[15]) == 1;
        //shape_file_path = argv[16];
    //}


    if (is_verbose) printf("beam size: %d\n", beamsize);

    // variables for decoding
    int num=0, total_len = 0;
    unsigned long long total_states = 0;
    double total_score = .0;
    double total_time = .0;

    int seq_index = 0;
    string bpp_file_index = "";
    string ThreshKnot_file_index = "";
    string MEA_file_index = "";

    //for (string seq; getline(cin, seq);) 
    {
        if (seq.length() == 0)
            return 1;//continue;

        if (seq[0] == ';' || seq[0] == '>') {
            printf("%s\n", seq.c_str());
            if (!bpp_file.empty()) {
                FILE *fptr = fopen(bpp_file.c_str(), "a"); 
                if (fptr == NULL) { 
                    printf("Could not open file!\n"); 
                    return 0; 
                }
                fprintf(fptr, "%s\n", seq.c_str());
                fclose(fptr); 
            }
            return 1;//continue;
        }

        if (!isalpha(seq[0])){
            printf("Unrecognized sequence: %s\n", seq.c_str());
            return 1;//continue;
        }

        seq_index ++;
        if (!bpp_prefix.empty()) bpp_file_index = bpp_prefix + to_string(seq_index);

        if (!ThresKnot_prefix.empty()) ThreshKnot_file_index = ThresKnot_prefix + to_string(seq_index);

        if (!MEA_prefix.empty()) MEA_file_index = MEA_prefix + to_string(seq_index);
        
        // convert to uppercase
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);

        // convert T to U
        replace(seq.begin(), seq.end(), 'T', 'U');

        // lhuang: moved inside loop, fixing an obscure but crucial bug in initialization
        BeamCKYParser parser(beamsize, !sharpturn, is_verbose, bpp_file, bpp_file_index, pf_only, bpp_cutoff, forest_file, mea, MEA_gamma, MEA_file_index, MEA_bpseq, ThreshKnot, ThreshKnot_threshold, ThreshKnot_file_index, shape_file_path);

        // BeamCKYParser::DecoderResult result = parser.parse(seq);
        parser.parse(seq,RR_list);
    }

    //gettimeofday(&total_endtime, NULL);
    //double total_elapsed_time = total_endtime.tv_sec - total_starttime.tv_sec + (total_endtime.tv_usec-total_starttime.tv_usec)/1000000.0;

    //if(is_verbose) fprintf(stderr,"Total Time: %.2f seconds.\n", total_elapsed_time);

    return 0;
}
