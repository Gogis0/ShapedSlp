#include <criterion/assert.h>
#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <cstdlib> 
#include "Common.hpp"
#include "PlainSlp.hpp"
#include "FixedBitLenCode.hpp"
using namespace std;

Test(misc, powmod) {
    cr_assert(pow_mod_mersenne(18014398509481984, 2305843009213693949) == 128);
}

Test(misc, tiny_slp) {
    const string in = "../data/tiny";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(in.data(), false);
    slp.makeBinaryTree();
    //slp.printStatus(false);
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);
    cr_assert(pslp.getHashWhole() == 2225879177060430003);
}

Test(misc, pattern_hash_small) {
    string p = "GGGAACTTCTTCT";
    vector<uint64_t> hashes, powers;
    const int base = 256;
    precompute_pattern_hashes(p, hashes, powers, base);
    cr_assert(hashes[p.size()-1] == 1464872621604202849);
    cr_assert(compute_subpattern_hash(hashes, powers, 0, p.size()-1) == 1464872621604202849);
    cr_assert(compute_subpattern_hash(hashes, powers, 3, 3) == 65);
    cr_assert(compute_subpattern_hash(hashes, powers, 3, 5) == 4407617);
}

Test(misc, pattern_hash_larger) {
    const string filename = "../data/larger";
    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    
    vector<uint64_t> hashes, powers;
    const int base = 256;
    precompute_pattern_hashes(pattern, hashes, powers, base);
    for (int i = 0; i < m; i++) {
        cr_assert(compute_subpattern_hash(hashes, powers, i, i) ==  hash_char(pattern[i]));
    }
}


Test(misc, save_load_bigrepair) {
    const string filename = "../data/tiny";
    const string filename_serialized = filename + ".plainslp.tmp";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(filename.data(), false);
    slp.makeBinaryTree();
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);
    uint64_t hash_whole = pslp.getHashWhole();
    ofstream out (filename_serialized);
    pslp.serialize(out);
    out.flush();
    out.close();

    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp_new;
    ifstream in (filename_serialized);
    pslp_new.load(in);
    //pslp_new.printStatus(true);
    in.close();

    cr_assert(pslp_new.getHashWhole() == hash_whole);
}

Test(misc, tiny_slp_mlq) {
    const string in = "../data/tiny";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(in.data(), false);
    slp.makeBinaryTree();
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);


    string pattern = "ATTCGGATTAGGATTAGG"; // p is actually whole string in the slp
    int m = pattern.length();
    pslp.precompute_pattern(pattern);
    cr_assert(match_length_query(pslp, 0, 6) == 3);
    cr_assert(match_length_query(pslp, 4, 10) == 8);
    cr_assert(match_length_query(pslp, 3, 5) == 0);
    cr_assert(match_length_query(pslp, 1, 7) == 2);
    cr_assert(match_length_query(pslp, 3, 8) == 0);
    cr_assert(match_length_query(pslp, 0, 0) == m);
    cr_assert(match_length_query(pslp, m-1, m-1) == 1);
}


Test(misc, larger_slp_mlq) {
    const string filename = "../data/larger";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(filename.data(), false);
    slp.makeBinaryTree();
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);


    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    pslp.precompute_pattern(pattern);
    
    //cout << pattern[38] << endl;
    cr_assert(match_length_query(pslp, 112, 38) == 1);
    cr_assert(match_length_query(pslp, 121, 54) == 1);
}

Test(misc, chr19_slp_seq) {
    const string in = "../data/chr19.1.fa";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(in.data(), false);
    slp.makeBinaryTree();
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);

    std::string pattern("GGGAACTTCTTCT");
    pslp.precompute_pattern(pattern);
    cr_assert(match_length_query(pslp, 90000, 3) == 9);
    cr_assert(match_length_query(pslp, 90001, 4) == 8);
    cr_assert(match_length_query(pslp, 90002, 5) == 7);
}
Test(misc, larger_randomized) {
    const string filename = "../data/larger";
    NaiveSlp<uint32_t> slp;
    slp.load_Bigrepair(filename.data(), false);
    slp.makeBinaryTree();
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> pslp;
    pslp.init(slp);

    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    pslp.precompute_pattern(pattern);
    for (int i = 0; i < 100; i++) {
        const size_t pos1 = rand() % m;
        const size_t pos2 = rand() % m;
        const size_t lce_regular = lceToR(pslp, pos1, pos2);
        const size_t lce_mlq = match_length_query(pslp, pos1, pos2);
        //std::cout << pos1 << " " << pos2 << " " << lce_regular << " " << lce_mlq << std::endl;
        cr_assert(lce_regular == lce_mlq);
    }
}