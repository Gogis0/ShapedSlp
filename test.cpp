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

void load_grammar(
    const string filename,
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> &slp)
{
    NaiveSlp<uint32_t> nslp;
    nslp.load_Bigrepair(filename.data(), false);
    nslp.makeBinaryTree();
    slp.init(nslp);
    //slp.printStatus(false);
}

Test(misc, powmod) {
    cr_assert(pow_mod_mersenne(18014398509481984, 2305843009213693949) == 128);
}

Test(misc, tiny_slp) {
    const string filename = "../data/tiny";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);
    cr_assert(slp.getHashWhole() == 2225879177060430003);
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
    const string filename_slp_serialized = filename + ".plainslp.tmp";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename,slp);

    uint64_t hash_whole = slp.getHashWhole();
    ofstream out (filename_slp_serialized);
    slp.serialize(out);
    out.flush();
    out.close();

    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp_new;
    ifstream in (filename_slp_serialized);
    slp_new.load(in);
    //pslp_new.printStatus(true);
    in.close();

    cr_assert(slp_new.getHashWhole() == hash_whole);
}

Test(misc, tiny_slp_mlq) {
    const string filename = "../data/tiny";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);


    string pattern = "ATTCGGATTAGGATTAGG"; // p is actually whole string in the slp
    int m = pattern.length();
    slp.precompute_pattern(pattern);
    cr_assert(match_length_query(slp, 0, 6) == 3);
    cr_assert(match_length_query(slp, 4, 10) == 8);
    cr_assert(match_length_query(slp, 3, 5) == 0);
    cr_assert(match_length_query(slp, 1, 7) == 2);
    cr_assert(match_length_query(slp, 3, 8) == 0);
    cr_assert(match_length_query(slp, 0, 0) == m);
    cr_assert(match_length_query(slp, m-1, m-1) == 1);
}


Test(misc, larger_mlq_special) {
    const string filename = "../data/larger";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);
    

    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    slp.precompute_pattern(pattern);
    
    //cout << pattern[38] << endl;
    cr_assert(match_length_query(slp, 112, 38) == 1);
    cr_assert(match_length_query(slp, 121, 54) == 1);
}

Test(misc, chr19_slp_seq) {
    const string filename = "../data/chr19.1.fa";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);

    std::string pattern("GGGAACTTCTTCT");
    slp.precompute_pattern(pattern);
    cr_assert(match_length_query(slp, 90000, 3) == 9);
    cr_assert(match_length_query(slp, 90001, 4) == 8);
    cr_assert(match_length_query(slp, 90002, 5) == 7);
}
Test(misc, larger_randomized) {
    const string filename = "../data/larger";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);

    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    slp.precompute_pattern(pattern);
    for (int i = 0; i < 1000; i++) {
        const size_t pos1 = rand() % m;
        const size_t pos2 = rand() % m;
        const size_t lce_regular = lceToR(slp, pos1, pos2);
        const size_t lce_mlq = match_length_query(slp, pos1, pos2);
        //std::cout << pos1 << " " << pos2 << " " << lce_regular << " " << lce_mlq << std::endl;
        cr_assert(lce_regular == lce_mlq);
    }
}

Test(misc, chr19_randomized) {
    const string filename = "../data/chr19.1.fa";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);
    
    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    std::string pattern = ss.str();
    const size_t m = pattern.length();
    slp.precompute_pattern(pattern);
    for (int i = 0; i < 100000; i++) {
        const size_t pos1 = rand() % m;
        const size_t pos2 = rand() % m;
        const size_t lce_regular = lceToR(slp, pos1, pos2);
        const size_t lce_mlq = match_length_query(slp, pos1, pos2);
        //std::cout << pos1 << " " << pos2 << " " << lce_regular << " " << lce_mlq << std::endl;
        cr_assert(lce_regular == lce_mlq);
    }
}