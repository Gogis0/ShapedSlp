#include <criterion/assert.h>
#include <criterion/criterion.h>
#include <criterion/internal/assert.h>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <queue>
#include <cstdlib> 
#include "Common.hpp"
#include "PlainSlp.hpp"
#include "FixedBitLenCode.hpp"
using namespace std;

string load_pattern(const string filename) {
    ifstream in(filename);
    stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

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
    cr_assert(pow_mod_mersenne(18014398509481984, 0) == 1);
    cr_assert(pow_mod_mersenne(18014398509481984, 1) == 18014398509481984);
    cr_assert(pow_mod_mersenne(18014398509481984, 123) == 18014398509481984);
    cr_assert(pow_mod_mersenne(18014398509481984, 123456) == 72057594037927936);
    cr_assert(pow_mod_mersenne(18014398509481984, 2305843009213693949) == 128);
}

Test(misc, load_grammar) {
    const string filename = "../data/tiny";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);
    cr_assert(slp.getLen() == 18);
    cr_assert(slp.getHashWhole() == 1885666546310520423);
}

Test(misc, pattern_hash_small) {
    const string p = "GGGAACTTCTTCT";
    const size_t m = p.length();
    vector<uint64_t> hashes, powers;
    const int base = 256;
    precompute_pattern_hashes(p, hashes, powers, base);
    cr_assert(hashes[m-1] == 1464872621604202849);
    cr_assert(compute_subpattern_hash(hashes, powers, 2, 3) == 16711);
    cr_assert(compute_subpattern_hash(hashes, powers, 3, 5) == 4407617);
    cr_assert(compute_subpattern_hash(hashes, powers, m-2, m-1) == 21571);
    cr_assert(compute_subpattern_hash(hashes, powers, 0, m-1) == 1464872621604202849);
}

Test(misc, pattern_hash_single_chars) {
    const std::string pattern = load_pattern("../data/pattern1");
    const size_t m = pattern.length();
    
    vector<uint64_t> hashes, powers;
    const int base = 256;
    precompute_pattern_hashes(pattern, hashes, powers, base);
    for (int i = 0; i < m; i++) {
        cr_assert(compute_subpattern_hash(hashes, powers, i, i) == hash_char(pattern[i]));
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


Test(misc, mlq_special_cases) {
    const string filename = "../data/pattern1";
    const std::string pattern = load_pattern(filename);
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);

    slp.precompute_pattern(pattern);
    
    //cout << pattern[38] << endl;
    cr_assert(match_length_query(slp, 112, 38) == 1);
    cr_assert(match_length_query(slp, 121, 54) == 1);
}

Test(misc, chr19_16_mlq_agains_naive) {
    std::string pattern = load_pattern("../data/pattern2");

    const string grammar_file = "../data/chr19.16.fa.plain.slp";
    ifstream in(grammar_file);
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    slp.load(in);
    
    const size_t n = slp.getLen();
    const size_t m = pattern.length();
    slp.precompute_pattern(pattern);
    
    //cases that were causing touble before
    cr_assert(match_length_query(slp, 699280514, 92) == 0);
    cr_assert(match_length_query(slp, 102524317, 136) == naive_MLQ(pattern, slp, 102524317, 136));

    for (int i = 0; i < 100000; i++) {
        const size_t pos1 = rand() % n;
        const size_t pos2 = rand() % m;
        const uint64_t lce_naive = naive_MLQ(pattern, slp, pos1, pos2);
        const size_t lce_mlq = match_length_query(slp, pos1, pos2);
        //std::cout << pos1 << " " << pos2 << " naive: " << lce_naive << " mlq: " << lce_mlq << std::endl;
        cr_assert(lce_mlq == lce_naive);
    }

}

Test(misc, chr19slp_small_pattern) {
    const string filename = "../data/chr19.1.fa";
    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);

    std::string pattern("GGGAACTTCTTCT");
    slp.precompute_pattern(pattern);
    cr_assert(match_length_query(slp, 90000, 3) == 9);
    cr_assert(match_length_query(slp, 90001, 4) == 8);
    cr_assert(match_length_query(slp, 90002, 5) == 7);
}

Test(misc, randomized_pos_small) {
    const string filename = "../data/pattern1";
    std::string pattern = load_pattern(filename);
    const size_t m = pattern.length();

    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);


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

Test(misc, randomized_pos_chr19) {
    const string filename = "../data/chr19.1.fa";

    std::string pattern = load_pattern(filename);
    const size_t m = pattern.length();

    PlainSlp<uint32_t, FixedBitLenCode<>, FixedBitLenCode<>> slp;
    load_grammar(filename, slp);
    
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