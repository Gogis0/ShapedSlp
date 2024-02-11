#ifndef INCLUDE_GUARD_Common
#define INCLUDE_GUARD_Common

#include <cstdint>
#include <stdint.h> // include uint64_t etc.
#include <cassert>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <stack>
#include <vector>

template<typename var_t>
struct PairT
{
  var_t left, right;

  bool operator<(const PairT & another) const
  {
    return (this->left < another.left) || ((this->left == another.left) && this->right < another.right);
  };
};


///////////////////////////////////////////// from Dominik
// source: https://github.com/dominikkempa/lz77-to-slp/blob/main/src/karp_rabin_hashing.cpp

const uint64_t mersenne_prime_exponent = 61;
const uint64_t prime = ((uint64_t)1 << mersenne_prime_exponent) - 1; // 2^61 - 1
uint64_t base = 411910476928516559; // randomly generated

std::uint64_t mod_mersenne(std::uint64_t a) {
  const std::uint64_t k = mersenne_prime_exponent;
  std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  if (k < 32) {

    // We need to check if a <= 2^(2k).
    const std::uint64_t threshold = ((std::uint64_t)1 << (k << 1));
    if (a <= threshold) {
      a = (a & p) + (a >> k);
      a = (a & p) + (a >> k);
      return a == p ? 0 : a;
    } else return a % p;
  } else {

    // We are guaranteed that a < 2^(2k)
    // because a < 2^64 <= 2^(2k).
    a = (a & p) + (a >> k);
    a = (a & p) + (a >> k);
    return a == p ? 0 : a;
  }
}


std::uint64_t mul_mod_mersenne(
    const std::uint64_t a,
    const std::uint64_t b) {
  const uint32_t k = mersenne_prime_exponent;
  const std::uint64_t p = ((std::uint64_t)1 << k) - 1;
  __extension__ const unsigned __int128 ab =
    (unsigned __int128)a *
    (unsigned __int128)b;
  std::uint64_t lo = (std::uint64_t)ab;
  const std::uint64_t hi = (ab >> 64);
  lo = (lo & p) + ((lo >> k) + (hi << (64 - k)));
  lo = (lo & p) + (lo >> k);
  return lo == p ? 0 : lo;
}

// takes O(log n) time
std::uint64_t pow_mod_mersenne(
    const std::uint64_t a,
    std::uint64_t n) {
  const uint32_t k = mersenne_prime_exponent;
  std::uint64_t pow = mod_mersenne(a);
  std::uint64_t ret = mod_mersenne(1);
  while (n > 0) {
    //std::cout << n << std::endl;
    if (n & 1)
      ret = mul_mod_mersenne(ret, pow);
    pow = mul_mod_mersenne(pow, pow);
    n >>= 1;
  }
  return ret;
}

std::uint64_t concat(
    const std::uint64_t left_hash,
    const std::uint64_t right_hash,
    const std::uint64_t right_len,
    const std::uint64_t base) {
  return mod_mersenne(mul_mod_mersenne(pow_mod_mersenne(base, right_len), left_hash)
                      + right_hash);
}

std::uint64_t hash_char(const char c) {
  return mod_mersenne((std::uint64_t)c);
}

//////////////////////////////////////////// until here

void precompute_pattern_hashes(const std::string& pattern, std::vector<uint64_t>& hashes, std::vector<uint64_t>& powers, uint64_t base) {
    int m = pattern.size();
    hashes.resize(m);
    powers.resize(m);

    uint64_t hash = 0;
    uint64_t power = 1;

    // precompute: powers[i] = base^i
    for (std::uint64_t i = 0; i < m; i++) {
        powers[i] = power;
        power = mul_mod_mersenne(power, base);
    }

    hashes[0] = hash_char(pattern[0]);
    // precompute: hashes[i] = P[0] * x^(m-1) + P[1] * x^(m-2) + ... + P[i-1]
    for (std::uint64_t i = 1; i < m; i++) {
        hashes[i] = mod_mersenne(mul_mod_mersenne(hashes[i-1], base) + hash_char(pattern[i]));
    }
}

uint64_t compute_subpattern_hash(const std::vector<uint64_t> &prefix_hashes,
                                 const std::vector<uint64_t> &base_powers,
                                 const int i,
                                 const int j) {
    if (i == 0) {
        return prefix_hashes[j];
    }
    
    const uint64_t pref_i = prefix_hashes[i - 1];
    const uint64_t pref_j = prefix_hashes[j];
    const uint64_t base_pow = base_powers[j - i + 1];
    
    long long int hash_diff = pref_j;
    if (pref_i > 0) {
        hash_diff = (hash_diff - mul_mod_mersenne(pref_i, base_pow) + prime) % prime;
    }
    
    return hash_diff;
}

void padVLine
(
 const uint64_t pad
 ) {
  for (uint64_t i = 0; i < pad; ++i) {
    std::cout << "|";
  }
}


uint32_t ceilLog2(uint64_t x) {
  if (x == 0) {
    return 1;
  }
  return 64 - __builtin_clzll(x);
}


template<class type>
void printArray(type arr, uint64_t n, std::string delimiter = " ")
{
  for (uint64_t i = 0; i < n; ++i) {
    std::cout << arr[i] << delimiter;
  }
  std::cout << std::endl;
}


template<class Container>
void printVec(const Container & vec)
{
  for (uint64_t i = 0; i < vec.size(); ++i) {
    std::cout << "(" << i << ":" << vec[i] << ") ";
  }
  std::cout << std::endl;
}


template<class SlpT>
void decompressByCharAt
(
 const SlpT & slp,
 const std::string & ofile
) {
  std::ofstream ofs(ofile);
  for (uint64_t i = 0; i < slp.getLen(); ++i) {
    char c = slp.charAt(i);
    ofs.write(&c, 1);
  }
}


template<class SlpT>
void decompressByRecurse
(
 const SlpT & slp,
 const std::string & ofile
) {
  using var_t = typename SlpT::var_t;
  const uint64_t alphSize = slp.getAlphSize();
  std::ofstream ofs(ofile);
  for (uint64_t i = 0; i < slp.getLenSeq(); ++i) {
    const var_t var = slp.getSeq(i);
    std::stack<std::pair<var_t, uint64_t> > st;
    st.push(std::make_pair(var, 0));
    while (!st.empty()) {
      auto & e = st.top();
      if (e.first < alphSize) { // leaf
        char c = slp.getChar(e.first);
        ofs.write(&c, 1);
        st.pop();
      } else if (e.second == 0) {
        e.second = 1;
        st.push(std::make_pair(slp.getLeft(e.first - alphSize), 0));
      } else if (e.second == 1) {
        e.second = 2;
        st.push(std::make_pair(slp.getRight(e.first - alphSize), 0));
      } else {
        st.pop();
      }
    }
  }
}


template<class SlpT>
void printDerivationTree
(
 const SlpT & slp
) {
  std::stack<typename SlpT::nodeT> path;
  path.push(slp.getRootNode());
  if (path.size() == 0) {
    return;
  }
  while (true) {
    auto node = path.top();
    padVLine(path.size() - 1);
    std::cout << " (" << std::get<0>(node) << ", " << std::get<1>(node) << ", " << std::get<2>(node) << ")";
    if (std::get<0>(node) == 1) {
      std::cout << " " << (char)(std::get<1>(node));
    }
    std::cout << std::endl;
    if (std::get<0>(node) == 1) {
      if (!proceedPrefixPath(slp, path)) {
        return;
      }
    } else {
      descentPrefixPath(slp, path, std::get<0>(node) - 1);
    }
  }
}



template<class SlpT>
void getPrefixPath
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 uint64_t pos
) {
  if (pos >= slp.getLen()) {
    return;
  }
  path.push(slp.getRootNode());
  if (pos) {
    path.push(slp.getChildNodeForPos_Root(pos));
  }
  while (pos) {
    path.push(slp.getChildNodeForPos(path.top(), pos)); // pos is modified to relative pos in a node
  }
}

template<class SlpT>
void completePrefixPath
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 std::stack<uint64_t> & positions
 ) {
  auto n = path.top();
  while (std::get<0>(path.top()) > 1) {
    std::cout << "n.size = " << std::get<0>(n) << std::endl;
    n = slp.getChildNode(path.top(), 0);
    path.push(n);
    positions.push(0);
  }
}

/*!
 * modify the stack 'path' to point the highest node that is adjacent to the node path.top()
 * return false when such a node does not exist
 */
template<class SlpT>
bool proceedPrefixPath
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path
 ) {
  if (path.size() <= 1) {
    return false;
  }
  typename SlpT::nodeT n;
  do {
    n = path.top();
    path.pop();
  } while (path.size() > 1 and std::get<2>(n) == 1);
  if (path.size() > 1) {
    path.push(slp.getChildNode(path.top(), 1));
  } else { // add (std::get<2>(n) + 1)th (0base) child of root
    if (std::get<2>(n) + 1 < slp.getLenSeq()) {
      path.push(slp.getChildNode_Root(std::get<2>(n) + 1));
    } else {
      return false;
    }
  }
  return true;
}


template<class SlpT>
void descentPrefixPath
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 const uint64_t len
 ) {
  auto n = (path.size() == 1) ? slp.getChildNode_Root(0) : slp.getChildNode(path.top(), 0);
  path.push(n);
  while (std::get<0>(n) > len) {
    n = slp.getChildNode(path.top(), 0);
    path.push(n);
  }
}


template<class SlpT>
uint64_t lceToR
(
 const SlpT & slp,
 const uint64_t p1,
 const uint64_t p2
) {
  std::stack<typename SlpT::nodeT> path1, path2;

  getPrefixPath(slp, path1, p1);
  getPrefixPath(slp, path2, p2);

  uint64_t l = 0;
  while (true) {
    auto n1 = path1.top();
    auto n2 = path2.top();
    while (std::get<0>(n1) != std::get<0>(n2)) {
      if (std::get<0>(n1) > std::get<0>(n2)) {
        descentPrefixPath(slp, path1, std::get<0>(n2));
        n1 = path1.top();
        // std::cout << "descent n1: " << std::get<0>(n1) << ", " << std::get<1>(n1) << ", " << std::get<2>(n1) << std::endl;
      } else {
        descentPrefixPath(slp, path2, std::get<0>(n1));
        n2 = path2.top();
        // std::cout << "descent n2: " << std::get<0>(n2) << ", " << std::get<1>(n2) << ", " << std::get<2>(n2) << std::endl;
      }
    }
    if (std::get<1>(n1) == std::get<1>(n2)) { // match
      l += std::get<0>(n1);
      if (!(proceedPrefixPath(slp, path1))) {
        break;
      }
      if (!(proceedPrefixPath(slp, path2))) {
        break;
      }
    } else if (std::get<0>(n1) > 1) { // mismatch with non-terminal
      descentPrefixPath(slp, path1, std::get<0>(n1) - 1);
      descentPrefixPath(slp, path2, std::get<0>(n1) - 1);
    } else { // lce ends with mismatch char
      break;
    }
  }
  return l;
}

// returns the path to the node X such that exp(X) prefixes T[pos..n]
template<class SlpT>
void getPrefixPathWithPositions
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 uint64_t pos
) {
  if (pos >= slp.getLen()) {
    return;
  }
  auto n = slp.getChildNodeForPos_Root(pos);

  while (std::get<0>(n) > 1) {
    auto x = slp.getChildNodeForPos(n, pos);
    if (std::get<2>(x) == 0) {
      path.push(n);
      //std::cout << "node len: " << std::get<0>(n) << std::endl;
    }
    n = x;
  }
  path.push(n);
}

template<class SlpT>
void getSuffixPathWithPositions
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 uint64_t pos
) {
  if (pos >= slp.getLen()) {
    return;
  }
  auto n = slp.getChildNodeForPos_Root(pos);

  while (std::get<0>(n) > 1) {
    auto x = slp.getChildNodeForPos(n, pos);
    if (std::get<2>(x) == 1) {
      path.push(n);
      //std::cout << "node len: " << std::get<0>(n) << std::endl;
    }
    n = x;
  }
  path.push(n);
}

// compute LCP(i, j) w.h.p. in O(grammar_height) time
template<class SlpT>
uint64_t LCP
(
    const SlpT & slp,
    const uint64_t pos_s, // position in the string
    const uint64_t pos_p // position in the pattern
)
{
    std::stack<typename SlpT::nodeT> path;
    getPrefixPathWithPositions(slp, path, pos_s);

    auto x = path.top();
    path.pop();

    const size_t m = slp.getPatternLen();
    const size_t alpha = slp.getAlphSize();
    uint64_t node = std::get<1>(x);
    uint64_t last_from_left = node;
    uint64_t pattern_prefix_hash = slp.subpattern_hash(pos_p, pos_p); // hash of the first prefix character
    uint64_t subpattern_hash = 0; // hash of the first T prefix character
    uint64_t pattern_prefix_len = 1;

    if (node != pattern_prefix_hash) return 0;

    // ASCENT
    while (!path.empty()) {
        last_from_left = node;
        x = path.top();
        path.pop();
        node = std::get<1>(x);

        const uint64_t right_child = slp.getRight(node - alpha);
        const uint64_t right_child_len = slp.getNodeLen(right_child);
        const uint64_t new_pattern_prefix_len = pattern_prefix_len + right_child_len;

        if ((pos_p + new_pattern_prefix_len) >= m) break;
        const uint64_t right_hash = slp.getNodeHash(right_child);
        subpattern_hash = slp.subpattern_hash(pos_p+pattern_prefix_len, pos_p+new_pattern_prefix_len-1);
        if (right_hash != subpattern_hash) break;
        pattern_prefix_len += right_child_len;
    }

    if (node < alpha) return pattern_prefix_len;
    node = slp.getRight(node - alpha);

    // RE-DESCENT
    while (true) {
      if (node < alpha) {
        const uint64_t last_char_hash = slp.subpattern_hash(pos_p+pattern_prefix_len, pos_p+pattern_prefix_len);
        pattern_prefix_len += (node == last_char_hash);
        break;
      }
      const uint64_t left_child = slp.getLeft(node - alpha);
      const uint64_t right_child = slp.getRight(node - alpha);
      const uint64_t left_size = slp.getNodeLen(left_child);
      const uint64_t left_hash = slp.getNodeHash(left_child);
      const uint64_t new_pattern_prefix_len = pattern_prefix_len+left_size;
      if (((pos_p + new_pattern_prefix_len - 1) >= m) ||
          (left_hash != slp.subpattern_hash(pos_p+pattern_prefix_len, pos_p+new_pattern_prefix_len-1))) {
          node = left_child;
      } else {
        node = right_child;
        pattern_prefix_len += left_size;
      }
    }
    return pattern_prefix_len;
}

// compute LCS(i, j) w.h.p. in O(grammar_height) time
template<class SlpT>
uint64_t LCS
(
    const SlpT & slp,
    const uint64_t pos_s, // position in the string
    const uint64_t pos_p // position in the pattern
)
{
    std::stack<typename SlpT::nodeT> path;
    getSuffixPathWithPositions(slp, path, pos_s);

    auto x = path.top();
    path.pop();

    const size_t alpha = slp.getAlphSize();
    uint64_t node = std::get<1>(x);
    uint64_t last_from_right = node;
    uint64_t pattern_suffix_hash = slp.subpattern_hash(pos_p, pos_p); // hash of the first prefix character
    uint64_t subpattern_hash = 0; // hash of the first T prefix character
    uint64_t pattern_suffix_len = 1;

    if (node != pattern_suffix_hash) return 0;

    // ASCENT
    while (!path.empty()) {
        last_from_right = node;
        x = path.top();
        path.pop();
        node = std::get<1>(x);

        const uint64_t left_child = slp.getLeft(node - alpha);
        const uint64_t left_child_len = slp.getNodeLen(left_child);
        const uint64_t new_pattern_suffix_len = pattern_suffix_len + left_child_len;

        if ((long long)(pos_p - new_pattern_suffix_len+1) < 0) break;
        const uint64_t left_hash = slp.getNodeHash(left_child);
        subpattern_hash = slp.subpattern_hash(pos_p-new_pattern_suffix_len+1, pos_p-pattern_suffix_len);
        if (left_hash != subpattern_hash) break;
        pattern_suffix_len += left_child_len;
    }

    if (node < alpha) return pattern_suffix_len;
    node = slp.getLeft(node - alpha);

    // RE-DESCENT
    while (true) {
      if (node < alpha) {
        const uint64_t last_char_hash = slp.subpattern_hash(pos_p-pattern_suffix_len, pos_p-pattern_suffix_len);
        pattern_suffix_len += (node == last_char_hash);
        break;
      }
      const uint64_t left_child = slp.getLeft(node - alpha);
      const uint64_t right_child = slp.getRight(node - alpha);
      const uint64_t right_size = slp.getNodeLen(right_child);
      const uint64_t right_hash = slp.getNodeHash(right_child);
      const uint64_t new_pattern_prefix_len = pattern_suffix_len + right_size;
      if (((long long)(pos_p - new_pattern_prefix_len+1) < 0) ||
          (right_hash != slp.subpattern_hash(pos_p-new_pattern_prefix_len+1, pos_p-pattern_suffix_len))) {
          node = right_child;
      } else {
        node = left_child;
        pattern_suffix_len += right_size;
      }
    }
    return pattern_suffix_len;
}

// compute LCP(i, j) in O(LCP(i, j) * grammar_height) time
template<class SlpT>
uint64_t naive_LCP(
    const std::string &pattern,
    SlpT &slp,
    size_t i,
    size_t j)
{
    const size_t n = slp.getLen();
    const size_t m = pattern.length();
    uint64_t ans = 0;
    while ((i < n) && (j < m) && (slp.charAt(i) == pattern[j])) {
        i++;
        j++;
        ans++;
    }
    return ans; 
}

// compute LCS(i, j) in O(LCS(i, j) * grammar_height) time
template<class SlpT>
uint64_t naive_LCS(
    const std::string &pattern,
    SlpT &slp,
    size_t i,
    size_t j)
{
    uint64_t ans = 0;
    while ((i >= 0) && (j >= 0) && (slp.charAt(i) == pattern[j])) {
        i--;
        j--;
        ans++;
    }
    return ans; 
}

template<class SlpT>
uint64_t lceToRBounded
(
 const SlpT & slp,
 const uint64_t p1,
 const uint64_t p2,
 const uint64_t upperbound
) {
  std::stack<typename SlpT::nodeT> path1, path2;

  getPrefixPath(slp, path1, p1);
  getPrefixPath(slp, path2, p2);

  uint64_t l = 0;
  while (true) {
    auto n1 = path1.top();
    auto n2 = path2.top();
    while (std::get<0>(n1) != std::get<0>(n2)) {
      if (std::get<0>(n1) > std::get<0>(n2)) {
        descentPrefixPath(slp, path1, std::get<0>(n2));
        n1 = path1.top();
        // std::cout << "descent n1: " << std::get<0>(n1) << ", " << std::get<1>(n1) << ", " << std::get<2>(n1) << std::endl;
      } else {
        descentPrefixPath(slp, path2, std::get<0>(n1));
        n2 = path2.top();
        // std::cout << "descent n2: " << std::get<0>(n2) << ", " << std::get<1>(n2) << ", " << std::get<2>(n2) << std::endl;
      }
    }
    if (std::get<1>(n1) == std::get<1>(n2)) { // match
      l += std::get<0>(n1);
      if(l >= upperbound) { return l; }
      if (!(proceedPrefixPath(slp, path1))) {
        break;
      }
      if (!(proceedPrefixPath(slp, path2))) {
        break;
      }
    } else if (std::get<0>(n1) > 1) { // mismatch with non-terminal
      descentPrefixPath(slp, path1, std::get<0>(n1) - 1);
      descentPrefixPath(slp, path2, std::get<0>(n1) - 1);
    } else { // lce ends with mismatch char
      break;
    }
  }
  return l;
}


template<class SlpT>
uint64_t lceToR_Naive
(
 const SlpT & slp,
 const uint64_t p1,
 const uint64_t p2
) {
  std::stack<typename SlpT::nodeT> path1, path2;

  getPrefixPath(slp, path1, p1);
  getPrefixPath(slp, path2, p2);

  uint64_t l = 0;
  while (true) {
    auto n1 = path1.top();
    auto n2 = path2.top();
    if (std::get<0>(n1) != 1) {
      descentPrefixPath(slp, path1, 1);
      n1 = path1.top();
    }
    if (std::get<0>(n2) != 1) {
      descentPrefixPath(slp, path2, 1);
      n2 = path2.top();
    }
    if (std::get<1>(n1) == std::get<1>(n2)) { // match char
      ++l;
      if (!(proceedPrefixPath(slp, path1))) {
        break;
      }
      if (!(proceedPrefixPath(slp, path2))) {
        break;
      }
    } else { // lce ends with mismatch char
      break;
    }
  }
  return l;
}

template<class SlpT>
uint64_t lceToR_NaiveBounded
(
 const SlpT & slp,
 const uint64_t p1,
 const uint64_t p2,
 const uint64_t upperbound
) {
  std::stack<typename SlpT::nodeT> path1, path2;

  getPrefixPath(slp, path1, p1);
  getPrefixPath(slp, path2, p2);

  uint64_t l = 0;
  while (true) {
    auto n1 = path1.top();
    auto n2 = path2.top();
    if (std::get<0>(n1) != 1) {
      descentPrefixPath(slp, path1, 1);
      n1 = path1.top();
    }
    if (std::get<0>(n2) != 1) {
      descentPrefixPath(slp, path2, 1);
      n2 = path2.top();
    }
    if (std::get<1>(n1) == std::get<1>(n2)) { // match char
      ++l;
      if(l >= upperbound) { return l; }
      if (!(proceedPrefixPath(slp, path1))) {
        break;
      }
      if (!(proceedPrefixPath(slp, path2))) {
        break;
      }
    } else { // lce ends with mismatch char
      break;
    }
  }
  return l;
}


// template<class SlpT>
// uint64_t lceToR
// (
//  const SlpT & slp,
//  const uint64_t p1,
//  const uint64_t p2
// ) {
//   std::stack<typename SlpT::nodeT> s1, s2;
//   slp.init_LceToR(s1, p1);
//   slp.init_LceToR(s2, p2);

//   // {
//   //   std::stack<typename SlpT::nodeT> s;
//   //   slp.init_LceToR(s, p1);
//   //   std::cout << "p1 = " << p1 << std::endl;
//   //   while (!s.empty()) {
//   //     std::cout << std::get<0>(s.top()) << ", " << std::get<1>(s.top()) << ", " << std::get<2>(s.top()) << std::endl;
//   //     s.pop();
//   //   }
//   // }
//   // {
//   //   std::stack<typename SlpT::nodeT> s;
//   //   slp.init_LceToR(s, p2);
//   //   std::cout << "p2 = " << p2 << std::endl;
//   //   while (!s.empty()) {
//   //     std::cout << std::get<0>(s.top()) << ", " << std::get<1>(s.top()) << ", " << std::get<2>(s.top()) << std::endl;
//   //     s.pop();
//   //   }
//   // }

//   uint64_t l = 0;
//   while (true) {
//     auto n1 = s1.top();
//     auto n2 = s2.top();
//     while (std::get<0>(n1) != std::get<0>(n2)) {
//       if (std::get<0>(n1) > std::get<0>(n2)) {
//         slp.descent_LceToR(s1, std::get<0>(n2));
//         n1 = s1.top();
//         // std::cout << "descent n1: " << std::get<0>(n1) << ", " << std::get<1>(n1) << ", " << std::get<2>(n1) << std::endl;
//       } else {
//         slp.descent_LceToR(s2, std::get<0>(n1));
//         n2 = s2.top();
//         // std::cout << "descent n2: " << std::get<0>(n2) << ", " << std::get<1>(n2) << ", " << std::get<2>(n2) << std::endl;
//       }
//     }
//     if (std::get<1>(n1) == std::get<1>(n2)) { // match
//       l += std::get<0>(n1);
//       if (!(slp.next_LceToR(s1))) {
//         break;
//       }
//       if (!(slp.next_LceToR(s2))) {
//         break;
//       }
//     } else if (std::get<0>(n1) > 1) {
//       slp.descent_LceToR(s1, std::get<0>(n1) - 1);
//       slp.descent_LceToR(s2, std::get<0>(n1) - 1);
//     } else { // lce ends with mismatch char
//       break;
//     }
//   }
//   return l;
// }


template<typename ArrayT>
class PackedArrayTypeValRef
{
  friend ArrayT;


private:
  ArrayT * const obj_;
  const uint64_t idx_;


  PackedArrayTypeValRef
  (
   ArrayT * obj,
   uint64_t idx
   ) :
    obj_(obj),
    idx_(idx)
  {}


public:
  uint64_t operator=
  (
   uint64_t val
   ) {
    obj_->write(val, idx_);
    return val;
  }


  operator uint64_t() const {
    return obj_->read(idx_);
  }
};




/*!
 * @tparam kBucketWidth: Bitwidth of bucket size
 * @tparam elem_t: type of element to be sorted
 * @tparam
 *   Func: Fuction that returns kBucketWidth width integer from an element
 */
template<uint8_t kBucketWidth = 8, class elem_t, class Func>
void my_bucket_sort
(
 elem_t * earray, //!< given array to be sorted by some criterion specified by func
 uint64_t n, //!< length of array
 Func func
 ) {
  const uint64_t kBS = UINT64_C(1) << kBucketWidth; // bucket size
  std::queue<elem_t> bucket[kBS];
  for (uint64_t i = 0; i < n; ++i) {
    auto e = earray[i];
    bucket[func(e)].push(e);
  }
  uint64_t i = 0;
  for (uint64_t k = 0; k < kBS; ++k) {
    while (!bucket[k].empty()) {
      elem_t e = bucket[k].front();
      bucket[k].pop();
      earray[i++] = e;
    }
  }
}


/*!
 * @tparam kBucketWidth: Bitwidth of bucket size
 * @tparam elem_t: type of element to be sorted
 * @tparam
 *   Func: Fuction that returns kBucketWidth width integer from an element
 */
template<uint8_t kBucketWidth = 8, class elem_t, class keys_t>
void my_radix_sort
(
 elem_t * earray, //!< given array to be sorted by some criterion specified by func
 keys_t * keys, //!< i \in [0..n) is sorted based on keys[i]
 uint64_t n, //!< length of array
 uint8_t keyWidth
 ) {
  const uint64_t mask = (UINT64_C(1) << kBucketWidth) - 1; // bucket size
  for (uint64_t k = 0; k < (keyWidth + kBucketWidth - 1) / kBucketWidth; ++k) {
    my_bucket_sort<8>
      (earray, n, 
       [keys, k](uint64_t i){
         return (keys[i] >> (kBucketWidth * k)) & mask;
       }
       );
  }
}



#endif
