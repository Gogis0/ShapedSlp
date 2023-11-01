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

template<typename var_t>
struct PairT
{
  var_t left, right;

  bool operator<(const PairT & another) const
  {
    return (this->left < another.left) || ((this->left == another.left) && this->right < another.right);
  };
};

const uint64_t mersenne_prime_exponent = 61;
const uint64_t prime = ((uint64_t)1 << mersenne_prime_exponent) - 1; // 2^61 - 1

///////////////////////////////////////////// from Dominik

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

std::uint64_t pow_mod_mersenne(
    const std::uint64_t a,
    std::uint64_t n) {
  const uint32_t k = mersenne_prime_exponent;
  std::uint64_t pow = mod_mersenne(a);
  std::uint64_t ret = mod_mersenne(1);
  while (n > 0) {
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
  return mod_mersenne(left_hash + mul_mod_mersenne(pow_mod_mersenne(base, right_len), right_hash));
  const std::uint64_t pow = pow_mod_mersenne(
      base, right_len);
  const std::uint64_t tmp = mul_mod_mersenne(
      left_hash, pow);
  const std::uint64_t ret = mod_mersenne(
      tmp + right_hash);
  return ret;
}

//////////////////////////////////////////// until here

std::uint64_t hash_char(const char c) {
  return mod_mersenne((std::uint64_t)c);
}

void precompute_pattern_hashes(
          const std::string &pattern,
          std::vector<uint64_t> &hashes,
          std::vector<uint64_t> &powers,
          uint64_t base) {
    int m = pattern.size();
    hashes.resize(m);
    powers.resize(m);

    uint64_t h = 0, b = 1;
    for (std::uint64_t i = 0; i < m; i++) {
      h = mod_mersenne(h + mul_mod_mersenne(b, pattern[i]));
      hashes[i] = h;
      powers[i] = b;
      b = mul_mod_mersenne(b, base);
    }
    //std::cout << "final pattern hash = " << hashes[m-1] << std::endl;
}

uint64_t compute_subpattern_hash(const std::vector<uint64_t> &prefix_hashes,
                                 const std::vector<uint64_t> &base_powers,
                                 const int i,
                                 const int j) {
    if (i == 0) return prefix_hashes[j];
    const uint64_t pref_i = prefix_hashes[i - 1];
    uint64_t pref_j = prefix_hashes[j];
    const uint64_t base_pow_inv = pow_mod_mersenne(base_powers[i], prime-2);
    const uint64_t hash_diff = ((long long)(pref_j - pref_i) + ((pref_i > pref_j) ? prime : 0)) % prime;
    const uint64_t ans = mul_mod_mersenne(hash_diff, base_pow_inv);
    return ans;
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

template<class SlpT>
void getPrefixPathWithPositions
(
 const SlpT & slp,
 std::stack<typename SlpT::nodeT> & path,
 std::stack<uint64_t> & positions,
 uint64_t pos
) {
  if (pos >= slp.getLen()) {
    return;
  }
  //path.push(slp.getRootNode());
  //positions.push(pos);
  //if (pos) {
  //path.push(slp.getRootNode());
  //positions.push(pos);
  path.push(slp.getChildNodeForPos_Root(pos));
  positions.push(pos);
  //}
  //while (pos) {
  //  path.push(slp.getChildNodeForPos(path.top(), pos)); // pos is modified to relative pos in a node
  //  std::cout << "node id: " << std::get<1>(path.top()) << std::endl;
  //  positions.push(pos);
  //}
  auto n = path.top();
  while (std::get<0>(path.top()) > 1) {
    n = slp.getChildNodeForPos(path.top(), pos);
    path.push(n);
    positions.push(pos);
  }
  assert(std::get<0>(path.top()) == pos);
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

template<class SlpT>
uint64_t match_length_query
(
    const SlpT & slp,
    const uint64_t pos_s, // position in the string
    const uint64_t pos_p // position in the pattern
)
{
    std::stack<typename SlpT::nodeT> path;
    std::stack<uint64_t> positions;
    getPrefixPathWithPositions(slp, path, positions, pos_s);

    assert(path.size() == positions.size());
    assert(positions.top() == 0);

    // let's go UP as much as we can
    auto x = path.top();
    path.pop();
    auto p = positions.top();
    positions.pop();

    uint64_t node = std::get<1>(x), last_from_left = node;
    bool direction = true;
    uint64_t prefix_hash = slp.subpattern_hash(pos_p, pos_p); // hash of the first prefix character
    //std :: cout << "ph: " << prefix_hash << std:: endl;
    uint64_t subpattern_hash = 0; // hash of the first prefix character
    uint64_t pattern_prefix_len = 1;

    //std :: cout << "node: " << node << " ph: " << prefix_hash << std:: endl;
    if (node != prefix_hash) return 0;
    //for (int i = 0; i < 10; i++) std::cout << slp.charAt(pos_s+i);
    //std::cout << std::endl;

    // ASCENT
    while (!path.empty()) {
        //std::cout << " pos: " << p  << " left len: " << slp.getNodeLen(left_sibling) << " right len: " << slp.getNodeLen(right_sibling) << std::endl; 
        //std::cout << "subpattern hash: " << p_hash << std::endl; 
        ///std::cout << "len: " << std::get<0>(x) << " hash: " << slp.get_hash(std::get<1>(x)) << std::endl;
        if (!direction) { // I come from the left child, the prefix can be extended
          last_from_left = node;
          uint64_t right_child = slp.getRight(node - slp.getAlphSize());
          uint64_t right_child_len = slp.getNodeLen(right_child);
          uint64_t right_hash = slp.getNodeHash(right_child);

          if ((pos_p + pattern_prefix_len + right_child_len) >= slp.getPatternLen()) break;
          //std::cout << "came from left: " << std::endl;
          prefix_hash = concat(prefix_hash, right_hash, pattern_prefix_len, slp.getAlphSize());
          subpattern_hash = slp.subpattern_hash(pos_p, pos_p+pattern_prefix_len+right_child_len-1);
          if (prefix_hash != subpattern_hash) break;
          pattern_prefix_len += right_child_len;

          //std::cout << "rl:" << right_child_len << "rh: " << right_hash << std::endl;
          //std::cout << "pattern pref len: " << pattern_prefix_len<< std::endl;
          //std::cout << " slp_prefix_hash: " << prefix_hash << " subpattern hash: " << subpattern_hash << std::endl;
          assert(prefix_hash == subpattern_hash);
        }

        //std::cout << "act len: " << std::get<0>(x) << " direction: " << direction
        //<< " pattern_prefix_len: " << pattern_prefix_len << std::endl;
        direction = std::get<2>(x);
        x = path.top();
        path.pop();
        p = positions.top();
        positions.pop();
        node = std::get<1>(x);
    }
    //std::cout << "act len: " << std::get<0>(x) << " dir: " << direction <<  " preflen:" << pattern_prefix_len << std::endl;

    if (last_from_left < slp.getAlphSize()) return pattern_prefix_len;
    if (!direction) node = slp.getRight(node - slp.getAlphSize());
    else            node = slp.getRight(last_from_left - slp.getAlphSize());
    // RE-DESCENT
    while (true) {
      if (node < slp.getAlphSize()) {
        uint64_t last_char = slp.subpattern_hash(pos_p+pattern_prefix_len, pos_p+pattern_prefix_len);
        pattern_prefix_len += (node == last_char);
        break;
      }
      //std::cout << "act len: " << slp.getNodeLen(node) <<  " preflen:" << pattern_prefix_len << std::endl;
      uint64_t left_child = slp.getLeft(node - slp.getAlphSize());
      uint64_t right_child = slp.getRight(node - slp.getAlphSize());
      uint64_t left_size = slp.getNodeLen(left_child);
      uint64_t left_hash = slp.getNodeHash(left_child);
      uint64_t candidate_hash = slp.subpattern_hash(pos_p+pattern_prefix_len, pos_p+pattern_prefix_len+left_size-1);
      //std::cout << "left_size: " << left_size << std::endl;
      //std::cout << "left_hash: " << left_hash << " prefix_hash: " << candidate_hash << std::endl;
      if (((pos_p + pattern_prefix_len + left_size - 1) > slp.getPatternLen()) || (left_hash != candidate_hash)) {
          //std::cout << "going left: " << std::endl;
          node = left_child;
      } else {
        node = right_child;
          //std::cout << "going right: " << std::endl;
        pattern_prefix_len += left_size;
      }
    }

    //std::cout << "ans: " << pattern_prefix_len << std::endl;
    return pattern_prefix_len;
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
