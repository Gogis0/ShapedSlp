#ifndef INCLUDE_GUARD_PlainSlp
#define INCLUDE_GUARD_PlainSlp

#include <stdint.h> // include uint64_t etc.
#include <iostream>
#include <string>
#include <map>
#include <set>
#include <stack>
#include "Common.hpp"
#include "NaiveSlp.hpp"
#include "fbtree.h"


template
<
    typename tparam_var_t,
    class VarVecT,
    class LenVecT
    >
class PlainSlp
{
public:
    //// Public constant, alias etc.
    using var_t = tparam_var_t;
    using nodeT = std::tuple<uint64_t, var_t, var_t>; // (expansion_length, node_id, child_rank): stack of nodes indicates a path

private:
    size_t len_;
    size_t numRules_;
    uint64_t hash_;
    var_t startVar_;
    std::vector<char> alph_;
    VarVecT left_;
    VarVecT right_;
    LenVecT expLen_; // expansion lengths
    LenVecT kr_hash_; // KR fingerprints
    std::vector<uint64_t> prefix_hash_; // hash for the pattern suffixes
    std::vector<uint64_t> base_powers_; // prime powers

public:
    PlainSlp () : len_(0),
        numRules_(0),
        startVar_(0)
    {}


    ~PlainSlp()
    {}

    uint64_t getHashWhole() const {
      return kr_hash_[kr_hash_.size() - 1];
    }

    uint64_t getPatternLen() const {
      return prefix_hash_.size();
    }

    uint64_t getBasePower(const size_t i) const {
      return base_powers_[i];
    }

    uint64_t getNodeHash(var_t varId) const {
      return (varId < getAlphSize()) ? hash_char(varId) : kr_hash_[varId - getAlphSize()];
    }

    uint64_t getNodeLen(var_t varId) const {
      return (varId < getAlphSize()) ? 1 : expLen_[varId - getAlphSize()];
    }

    size_t getLen() const
    {
        return len_;
    }


    size_t getLenSeq() const
    {
        return 1;
    }


    size_t getNumRules() const
    {
        return numRules_;
    }


    size_t getAlphSize() const
    {
        return alph_.size();
    }


    char getChar(uint64_t i) const
    {
        return alph_[i];
    }


    var_t getLeft
    (
        const uint64_t ruleId
    ) const
    {
        return left_[ruleId];
    }


    var_t getRight
    (
        const uint64_t ruleId
    ) const
    {
        return right_[ruleId];
    }

    uint64_t get_hash
    (
        const uint64_t ruleId
    ) const
    {
        return kr_hash_[ruleId];
    }


    var_t getSeq(uint64_t i) const
    {
        return startVar_;
    }


    size_t getNumRulesOfSlp() const
    {
        return numRules_;
    }


    char charAt
    (
        const uint64_t pos //!< 0-based position
    ) const
    {
        return charAt(pos, startVar_);
    }

    char charAt
    (
        const uint64_t pos, //!< relative position in a variable
        const uint64_t varId //!< varId
    ) const
    {
        // std::cout << "pos = " << pos << ", len = " << len << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;

        if (varId < getAlphSize()) {
            return getChar(varId);
        }
        const uint64_t leftVarId = getLeft(varId - getAlphSize());
        const uint64_t leftLen = (leftVarId < getAlphSize()) ? 1 : expLen_[leftVarId - getAlphSize()];
        if (pos < leftLen) {
            return charAt(pos, leftVarId);
        } else {
            return charAt(pos - leftLen, getRight(varId - getAlphSize()));
        }
    }

    void expandSubstr
    (
        const uint64_t pos, //!< beginning position
        const uint64_t lenExpand, //!< length to expand
        char * str //!< [out] must have length at least 'len'
    ) const
    {
        assert(pos < getLen());
        assert(lenExpand > 0);
        assert(lenExpand <= getLen() - pos);

        expandSubstr(pos, lenExpand, str, startVar_);
    }


    void expandSubstr
    (
        const uint64_t pos, //!< relative position in a variable
        const uint64_t lenExpand, //!< length to expand
        char * str, //!< [out] must have length at least 'len'
        const uint64_t varId //!< varId
    ) const
    {
        // std::cout << "pos = " << pos << ", len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
        assert(lenExpand > 0);

        if (varId < getAlphSize()) {
            *str = getChar(varId);
            return;
        }
        const uint64_t leftVarId = getLeft(varId - getAlphSize());
        const uint64_t leftLen = (leftVarId < getAlphSize()) ? 1 : expLen_[leftVarId - getAlphSize()];
        if (pos < leftLen) {
            expandSubstr(pos, lenExpand, str, leftVarId);
            if (leftLen - pos < lenExpand) {
                expandPref(lenExpand - (leftLen - pos), str + (leftLen - pos), getRight(varId - getAlphSize()));
            }
        } else {
            expandSubstr(pos - leftLen, lenExpand, str, getRight(varId - getAlphSize()));
        }
    }


    void expandPref
    (
        const uint64_t lenExpand, //!< length to expand
        char * str, //!< [out] must have length at least 'len'
        const uint64_t varId //!< varId
    ) const
    {
        // std::cout << "len = " << len << ", varLen = " << varLen << ", stgOffset = " << stgOffset << ", slpOffset = " << slpOffset << std::endl;
        assert(lenExpand > 0);

        if (varId < getAlphSize()) {
            *str = getChar(varId);
            return;
        }
        const uint64_t leftVarId = getLeft(varId - getAlphSize());
        const uint64_t leftLen = (leftVarId < getAlphSize()) ? 1 : expLen_[leftVarId - getAlphSize()];
        expandPref(lenExpand, str, leftVarId);
        if (lenExpand > leftLen) {
            expandPref(lenExpand - leftLen, str + leftLen, getRight(varId - getAlphSize()));
        }
    }


    nodeT getRootNode() const
    {
        return std::forward_as_tuple(getLen(), 0, 0);
    }


    nodeT getChildNode_Root
    (
        const uint64_t idx
    ) const
    {
        return std::forward_as_tuple(getLen(), startVar_, idx);
    }


    nodeT getChildNode
    (
        const nodeT & node,
        const uint64_t idx
    ) const
    {
        assert(std::get<0>(node) > 1); // len > 1

        const uint64_t len = std::get<0>(node);
        const uint64_t varId = std::get<1>(node);
        const uint64_t newVarId = (idx == 0) ? getLeft(varId - getAlphSize()) : getRight(varId - getAlphSize());
        const uint64_t newLen = (newVarId < getAlphSize()) ? 1 : expLen_[newVarId - getAlphSize()];
        return std::forward_as_tuple(newLen, newVarId, idx);
    }


    nodeT getChildNodeForPos_Root
    (
        uint64_t & pos //! [in, out]
    ) const
    {
        return std::forward_as_tuple(getLen(), startVar_, 0);
    }


    nodeT getChildNodeForPos
    (
        const nodeT & node,
        uint64_t & pos //! [in, out]
    ) const
    {
        assert(std::get<0>(node) > 1); // len > 1

        const uint64_t len = std::get<0>(node);
        const uint64_t varId = std::get<1>(node);
        const uint64_t leftVarId = getLeft(varId - getAlphSize());
        const uint64_t leftLen = (leftVarId < getAlphSize()) ? 1 : expLen_[leftVarId - getAlphSize()];
        const uint64_t idx = (pos < leftLen) ? 0 : 1;
        const uint64_t newVarId = (idx == 0) ? leftVarId : getRight(varId - getAlphSize());
        const uint64_t newLen = (idx == 0) ? leftLen : len - leftLen;
        pos -= leftLen * idx;
        return std::forward_as_tuple(newLen, newVarId, idx);
    }


    void printStatus
    (
        const bool verbose = false
    ) const
    {
        std::cout << "PlainSlp object (" << this << ") " << __func__ << "(" << verbose << ") BEGIN" << std::endl;
        std::cout << "alphSize = " << getAlphSize() << ", len = " << getLen() << ", lenSeq = " << getLenSeq()
                  << ", numRulesOfSlp = " << getNumRulesOfSlp()
                  << std::endl;
        const size_t bytesAlph = sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
        const size_t bytesLeft = left_.calcMemBytes();
        const size_t bytesRight = right_.calcMemBytes();
        const size_t bytesExpLen = expLen_.calcMemBytes();
        const size_t bytesHashVec = kr_hash_.calcMemBytes();
        std::cout << "DS sizes (bytes)" << std::endl;
        std::cout << bytesAlph + bytesLeft + bytesRight + bytesExpLen << std::endl;
        std::cout << "| alph = " << bytesAlph << std::endl;
        std::cout << "| left = " << bytesLeft << std::endl;
        std::cout << "| right = " << bytesRight << std::endl;
        std::cout << "| expLen = " << bytesExpLen << std::endl;
        std::cout << "| hashVec: " << bytesHashVec << std::endl;
        if (verbose) {
            std::cout << "rules: " << std::endl;
            for (uint64_t i = 0; i < getNumRulesOfSlp(); ++i) {
                const auto l = getLeft(i);
                const auto r = getRight(i);
                std::cout << i + getAlphSize() << "(";
                if (l < getAlphSize()) {
                    std::cout << "'" << getChar(l) << "'";
                } else {
                    std::cout << l;
                }
                std::cout << ", ";
                if (r < getAlphSize()) {
                    std::cout << "'" << getChar(r) << "'";
                } else {
                    std::cout << r;
                }
                std::cout << ") ";
            }
            std::cout << std::endl;
            std::cout << "expLen: " << std::endl;
            for (uint64_t i = 0; i < getNumRulesOfSlp(); ++i) {
                std::cout << expLen_[i] << " ";
            }
            std::cout << std::endl;
        }
        std::cout << "PlainSlp object (" << this << ") " << __func__ << "(" << verbose << ") END" << std::endl;
    }


    size_t calcMemBytes() const
    {
        size_t ret = 0;
        ret += sizeof(std::vector<char>) + (sizeof(char) * alph_.size());
        ret += left_.calcMemBytes();
        ret += right_.calcMemBytes();
        ret += expLen_.calcMemBytes();
        return ret;
    }


    void load
    (
        std::istream & in
    )
    {
        in.read((char*) & len_, sizeof(len_));
        in.read((char*) & numRules_, sizeof(numRules_));
        in.read((char*) & startVar_, sizeof(startVar_));
        uint64_t alphSize = 0;
        in.read((char*) & alphSize, sizeof(alphSize));
        alph_.resize(alphSize);
        in.read((char*) alph_.data(), alphSize * sizeof(alph_[0]));
        left_.load(in);
        right_.load(in);
        expLen_.load(in);
        kr_hash_.load(in);
    }


    void serialize
    (
        std::ostream & out
    ) const
    {
        out.write((char*) & len_, sizeof(len_));
        out.write((char*) & numRules_, sizeof(numRules_));
        out.write((char*) & startVar_, sizeof(startVar_));
        uint64_t alphSize = getAlphSize();
        out.write((char*) & alphSize, sizeof(alphSize));
        out.write((char*) alph_.data(), alphSize * sizeof(alph_[0]));
        left_.serialize(out);
        right_.serialize(out);
        expLen_.serialize(out);
        kr_hash_.serialize(out);
    }

    void precompute_pattern(const std::string &pattern) {
        precompute_pattern_hashes(pattern, this->prefix_hash_, this->base_powers_, base);
    }

    uint64_t subpattern_hash(const uint64_t i, const uint64_t j) const {
        return compute_subpattern_hash(this->prefix_hash_, this->base_powers_, i, j);
    }


    void init
    (
        const NaiveSlp<var_t> & slp
    )
    {
        assert(slp.getLenSeq() == 1);

        const uint64_t alphSize = slp.getAlphSize();
        alph_.resize(alphSize);
        for (uint64_t i = 0; i < alphSize; ++i) {
            alph_[i] = slp.getChar(i);
        }
        numRules_ = slp.getNumRules();
        startVar_ = numRules_ - 1 + getAlphSize();
        {
            std::vector<uint64_t> left(slp.getNumRules());
            std::vector<uint64_t> right(slp.getNumRules());
            for (uint64_t i = 0; i < slp.getNumRules(); ++i) {
                left[i] = slp.getLeft(i);
                right[i] = slp.getRight(i);
            }
            left_.init(left);
            right_.init(right);
        }
        {
            std::vector<uint64_t> expLen(slp.getNumRules());
            slp.makeLenVec(expLen);
            len_ = expLen[slp.getNumRules() - 1];
            expLen_.init(expLen);

            std::vector<uint64_t> hashVec(slp.getNumRules());
            slp.makeHashVec(expLen, hashVec);
            hash_ = hashVec[slp.getNumRules() - 1];

            kr_hash_.init(hashVec);
        }
    }
};

#endif
