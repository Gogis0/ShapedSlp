#include <iostream>
#include <string>
#include <queue>
#include "cmdline.h"
#include "Common.hpp"
#include "PlainSlp.hpp"
#include "PoSlp.hpp"
#include "ShapedSlp_Status.hpp"
#include "FixedBitLenCode.hpp"
#include "SelectType.hpp"
#include "VlcVec.hpp"

using namespace std;

using namespace std::chrono;
using timer = std::chrono::high_resolution_clock;

using var_t = uint32_t;
template
<
  class SlpT
  >
void measure_PlainSlp
(
 const NaiveSlp<var_t> & slp,
 std::string out,
 std::string &pattern
) {
  NaiveSlp<var_t> temp(slp);
  auto start = timer::now();
  temp.makeBinaryTree();
  SlpT pslp;
  pslp.init(temp);
  auto stop = timer::now();
  cout << "time to encode (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;
  pslp.printStatus();

  start = timer::now();
  ofstream fs(out);
  pslp.serialize(fs);
  stop = timer::now();
  cout << "time to serialize (ms): " << duration_cast<milliseconds>(stop-start).count() << endl;

  pslp.precompute_pattern(pattern);
  match_length_query(pslp, 90000, 3);
}


int main(int argc, char* argv[]) {
  using Fblc = FixedBitLenCode<>;

  cmdline::parser parser;
  parser.add<string>("input", 'i', "input file name. <input>.C and <input>.R in Navarro's RePair format", true);
  parser.add<string>("output", 'o', "output file to which data structure is written", true);
  parser.add<string>("format", 'f', "format of input: NavarroRepair | Bigrepair | rrepair", true);
  parser.parse_check(argc, argv);
  const string in = parser.get<string>("input");
  const string out = parser.get<string>("output");
  const string format = parser.get<string>("format");

  NaiveSlp<var_t> slp;
  if (format.compare("NavarroRepair") == 0) {
    slp.load_NavarroRepair(in.data());
  } else if (format.compare("Bigrepair") == 0) {
    slp.load_Bigrepair(in.data(),false);
  } else if (format.compare("rrepair") == 0) {
    slp.load_Bigrepair(in.data(),true);
  } else if (format.compare("solca") == 0) {
    slp.load_solca(in.data());
  } else {
    cerr << "error: specify a valid format name: NavarroRepair. Bigrepair" << endl;
    exit(1);
  }

  std::ifstream instream("./data/chr19.1.fa");
  std::stringstream buffer;
  buffer << instream.rdbuf();
  std::string s(buffer.str());


  std::string smol("GGGAACTTCTTCT");
  std::cout << s.substr(90000, 1) << std::endl;
  
  measure_PlainSlp<PlainSlp<var_t, Fblc, Fblc>>(slp, out, smol);
  return 0;
}
