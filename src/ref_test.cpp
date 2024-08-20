#include <vector>    // for std::vector
#include <array>     // for std::array
#include <utility>   // for std::index_sequence
#include <cstddef>   // for std::size_t
#include <chrono>    // for std::chrono::high_resolution_clock
#include <iterator>  // for std::size
#include <print>     // for std::print
#include <execution>
#include <numeric>
#include <random>
#include <cmath>
#include <algorithm>
#include <string>
#include <stdexcept> // for std::runtime_error
#include <fstream>   // for std::ofstream
#include <format>    // for std::format
#include "codec.hpp"

extern "C" {
    #include <fec.h>
}

using uint = unsigned int;

template<
    uint dataSyms,
    uint paritySyms,
    uint symbolSize,
    uint firstConsecutiveRoot,
    uint galoisFieldPolynomial,
    uint primitiveElement,
    uint prepadSyms
>
using Codec = fec::rs::Codec<
    dataSyms,
    paritySyms,
    symbolSize,
    firstConsecutiveRoot,
    galoisFieldPolynomial,
    primitiveElement,
    prepadSyms>;

//- Struct to hold codec & test parameters
//
struct Params {
  uint symsize;
  uint genpoly;
  uint fcr;
  uint prim;
  uint nroots;
  uint ntrials;
};

//- Equality Function for checking parameter equality
bool operator==(const Params& lhs, const Params& rhs) {
  return lhs.symsize == rhs.symsize &&
         lhs.genpoly == rhs.genpoly &&
         lhs.fcr     == rhs.fcr     &&
         lhs.prim    == rhs.prim    &&
         lhs.nroots  == rhs.nroots  &&
         lhs.ntrials == rhs.ntrials;
}

//- Manual table of configurations for trialing RS Codecs
//
constexpr Params entries[] = {
  // symsize  genpoly  fcr  prim nroots ntrials
  //
  { 2,        0x7,       1,  1,   1,    1000 }, //  0 
  { 3,        0xb,       1,  1,   2,    1000 }, //  1 
  { 4,        0x13,      1,  1,   4,    1000 }, //  2 
  { 5,        0x25,      1,  1,   6,    1000 }, //  3 
  { 6,        0x43,      1,  1,   8,    1000 }, //  4 
  { 7,        0x89,      1,  1,  10,    1000 }, //  5 
  { 8,        0x11d,     1,  1,  32,    1000 }, //  6 
  { 8,        0x187,   112, 11,  32,    1000 }, //  7, Duplicates CCSDS codec
  { 9,        0x211,     1,  1,  32,    1000 }, //  8 
  {10,        0x409,     1,  1,  32,    1000 }, //  9 
  {11,        0x805,     1,  1,  32,    1000 }, //  A 
  {12,        0x1053,    1,  1,  32,    1000 }, //  B 
  {13,        0x201b,    1,  1,  32,    1000 }, //  C 
  {14,        0x4443,    1,  1,  32,    1000 }, //  D 
  {15,        0x8003,    1,  1,  32,    1000 }, //  E 
  {16,        0x1100b,   1,  1,  32,    1000 }  //  F 
};


//- Programmatically determine number of configurations
//
constexpr auto tab_size = std::size(entries);

//- Utility constexpr function to create compile-time array of configs
//
constexpr auto createStimulusTable() {
  std::array<Params, tab_size> arr{};

  for (std::size_t i = 0; i < tab_size; ++i) {
    arr[i] = entries[i];
  }

  return arr;
}

struct CodecMap {
  private:
    static constexpr auto mm0 = entries[0].symsize;
    static constexpr auto mm1 = entries[1].symsize;
    static constexpr auto mm2 = entries[2].symsize;
    static constexpr auto mm3 = entries[3].symsize;
    static constexpr auto mm4 = entries[4].symsize;
    static constexpr auto mm5 = entries[5].symsize;
    static constexpr auto mm6 = entries[6].symsize;
    static constexpr auto mm7 = entries[7].symsize;
    static constexpr auto mm8 = entries[8].symsize;
    static constexpr auto mm9 = entries[9].symsize;
    static constexpr auto mmA = entries[10].symsize;
    static constexpr auto mmB = entries[11].symsize;
    static constexpr auto mmC = entries[12].symsize;
    static constexpr auto mmD = entries[13].symsize;
    static constexpr auto mmE = entries[14].symsize;
    static constexpr auto mmF = entries[15].symsize;
    static constexpr auto nn0 = (1<<mm0) - 1;
    static constexpr auto nn1 = (1<<mm1) - 1;
    static constexpr auto nn2 = (1<<mm2) - 1;
    static constexpr auto nn3 = (1<<mm3) - 1;
    static constexpr auto nn4 = (1<<mm4) - 1;
    static constexpr auto nn5 = (1<<mm5) - 1;
    static constexpr auto nn6 = (1<<mm6) - 1;
    static constexpr auto nn7 = (1<<mm7) - 1;
    static constexpr auto nn8 = (1<<mm8) - 1;
    static constexpr auto nn9 = (1<<mm9) - 1;
    static constexpr auto nnA = (1<<mmA) - 1;
    static constexpr auto nnB = (1<<mmB) - 1;
    static constexpr auto nnC = (1<<mmC) - 1;
    static constexpr auto nnD = (1<<mmD) - 1;
    static constexpr auto nnE = (1<<mmE) - 1;
    static constexpr auto nnF = (1<<mmF) - 1;
    static constexpr auto p0  = entries[0].nroots;
    static constexpr auto p1  = entries[1].nroots;
    static constexpr auto p2  = entries[2].nroots;
    static constexpr auto p3  = entries[3].nroots;
    static constexpr auto p4  = entries[4].nroots;
    static constexpr auto p5  = entries[5].nroots;
    static constexpr auto p6  = entries[6].nroots;
    static constexpr auto p7  = entries[7].nroots;
    static constexpr auto p8  = entries[8].nroots;
    static constexpr auto p9  = entries[9].nroots;
    static constexpr auto pA  = entries[10].nroots;
    static constexpr auto pB  = entries[11].nroots;
    static constexpr auto pC  = entries[12].nroots;
    static constexpr auto pD  = entries[13].nroots;
    static constexpr auto pE  = entries[14].nroots;
    static constexpr auto pF  = entries[15].nroots;
    static constexpr auto fcr0 = entries[0].fcr;
    static constexpr auto fcr1 = entries[1].fcr;
    static constexpr auto fcr2 = entries[2].fcr;
    static constexpr auto fcr3 = entries[3].fcr;
    static constexpr auto fcr4 = entries[4].fcr;
    static constexpr auto fcr5 = entries[5].fcr;
    static constexpr auto fcr6 = entries[6].fcr;
    static constexpr auto fcr7 = entries[7].fcr;
    static constexpr auto fcr8 = entries[8].fcr;
    static constexpr auto fcr9 = entries[9].fcr;
    static constexpr auto fcrA = entries[10].fcr;
    static constexpr auto fcrB = entries[11].fcr;
    static constexpr auto fcrC = entries[12].fcr;
    static constexpr auto fcrD = entries[13].fcr;
    static constexpr auto fcrE = entries[14].fcr;
    static constexpr auto fcrF = entries[15].fcr;
    static constexpr auto gfp0 = entries[0].genpoly;
    static constexpr auto gfp1 = entries[1].genpoly;
    static constexpr auto gfp2 = entries[2].genpoly;
    static constexpr auto gfp3 = entries[3].genpoly;
    static constexpr auto gfp4 = entries[4].genpoly;
    static constexpr auto gfp5 = entries[5].genpoly;
    static constexpr auto gfp6 = entries[6].genpoly;
    static constexpr auto gfp7 = entries[7].genpoly;
    static constexpr auto gfp8 = entries[8].genpoly;
    static constexpr auto gfp9 = entries[9].genpoly;
    static constexpr auto gfpA = entries[10].genpoly;
    static constexpr auto gfpB = entries[11].genpoly;
    static constexpr auto gfpC = entries[12].genpoly;
    static constexpr auto gfpD = entries[13].genpoly;
    static constexpr auto gfpE = entries[14].genpoly;
    static constexpr auto gfpF = entries[15].genpoly;
    static constexpr auto prim0 = entries[0].prim;
    static constexpr auto prim1 = entries[1].prim;
    static constexpr auto prim2 = entries[2].prim;
    static constexpr auto prim3 = entries[3].prim;
    static constexpr auto prim4 = entries[4].prim;
    static constexpr auto prim5 = entries[5].prim;
    static constexpr auto prim6 = entries[6].prim;
    static constexpr auto prim7 = entries[7].prim;
    static constexpr auto prim8 = entries[8].prim;
    static constexpr auto prim9 = entries[9].prim;
    static constexpr auto primA = entries[10].prim;
    static constexpr auto primB = entries[11].prim;
    static constexpr auto primC = entries[12].prim;
    static constexpr auto primD = entries[13].prim;
    static constexpr auto primE = entries[14].prim;
    static constexpr auto primF = entries[15].prim;

  public:
    static constexpr auto params = createStimulusTable();
    static constexpr auto Codec0 = Codec<nn0 - p0, p0, mm0, fcr0, gfp0, prim0, 0>();
    static constexpr auto Codec1 = Codec<nn1 - p1, p1, mm1, fcr1, gfp1, prim1, 0>();
    static constexpr auto Codec2 = Codec<nn2 - p2, p2, mm2, fcr2, gfp2, prim2, 0>();
    static constexpr auto Codec3 = Codec<nn3 - p3, p3, mm3, fcr3, gfp3, prim3, 0>();
    static constexpr auto Codec4 = Codec<nn4 - p4, p4, mm4, fcr4, gfp4, prim4, 0>();
    static constexpr auto Codec5 = Codec<nn5 - p5, p5, mm5, fcr5, gfp5, prim5, 0>();
    static constexpr auto Codec6 = Codec<nn6 - p6, p6, mm6, fcr6, gfp6, prim6, 0>();
    static constexpr auto Codec7 = Codec<nn7 - p7, p7, mm7, fcr7, gfp7, prim7, 0>();
    static constexpr auto Codec8 = Codec<nn8 - p8, p8, mm8, fcr8, gfp8, prim8, 0>();
    static constexpr auto Codec9 = Codec<nn9 - p9, p9, mm9, fcr9, gfp9, prim9, 0>();
    static constexpr auto CodecA = Codec<nnA - pA, pA, mmA, fcrA, gfpA, primA, 0>();
    static constexpr auto CodecB = Codec<nnB - pB, pB, mmB, fcrB, gfpB, primB, 0>();
    static constexpr auto CodecC = Codec<nnC - pC, pC, mmC, fcrC, gfpC, primC, 0>();
    static constexpr auto CodecD = Codec<nnD - pD, pD, mmD, fcrD, gfpD, primD, 0>();
    static constexpr auto CodecE = Codec<nnE - pE, pE, mmE, fcrE, gfpE, primE, 0>();
    static constexpr auto CodecF = Codec<nnF - pF, pF, mmF, fcrF, gfpF, primF, 0>();

    static inline auto getIndex(const Params& p) {
      for (int i = 0; i < params.size(); ++i) {
        if (params[i] == p) {
          return i;
        }
      }

      throw std::runtime_error("No matching parameter!");
    }

    // NOTE: Don't actually use this lmao
    template<typename Container>
    static inline auto encode(const Params& p, Container& data) {
      auto idx = getIndex(p);

      return idx ==  0 ? Codec0.encode(data) :
             idx ==  1 ? Codec1.encode(data) :
             idx ==  2 ? Codec2.encode(data) :
             idx ==  3 ? Codec3.encode(data) :
             idx ==  4 ? Codec4.encode(data) :
             idx ==  5 ? Codec5.encode(data) :
             idx ==  6 ? Codec6.encode(data) :
             idx ==  7 ? Codec7.encode(data) :
             idx ==  8 ? Codec8.encode(data) :
             idx ==  9 ? Codec9.encode(data) :
             idx == 10 ? CodecA.encode(data) :
             idx == 11 ? CodecB.encode(data) :
             idx == 12 ? CodecC.encode(data) :
             idx == 13 ? CodecD.encode(data) :
             idx == 14 ? CodecE.encode(data) : CodecF.encode(data);
    }

    // NOTE: Don't actually use this lmao
    template<typename Container>
    static inline auto decode(const Params& p, Container& data) {
      auto idx = getIndex(p);

      return idx ==  0 ? Codec0.decode(data) :
             idx ==  1 ? Codec1.decode(data) :
             idx ==  2 ? Codec2.decode(data) :
             idx ==  3 ? Codec3.decode(data) :
             idx ==  4 ? Codec4.decode(data) :
             idx ==  5 ? Codec5.decode(data) :
             idx ==  6 ? Codec6.decode(data) :
             idx ==  7 ? Codec7.decode(data) :
             idx ==  8 ? Codec8.decode(data) :
             idx ==  9 ? Codec9.decode(data) :
             idx == 10 ? CodecA.decode(data) :
             idx == 11 ? CodecB.decode(data) :
             idx == 12 ? CodecC.decode(data) :
             idx == 13 ? CodecD.decode(data) :
             idx == 14 ? CodecE.decode(data) : CodecF.decode(data);
    }
};

template<typename SymbolType>
auto random_codeword(std::vector<SymbolType>& symbols, size_t numParity, size_t symbolSize) {
  // Only producing random data symbols
  auto numData = symbols.size() - numParity;

  // Create a bitmask to limit the random number to the specified symbol size in bits
  uint32_t bitmask = (1U << symbolSize) - 1;

  // Random number generator
  std::random_device rd;  // Non-deterministic random number generator
  std::mt19937 gen(rd()); // Mersenne Twister engine for random number generation

  // Fill the vector with random symbols by applying the bitmask
  for (size_t i = 0; i < numData; ++i) {
    symbols[i] = gen() & bitmask;
  }
}

auto rs_trial(const Params& p, std::vector<uint64_t>& durations) {
  // Setup
  //
  auto nn = (1<<p.symsize) - 1; // Codeword size
  auto kk = nn - p.nroots;      // data symbol count
  auto rs = init_rs_int(p.symsize,p.genpoly,p.fcr,p.prim,p.nroots,0);

  if(rs == nullptr){
    throw std::runtime_error("init_rs_int failed!");
  }

  std::vector<int32_t> codeword(nn, 0);
  random_codeword(codeword, p.nroots, p.symsize);

  for (auto i = 0; i < durations.size(); ++i) {
    // Generate fresh random codeword

    // Compute & append parity symbols
    // Provide codec, start of data, & start of parity
    encode_rs_int(rs, codeword.data(), codeword.data() + kk);

    // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    decode_rs_int(rs, codeword.data(), nullptr, 0);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
  }
};

auto rs_trial_ng(const Params& p, std::vector<uint64_t>& durations) {
  // Find the appropriate Codec and do the routine. Technically each Codec is a different type so
  // I can't simply do auto codec = findCodec() because I'm an idiot w/ a lack of foresight
  auto idx = CodecMap::getIndex(p);

  if (idx == 0) {
    auto& codec = CodecMap::Codec0;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 1) {
    auto& codec = CodecMap::Codec1;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 2) {
    auto& codec = CodecMap::Codec2;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 3) {
    auto& codec = CodecMap::Codec3;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 4) {
    auto& codec = CodecMap::Codec4;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 5) {
    auto& codec = CodecMap::Codec5;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 6) {
    auto& codec = CodecMap::Codec6;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 7) {
    auto& codec = CodecMap::Codec7;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 8) {
    auto& codec = CodecMap::Codec8;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 9) {
    auto& codec = CodecMap::Codec9;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 10) {
    auto& codec = CodecMap::CodecA;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 11) {
    auto& codec = CodecMap::CodecB;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 12) {
    auto& codec = CodecMap::CodecC;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 13) {
    auto& codec = CodecMap::CodecD;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 14) {
    auto& codec = CodecMap::CodecE;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
  else if (idx == 15) {
    auto& codec = CodecMap::CodecF;

    using symbol_t = typename decltype(codec.alpha_to)::value_type;

    std::vector<symbol_t> codeword(codec.NN, 0);
    random_codeword(codeword, codec.P, codec.MM);

    for (auto i = 0; i < durations.size(); ++i) {
      // Compute & append parity symbols
      codec.encode(codeword);

      // Perform timed decode and record duration
    auto start = std::chrono::high_resolution_clock::now();
    codec.decode(codeword);
    auto stop = std::chrono::high_resolution_clock::now();

    durations[i] = std::chrono::duration_cast<std::chrono::microseconds>(stop - start).count();
    }
  }
};

auto metrics(const std::vector<uint64_t>& durations) {
    auto [min_iter, max_iter] = std::minmax_element(
        std::execution::par_unseq,
        durations.begin(),
        durations.end()
    );
    auto sum = std::reduce(
        std::execution::par_unseq, 
        durations.begin(), 
        durations.end(), 
        0
    );
    auto mean = static_cast<float>(sum) / durations.size();
    auto variance = std::transform_reduce(
        std::execution::par_unseq, 
        durations.begin(), 
        durations.end(), 
        0.0, 
        std::plus<>(), 
        [mean](uint64_t x){return (x - mean) * (x - mean);}
    );
    auto stddev = std::sqrt(variance);

    std::print("Min Duration:     {} us\n"
               "Max Duration:     {} us\n"
               "Average Duration: {:.4f} us\n"
               "Std. Dev.:        {:.4f} us\n", *min_iter, *max_iter, mean, stddev);
}

inline void toBinFile(std::vector<std::uint64_t> const & data, int N, int K, std::string const & prefix) {
  std::string filename = std::format("durations_{}_{}_{}.bin", N, K, prefix);

  std::ofstream outFile(filename, std::ios::binary);

  outFile.write(reinterpret_cast<const char*>(data.data()), data.size() * sizeof(uint64_t));

  outFile.close();
}


int main() {
  constexpr auto codecs = CodecMap();

  constexpr uint mm    = 8;
  constexpr uint nn    = ((1<<mm) - 1);
  constexpr uint p     = 32;
  constexpr uint fcr   = 112;
  constexpr uint gfp   = 0x187;
  constexpr uint prim  = 11;
  constexpr auto codec = Codec<
      nn - p, // dataSyms (max codeword size - numroots)
      p,      // paritySyms
      mm,     // symbolSize
      fcr,    // firstConsecutiveRoot
      gfp,    // galoisFieldPolynomial
      prim,   // primitiveElement
      0       // prepadSyms
  >();

  auto entryNum = 0;

  std::vector<std::uint64_t> durations;

  for (const auto& entry : codecs.params) {
    std::print("---------- TEST PARAM #{} ----------\n"
        "Symbol Size (bits):     {}\n"
        "Generator Polynomial:   0x{:x}\n"
        "First Consecutive Root: {}\n"
        "Primitive Element:      {}\n"
        "Number of Roots:        {}\n"
        "Number of Trials:       {}\n"
        "------------------------------------\n",
        entryNum,
        entry.symsize,
        entry.genpoly,
        entry.fcr,
        entry.prim,
        entry.nroots,
        entry.ntrials);
    int nn = (1 << entry.symsize) - 1;
    int kk = nn - entry.nroots;
    std::print("Testing ({},{}) code...\n", nn, kk);

    //- Perform rs_trial(), will implicitly do ntrials and record durations
    //
    durations.resize(entry.ntrials);
    rs_trial(entry, durations);

    std::print("**********************\n");
    std::print("*** LEGACY METRICS ***\n");
    std::print("**********************\n");
    metrics(durations);
    std::print("**********************\n");

    //- Write durations to a file for follow-on histogram plotting
    //
    // TODO
    toBinFile(durations, nn, kk, "legacy");

    //- Perform rs_trial_ng entry.ntrials number of times, measure & collect each duration
    //
    rs_trial_ng(entry, durations);

    std::print("**********************\n");
    std::print("***   NG METRICS   ***\n");
    std::print("**********************\n");
    metrics(durations);
    std::print("**********************\n");


    //- Write durations to a file for follow-on histogram plotting
    //
    // TODO
    toBinFile(durations, nn, kk, "ng");
    entryNum++;
  }

  return 0;
}

