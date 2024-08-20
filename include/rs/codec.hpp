#ifndef CODEC_HPP
#define CODEC_HPP

#include <cstddef>
#include <type_traits>
#include <cstdint>
#include <array>
#include <stdexcept>
#include <type_traits>
#include <algorithm>
#include <execution>

void encode();
void decode();

namespace fec::rs {

  /*
  inline constexpr unsigned int modnn(unsigned int x, unsigned int NN, unsigned int MM) {
    while (x >= NN) {
      x -= NN;
      x  = (x >> MM) + (x & NN);
    }

    return x;
  }
  */

  template<
    unsigned int dataSyms, 
    unsigned int paritySyms, 
    unsigned int symbolSize,
    unsigned int firstConsecutiveRoot,
    unsigned int galoisFieldPolynomial,
    unsigned int primitiveElement,
    unsigned int prepadSyms
  >
  class Codec {
    private:
      static constexpr auto find_iprim(unsigned int pe, unsigned int nn)
      {
        // Find the prim-th root of 1, used in decoding
        unsigned int ipe = 1;
        for (ipe = 1; (ipe % pe) != 0; ipe += nn)
          ;

        return ipe/pe;
      }

    public:
      // Codec Parameters
      static constexpr std::size_t  MM      = symbolSize;

      // Determine the symbol type based on the number of bits per symbol
      using symbol_t = typename std::conditional<
                           (MM <= 8),
                           std::uint8_t,
                           typename std::conditional<
                               (MM <= 16),
                               std::uint16_t,
                               typename std::conditional<
                                   (MM <= 32),
                                   std::uint32_t,
                                   typename std::conditional<
                                       (MM <= 64),
                                       std::uint64_t,
                                       std::intmax_t
                                   >::type // Unreasonably large, lol
                               >::type     // Word sized
                           >::type         // Half-word sized
                       >::type;            // Byte sized

      static constexpr symbol_t     KK      = dataSyms;
      static constexpr symbol_t     P       = paritySyms;
      static constexpr symbol_t     maxSyms = (1UL << MM) - 1;
      static constexpr symbol_t     NN      = maxSyms;
      static constexpr symbol_t     fcr     = firstConsecutiveRoot;
      static constexpr symbol_t     gfPoly  = galoisFieldPolynomial;
      static constexpr symbol_t     prim    = primitiveElement;
      static constexpr symbol_t     iprim   = find_iprim(prim, NN);
      static constexpr symbol_t     pad     = prepadSyms;

      // Ensure Parameter Validity
      static_assert(KK > 0);              // Positive number of data symbols
      static_assert(P  > 0);              // Positive number of parity symbols
      static_assert((KK + P) <= maxSyms); // Symbol count cannot exceed max block size
      static_assert(pad < KK);            // Padding symbols can't exceed data
      static_assert(MM <= sizeof(std::intmax_t) * 8); // (Unreasonable) Max symbol size
      static_assert(fcr <= NN);
      static_assert(prim <= NN);


      // Look-up Tables
      /* static constexpr */ std::array<symbol_t, NN+1> alpha_to/*{}*/;
      /* static constexpr */ std::array<symbol_t, NN+1> index_of/*{}*/;
      /* static constexpr */ std::array<symbol_t,  P+1> genpoly/*{}*/;


      constexpr Codec<
        dataSyms,
        paritySyms,
        symbolSize,
        firstConsecutiveRoot,
        galoisFieldPolynomial,
        primitiveElement,
        prepadSyms>() 
        : alpha_to{}
        , index_of{}
        , genpoly{} {
          /* Generate Galois Field Lookup Tables */
          index_of[0] = NN;
          alpha_to[NN] = 0;

          symbol_t sr = 1;

          for (symbol_t i = 0; i < NN; ++i) {
            index_of[sr] = i;
            alpha_to[i] = sr;
            sr <<= 1;
            if (sr & (NN + 1)) {
              sr ^= gfPoly;
            }
            sr &= NN;
          }

          // TODO: Figure out way to detect this at compile time
          // static_assert(sr == 1);

          /* Form RS code generator polynomial from its roots */
          for (symbol_t i = 0, root = fcr*prim; i < P; ++i, root += prim) {
            genpoly[i+1] = 1;

            /* Multiply genpoly[] by 2**(root + x) */
            for (symbol_t j = i; j > 0; --j) {
              if (genpoly[j]) {
                genpoly[j] = genpoly[j-1] ^ alpha_to[modnn(index_of[genpoly[j]] + root)];
              }
              else {
                genpoly[j] = genpoly[j-1];
              }
            }

            /* genpoly[0] can never be zero */
            genpoly[0] = alpha_to[modnn(index_of[genpoly[0]] + root)];
          }

          /* Convert genpoly[] to index form for quicker encoding */
          for (symbol_t i = 0; i <= P; ++i) {
            genpoly[i] = index_of[genpoly[i]];
          }
      }

      template<typename Container>
      void encode(Container& buffer) const {
        static_assert(std::is_same_v<typename Container::value_type, symbol_t>,
                      "Container must hold elements of symbol_t");

        if (buffer.size() < (KK + P)) {
          throw std::invalid_argument(
              "Insufficient buffer size to hold data + parity symbols"
          );
        }

        // Iterators to start of data/parity symbol portions of buffer
        auto data_begin = buffer.begin();
        auto parity_begin = buffer.begin() + KK;

        // Zero out parity symbols
        std::fill(parity_begin, buffer.end(), symbol_t{0});

        for (symbol_t i = 0; i < (NN - P - pad); ++i) {
          symbol_t feedback = index_of[data_begin[i] ^ parity_begin[0]];

          if (feedback != NN) { // feedback term is non-zero
            for (symbol_t j = 1; j < P; ++j) {
              parity_begin[j] = alpha_to[modnn(feedback + genpoly[P-j])];
            }
          }
          // Shift parity symbols
          std::shift_left(parity_begin, buffer.end(), 1);
          if (feedback != NN) {
            parity_begin[P-1] = alpha_to[modnn(feedback + genpoly[0])];
          }
          else {
            parity_begin[P-1] = 0;
          }
        }
      }

      template<typename Container>
      std::intmax_t decode(Container& buffer) const {
        // TODO: SHOULD also return symbol error locations
        // NOTE: leaving out erasure detection
        
        // Form the syndromes. Evalute buffer(x) @ roots of g(x)
        std::array<symbol_t, P> syndromes{};
        std::fill(syndromes.begin(), syndromes.end(), buffer[0]);

        for (symbol_t j = 1; j < (NN - pad); ++j) {
          for (symbol_t i = 0; i < P; ++i) {
            if (syndromes[i] == 0) {
              syndromes[i] = buffer[j];
            }
            else {
              syndromes[i]  = buffer[j];
              syndromes[i] ^= alpha_to[modnn(index_of[syndromes[i]] + (fcr+i)*prim)];
            }
          }
        }

        // Convert syndromes to index form, checking for nonzero condition
        symbol_t syn_error{0};
        for (symbol_t i = 0; i < P; ++i) {
          syn_error |= syndromes[i];
          syndromes[i] = index_of[syndromes[i]];
        }

        if (!syn_error) {
          // zero errors
          return 0;
        }
        std::array<symbol_t,P+1> lambda{}; // Error/Erasure polynomial
        std::fill(lambda.begin() + 1, lambda.end(), symbol_t{0});
        lambda[0] = 1;

        // NOTE: initializing lambda to be the erasure locator poly would go here

        std::array<symbol_t,P+1> beta{};
        std::array<symbol_t,P+1> tau{};
        std::array<symbol_t,P+1> omega{};
        for (symbol_t i = 0; i < (P + 1); ++i) {
          beta[i] = index_of[lambda[i]];
        }

        // Berlekamp-Massey
        // NOTE: erasure detection would slightly change some of this
        constexpr symbol_t no_eras = 0;
        symbol_t r  = 0; // step number
        symbol_t el = 0;

        while (++r <= P) {
          // Compute discrepancy at r'th step in polynomial form
          symbol_t discr_r = 0;
          for (symbol_t i = 0; i < r; ++i) {
            if ((lambda[i] != 0) && (syndromes[r-i-1] != NN)) {
              discr_r ^= alpha_to[modnn(index_of[lambda[i]] + syndromes[r-i-1])];
            }
          }

          discr_r = index_of[discr_r]; // Index form
          if (discr_r == NN) {
            std::shift_right(beta.begin(), beta.end(), 1);
            beta[0] = NN;
          }
          else {
            tau[0] = lambda[0];
            for (symbol_t i = 0; i < P; ++i) {
              if (beta[i] != NN) {
                tau[i+1] = lambda[i+1] ^ alpha_to[modnn(discr_r + beta[i])];
              }
              else {
                tau[i+1] = lambda[i+1];
              }
            }
            if ( (2*el) <= (r + no_eras - 1) ) {
              el = r + no_eras - el;

              for (symbol_t i = 0; i <= P; ++i) {
                if (lambda[i] == 0) {
                  beta[i] = NN;
                }
                else {
                  beta[i] = modnn(index_of[lambda[i]] - discr_r + NN);
                }
              }
            }
            else {
              std::shift_right(beta.begin(), beta.end(), 1);
              beta[0] = NN;
            }
            std::copy(
                std::execution::par_unseq,
                tau.begin(), tau.end(),
                lambda.begin()
            );
          }
        } // end of while loop

        // Convert lambda to index form and compute degree of lambda
        symbol_t deg_lambda = 0;
        for (symbol_t i = 0; i < (P + 1); ++i) {
          lambda[i] = index_of[lambda[i]];
          if (lambda[i] != NN) {
            deg_lambda = i;
          }
        }

        // Find roots of error poly via Chien Serch
        std::array<symbol_t, P+1> reg{};
        std::copy(
            std::execution::par_unseq,
            lambda.begin() + 1, lambda.end(),
            reg.begin() + 1
        );
        std::intmax_t count = 0; // Number of roots of lambda

        std::array<symbol_t, P> root{};
        std::array<symbol_t, P> loc{};

        for (symbol_t i = 1, k = iprim-1; i <= NN; ++i, k = modnn(k+iprim)) {
          symbol_t q = 1; // lambda[0] is always 0
          for (symbol_t j = deg_lambda; j > 0; --j) {
            if (reg[j] != NN) {
              reg[j] = modnn(reg[j] + j);
              q ^= alpha_to[reg[j]];
            }
          }
          if (q != 0) continue; // Not a root

          // Store root in index form and error loc number
          root[count] = i;
          loc[count] = k;

          // If max roots found, abort to save cycles
          if (++count == deg_lambda) break;
        }
        if (deg_lambda != count) {
          // uncorrectable error detected
          return -1;
        }
        // Compute error evaluator polynomial omega = syndrome*lambda (mod x**P)
        symbol_t deg_omega = deg_lambda-1;
        for (symbol_t i = 0; i <= deg_omega; ++i) {
          symbol_t tmp = 0;
          for (symbol_t j = i; j >= 0; --j) {
            if ((syndromes[i - j] != NN) && (lambda[j] != NN)) {     // TODO: RESOLVE INFINITE LOOP HERE, J NEVER GOES BELOW ZERO BECAUSE UNSIGNED
              tmp ^= alpha_to[modnn(syndromes[i-j] + lambda[j])];
            }
            // TODO: DEBUGGING UNSIGNED TO CATCH ROLL OVER
            if (j == 0) break;
          }
          omega[i] = index_of[tmp];
        }

        // Compute error values in polynomial form.
        // num1 = omega(inv(X(l)))
        // num2 = inv(X(l))**(firstConsecRoot-1)
        // den = lambda_pr(inv(X(l)))
        for (symbol_t j = count - 1; j >= 0; --j) {
          symbol_t num1 = 0;
          for (symbol_t i = deg_omega; i >= 0; --i) {
            if (omega[i] != NN) {
              num1 ^= alpha_to[modnn(omega[i] + i * root[j])];
            }
            // TODO: DEBUGGING UNSIGNED TO CATCH ROLL OVER
            if (i == 0) break;
          }
          symbol_t num2 = alpha_to[modnn(root[j] * (fcr - 1) + NN)];
          symbol_t den = 0;

          // lambda[i+1] for i even is the formal deriv. lambda_pr of lambda[i]
          for (symbol_t i = std::min(deg_lambda, static_cast<symbol_t>(P-1)) & ~1; i >= 0; i -=2) {
            if (lambda[i+1] != NN) {
              den ^= alpha_to[modnn(lambda[i+1] + i * root[j])];
            }
            // TODO: DEBUGGING UNSIGNED TO CATCH ROLL OVER
            if (i == 0) break;
          }
          // Apply error to data
          if (num1 != 0 && loc[j] >= pad) {
            buffer[loc[j]-pad] ^= alpha_to[
              modnn(index_of[num1] + index_of[num2] + NN - index_of[den])
            ];
          }
          // TODO: DEBUGGING UNSIGNED TO CATCH ROLL OVER
          if (j == 0) break;
        }

        /*
         * TODO: Implement if we actually want erasure logic
        if (eras_pos) {
          for (symbol_t i = 0; i < count; ++i) {
            eras_pos[i] = loc[i];
          }
        }
        */

        return count;
      }
    private:
      inline constexpr symbol_t modnn(symbol_t x) const{
        while(x >= NN) {
          x -= NN;
          x = (x >> MM) + (x & NN);
        }

        return x;
      }

  }; // End class Codec
} // End namespace fec::rs

#endif // CODEC_HPP

