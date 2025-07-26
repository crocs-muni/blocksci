#ifndef blocksci_coinjoin_utils_hpp
#define blocksci_coinjoin_utils_hpp

#include <blocksci/blocksci_export.h>
#include <limits.h>

#include <algorithm>
#include <blocksci/address/address.hpp>
#include <unordered_set>
#include <vector>

namespace blocksci {
    class BLOCKSCI_EXPORT CoinjoinUtils {
        // port from dumplings
        static inline std::unordered_set<long> compute_ww2_denominations() {
            long maxSatoshis = 134375000000;
            long minSatoshis = 5000;
            long denom = 1l;

            std::unordered_set<long> denominations;
            denominations.insert(minSatoshis);
            denominations.insert(maxSatoshis);

            // Powers of 2
            for (int i = 0; i < maxSatoshis; i++) {
                denom *= 2;
                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(denom);
            }

            denom = 3;
            // Powers of 3
            for (int i = 0; i < INT_MAX; i++) {
                denom *= 3;

                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(denom);
            }

            // Powers of 3 * 2
            denom = 3;
            for (int i = 0; i < INT_MAX; i++) {
                denom *= 3;
                auto next = denom * 2;
                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(next);
            }

            denom = 10;
            // Powers of 10 (1-2-5 series)
            for (int i = 0; i < INT_MAX; i++) {
                denom *= 10;

                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(denom);
            }

            denom = 10;
            // Powers of 10 * 2 (1-2-5 series)
            for (int i = 0; i < INT_MAX; i++) {
                denom *= 10;
                auto a = denom * 2;
                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(a);
            }

            denom = 10;
            // Powers of 10 * 5 (1-2-5 series)
            for (int i = 0; i < INT_MAX; i++) {
                denom *= 10;
                auto a = denom * 2;

                if (denom < minSatoshis) {
                    continue;
                }

                if (denom > maxSatoshis) {
                    break;
                }

                denominations.insert(a);
            }

            return denominations;
        }

        static inline std::vector<std::string> compute_ww1_coord_scripts() {
            return {"bc1qs604c7jv6amk4cxqlnvuxv26hv3e48cds4m0ew", "bc1qa24tsgchvuxsaccp8vrnkfd85hrcpafg20kmjw"};
        }

       public:
        static inline std::vector<std::string> ww1_coord_scripts = CoinjoinUtils::compute_ww1_coord_scripts();
        static inline std::unordered_set<long> ww2_denominations = CoinjoinUtils::compute_ww2_denominations();

        CoinjoinUtils() = default;
        CoinjoinUtils(const CoinjoinUtils &other) = default;
        CoinjoinUtils(CoinjoinUtils &&other) = default;
        CoinjoinUtils &operator=(const CoinjoinUtils &other) = default;
        CoinjoinUtils &operator=(CoinjoinUtils &&other) = default;
        ~CoinjoinUtils() = default;
        static inline int64_t FirstWasabi2Block = 741213;
        static inline int64_t FirstWasabiBlock = 530500;
        static inline int64_t FirstSamouraiBlock = 570000;
        static inline int64_t FirstWasabiNoCoordAddressBlock = 610000;
        static inline int64_t FirstAshigaruBlock = 890000;  // a guess, but yea
    };

}  // namespace blocksci

#endif
