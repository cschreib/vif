#ifndef VIF_INCLUDING_GENERIC_BITS
#error this file is not meant to be included separately, include "vif/utilty/generic.hpp" instead
#endif

namespace vif {
    // Generate linearly increasing values.
    namespace impl {
        template<typename T, typename ... Dims>
        auto indgen_(Dims&& ... ds) -> decltype(arr<T>(std::forward<Dims>(ds)...)) {
            auto v = arr<T>(std::forward<Dims>(ds)...);
            for (uint_t i = 0; i < v.size(); ++i) {
                v.safe[i] = i;
            }
            return v;
        }
    }

    template<typename ... Dims>
    auto findgen(Dims&& ... ds) -> decltype(impl::indgen_<float>(std::forward<Dims>(ds)...)) {
        return impl::indgen_<float>(std::forward<Dims>(ds)...);
    }

    template<typename ... Dims>
    auto dindgen(Dims&& ... ds) -> decltype(impl::indgen_<double>(std::forward<Dims>(ds)...)) {
        return impl::indgen_<double>(std::forward<Dims>(ds)...);
    }

    template<typename ... Dims>
    auto indgen(Dims&& ... ds) -> decltype(impl::indgen_<int_t>(std::forward<Dims>(ds)...)) {
        return impl::indgen_<int_t>(std::forward<Dims>(ds)...);
    }

    template<typename ... Dims>
    auto uindgen(Dims&& ... ds) -> decltype(impl::indgen_<uint_t>(std::forward<Dims>(ds)...)) {
        return impl::indgen_<uint_t>(std::forward<Dims>(ds)...);
    }
}
