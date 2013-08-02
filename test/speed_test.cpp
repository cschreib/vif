#include <phypp.hpp>

template<typename T>
void prof(const std::string& msg, T&& t) {
    double d = profile(std::forward<T>(t));
    print(msg, ": ", date_str(d));
}

int main(int argc, char* argv[]) {
    const std::string& std = " (std)";
    const std::string& phypp = " (phypp)";

    vec1i u = indgen(1e7);
    vec1i v1 = intarr(1e7);
    vec1i v2 = intarr(1e7);
    vec1i v3;
    uint_t n = u.size();
    prof("Simple loop"+std, [&]() {
        for (uint_t i = 0; i < n; ++i) {
            v1.data[i] = 2*u.data[i];
        }
    });
    prof("Simple loop"+phypp, [&]() {
        for (uint_t i = 0; i < n; ++i) {
            v2[i] = 2*u[i];
        }
    });
    prof("Vectorized loop"+phypp, [&]() {
        v3 = 2*u;
    });

    return 0;
}
