#include <vif.hpp>
#include <vif/unit_test.hpp>

using namespace vif;

int vif_main(int argc, char* argv[]) {
    vec2d v = gaussian_profile({{41,41}}, 4.0) + 0.1*gaussian_profile({{41,41}}, 10.0);
    vec2cd cv = fft(v);
    vec2d iv = ifft(cv)/v.size();

    for (uint_t i : range(v)) {
        check(v[i], iv[i]);
    }

    auto seed = make_seed(42);
    vec2d img = randomn(seed, 1000, 1000);
    vec2d psf = v;

    double st = now();
    vec2d cimg1 = convolve2d(img, psf);
    double fast = now() - st;
    st = now();
    vec2d cimg2 = convolve2d_naive(img, psf);
    double slow = now() - st;

    // TODO: investigate expected numerical precision of the FFT convolution algorithm
    for (uint_t i : range(v)) {
        check(cimg1[i], cimg2[i]);
    }

    print("fast version: ", fast);
    print("slow version: ", slow);

    return 0;
}
