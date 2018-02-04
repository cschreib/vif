#include <phypp.hpp>
#include <phypp/unit_test.hpp>

int main(int argc, char* argv[]) {
    vec1i v;

    auto get_and_check = [&]() {
        auto res = bounds(v, 4, 7);
        check(res[0], lower_bound(v, 4));
        check(res[1], upper_bound(v, 7));
        return res;
    };

    auto res = get_and_check();
    check(res[0], npos);
    check(res[1], npos);

    v = {0,1};

    res = get_and_check();
    check(res[0], 1);
    check(res[1], npos);

    v = {0,10};

    res = get_and_check();
    check(res[0], 0);
    check(res[1], 1);

    v = {10,20};

    res = get_and_check();
    check(res[0], npos);
    check(res[1], 0);

    v = {2,3,4,5,6,7,8,9};

    res = get_and_check();
    check(res[0], 2);
    check(res[1], 6);

    return 0;
}
