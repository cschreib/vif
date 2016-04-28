#ifndef PHYPP_COLOR_HPP
#define PHYPP_COLOR_HPP

// Structure holding color data
struct rgb {
    float r = 0, g = 0, b = 0;

    rgb() = default;
    rgb(float r_, float g_, float b_) : r(r_), g(g_), b(b_) {}

    #define OPERATOR(name) \
        void operator name (const rgb& c) { \
            r name c.r; g name c.g; b name c.b; \
        } \
        void operator name (float f) { \
            r name f; g name f; b name f; \
        }

    OPERATOR(+=)
    OPERATOR(-=)
    OPERATOR(*=)
    OPERATOR(/=)

    #undef OPERATOR
};

#define OPERATOR(name) \
    rgb operator name (const rgb& c1, const rgb& c2) { \
        return {c1.r name c2.r, c1.g name c2.g, c1.b name c2.b}; \
    } \
    rgb operator name (const rgb& c1, float f) { \
        return {c1.r name f, c1.g name f, c1.b name f}; \
    } \
    rgb operator name (float f, const rgb& c2) { \
        return {f name c2.r, f name c2.g, f name c2.b}; \
    }

OPERATOR(+)
OPERATOR(-)
OPERATOR(*)
OPERATOR(/)

#undef OPERATOR

#endif
