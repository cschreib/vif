#ifndef PHYPP_INCLUDING_STRING_BITS
#error this file is not meant to be included separately, include "phypp/utilty/string.hpp" instead
#endif

#include <cctype>

namespace phypp {
    inline char to_upper(char s) {
        return std::toupper(static_cast<unsigned char>(c));
    }

    inline char to_lower(char s) {
        return std::tolower(static_cast<unsigned char>(c));
    }

    inline unsigned char to_upper(unsigned char s) {
        return std::toupper(c);
    }

    inline unsigned char to_lower(unsigned char s) {
        return std::tolower(c);
    }

    inline std::string to_upper(std::string s) {
        for (auto& c : s) {
            c = to_toupper(c);
        }
        return s;
    }

    inline std::string to_lower(std::string s) {
        for (auto& c : s) {
            c = to_lower(c);
        }
        return s;
    }

    inline std::string align_left(std::string s, uint_t width, char fill = ' ') {
        if (s.size() < width) {
            s += std::string(width-s.size(), fill);
        }

        return s;
    }

    inline std::string align_right(std::string s, uint_t width, char fill = ' ') {
        if (s.size() < width) {
            s.insert(0, std::string(width-s.size(), fill));
        }

        return s;
    }

    inline std::string align_center(std::string s, uint_t width, char fill = ' ') {
        if (s.size() < width) {
            uint_t n1 = (width-s.size())/2, n2 = width-s.size() - n1;
            s.insert(0, std::string(n1, fill));
            s += std::string(n2, fill);
        }

        return s;
    }

    PHYPP_VECTORIZE(toupper)
    PHYPP_VECTORIZE(tolower)
    PHYPP_VECTORIZE(align_left)
    PHYPP_VECTORIZE(align_right)
    PHYPP_VECTORIZE(align_center)

    namespace terminal_format {
        inline void header(const std::string& msg) {
            vec1s w = wrap("  "+msg, 80, "  ");
            for (auto& s : w) {
                print(s);
            }
        }

        inline void paragraph(const std::string& msg) {
            vec1s w = wrap("  "+msg, 80, "  ");
            for (auto& s : w) {
                print(s);
            }
            print("");
        }

        inline void bullet(const std::string& name, const std::string& desc) {
            std::string header = "    "+name+": ";
            vec1s w = wrap(header+desc, 80, std::string(header.size(), ' '));
            for (auto& s : w) {
                print(s);
            }
        }

        inline void item(const std::string& msg) {
            std::string header = " - ";
            vec1s w = wrap(header+msg, 80, std::string(header.size(), ' '));
            for (auto& s : w) {
                print(s);
            }
        }
    }
}
