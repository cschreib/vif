#ifndef VIF_INCLUDING_STRING_BITS
#error this file is not meant to be included separately, include "vif/utilty/string.hpp" instead
#endif

namespace vif {
    inline uint_t find(const std::string& ts, const std::string& pattern) {
        auto p = ts.find(pattern);
        if (p != ts.npos) {
            return p;
        } else {
            return npos;
        }
    }

    inline std::string replace(std::string s, const std::string& pattern, const std::string& rep) {
        auto p = s.find(pattern);
        while (p != s.npos) {
            s.replace(p, pattern.size(), rep);
            p = s.find(pattern, p+rep.size());
        }

        return s;
    }

    inline std::string trim(std::string s, const std::string& chars = " \t\n\r") {
        std::size_t spos = s.find_first_of(chars);
        if (spos == 0) {
            std::size_t epos = s.find_first_not_of(chars);
            if (epos == s.npos) return "";
            s = s.substr(epos);
        }

        spos = s.find_last_of(chars);
        if (spos == s.size()-1) {
            std::size_t epos = s.find_last_not_of(chars);
            s = s.erase(epos+1, s.size() - epos+1);
        }

        return s;
    }

    inline bool begins_with(const std::string& s, const std::string& pattern) {
        if (s.size() < pattern.size()) return false;
        for (uint_t i = 0; i < pattern.size(); ++i) {
            if (s[i] != pattern[i]) return false;
        }

        return true;
    }

    inline bool ends_with(const std::string& s, const std::string& pattern) {
        if (s.size() < pattern.size()) return false;
        for (uint_t i = 1; i <= pattern.size(); ++i) {
            if (s[s.size()-i] != pattern[pattern.size()-i]) return false;
        }

        return true;
    }

    inline std::string erase_begin(std::string s, uint_t n) {
        if (n >= s.size()) {
            s.clear();
        } else {
            s.erase(0, n);
        }

        return s;
    }

    inline std::string erase_end(std::string s, uint_t n) {
        if (n >= s.size()) {
            s.clear();
        } else {
            s.erase(s.size()-n, n);
        }

        return s;
    }

    inline std::string erase_begin(std::string s, const std::string& pattern) {
        vif_check(begins_with(s, pattern), "unexpected string content: '"+s+"', "
            "should start with '"+pattern+"'");
        s.erase(0, pattern.size());
        return s;
    }

    inline std::string erase_end(std::string s, const std::string& pattern) {
        vif_check(ends_with(s, pattern), "unexpected string content: '"+s+"', "
            "should end with '"+pattern+"'");
        s.erase(s.size()-pattern.size(), pattern.size());
        return s;
    }

    template<typename T>
    std::string get_word(const std::string& str, int n, const T& chars) {
        if (n < 0) {
            n = -n;

            uint_t p1 = str.npos, p2 = str.npos;

            int nw = 0;
            while (n > nw) {
                p2 = str.find_last_not_of(chars, p1);
                vif_check(p2 != str.npos,
                    "not enough words (found ", nw, ", expected ", n, ")");

                p1 = str.find_last_of(chars, p2);
                vif_check(p1 != str.npos || nw == n-1,
                    "not enough words (found ", nw+1, ", expected ", n, ")");

                ++nw;
            }

            p1 = (p1 == str.npos ? 0 : p1 + 1);
            return str.substr(p1, p2-p1+1);
        } else {
            uint_t p1 = 0, p2 = 0;

            int nw = 0;
            while (n >= nw) {
                p1 = str.find_first_not_of(chars, p2);
                vif_check(p1 != str.npos,
                    "not enough words (found ", nw, ", expected ", n+1, ")");

                p2 = str.find_first_of(chars, p1);
                vif_check(p2 != str.npos || nw == n,
                    "not enough words (found ", nw+1, ", expected ", n+1, ")");

                ++nw;
            }

            p2 = (p2 == str.npos ? str.npos : p2-1);
            return str.substr(p1, p2-p1+1);
        }
    }

    template<typename T>
    std::string replace_block(const std::string& os, const std::string& begin,
        const std::string& end, T&& convert) {

        std::string s;

        uint_t p = os.find(begin);
        uint_t p0 = 0;
        while (p != npos) {
            s += os.substr(p0, p-p0);

            p += begin.size();

            p0 = os.find(end, p);
            vif_check(p0 != npos, "ill formed "+begin+"..."+end+" block");

            s += convert(os.substr(p, p0-p));
            p0 += end.size();

            p = os.find(begin, p0);
        }

        s += os.substr(p0);

        return s;
    }

    template<typename T>
    std::string replace_blocks(std::string s, const std::string& begin,
        const std::string& sep, const std::string& end, T&& convert) {

        std::string os = s;
        s.clear();

        uint_t p = os.find(begin);
        uint_t p0 = 0;
        while (p != npos) {
            s += os.substr(p0, p-p0);

            p += begin.size();

            vec1s b;

            p0 = os.find(sep, p);
            while (p0 != npos) {
                b.push_back(os.substr(p, p0-p));
                p0 += sep.size();
                p = p0;

                p0 = os.find(sep, p);
            }

            p0 = os.find(end, p);
            vif_check(p0 != npos, "ill formed "+begin+"..."+end+" command");

            b.push_back(os.substr(p, p0-p));
            s += convert(std::move(b));

            p0 += end.size();

            p = os.find(begin, p0);
        }

        s += os.substr(p0);

        return s;
    }

    VIF_VECTORIZE(find)
    VIF_VECTORIZE(replace)
    VIF_VECTORIZE(trim)
    VIF_VECTORIZE(erase_begin)
    VIF_VECTORIZE(erase_end)
    VIF_VECTORIZE(begins_with)
    VIF_VECTORIZE(ends_with)
}
