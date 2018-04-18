#ifndef PHYPP_INCLUDING_STRING_BITS
#error this file is not meant to be included separately, include "phypp/utilty/string.hpp" instead
#endif

namespace phypp {
    template<typename T>
    vec1s split(const std::string& ts, const T& pattern) {
        vec1s ret;
        std::size_t p = 0, op = 0;
        while ((p = ts.find(pattern, op)) != ts.npos) {
            ret.data.push_back(ts.substr(op, p - op));
            op = p + length(pattern);
        }

        ret.data.push_back(ts.substr(op));
        ret.dims[0] = ret.size();

        return ret;
    }

    template<typename T>
    vec1s split_any_of(const std::string& ts, const T& chars) {
        vec1s ret;
        std::size_t op = ts.find_first_not_of(chars);
        std::size_t p = op;
        while ((p = ts.find_first_of(chars, op)) != ts.npos) {
            ret.data.push_back(ts.substr(op, p - op));
            op = ts.find_first_not_of(chars, p);
        }

        if (op != ts.npos) {
            ret.data.push_back(ts.substr(op));
        }

        ret.dims[0] = ret.size();

        return ret;
    }

    inline vec1s cut(const std::string& ts, uint_t size) {
        uint_t ncut = floor(ts.size()/float(size));
        vec1s res(ncut);
        for (uint_t i = 0; i < ncut; ++i) {
            res.safe[i] = ts.substr(i*size, std::min(size, ts.size() - i*size));
        }

        return res;
    }

    inline vec1s wrap(const std::string& ts, uint_t width, const std::string& indent = "", bool ellipse = false) {
        vec1s ret;
        std::string s = ts;
        uint_t twidth = width;
        std::string header = "";
        while (s.size() > twidth) {
            uint_t i = twidth;
            while (i != npos && s[i] != ' ') --i;
            if (i == npos) {
                if (ellipse) {
                    ret.push_back(header+s.substr(0, twidth-3)+"...");
                    i = twidth+1;
                    while (i < s.size() && s[i] != ' ') ++i;
                    s.erase(0, i);
                    s = trim(s);
                } else {
                    i = twidth+1;
                    while (i < s.size() && s[i] != ' ') ++i;
                    ret.push_back(header+s.substr(0, i));
                    s.erase(0, i);
                    s = trim(s);
                }
            } else {
                ret.push_back(header+s.substr(0, i));
                s.erase(0, i);
                s = trim(s);
            }

            twidth = width - indent.size();
            header = indent;
        }

        if (!s.empty()) {
            ret.push_back(header+s);
        }

        return ret;
    }

    template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
        std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
    std::string collapse(const vec<Dim,Type>& v) {
        std::string r;
        for (auto& s : v) {
            r += s;
        }

        return r;
    }

    template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
        std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
    std::string collapse(const vec<Dim,Type>& v, const std::string& sep) {
        std::string r;
        bool first = true;
        for (auto& s : v) {
            if (!first) r += sep;
            r += s;
            first = false;
        }

        return r;
    }
}
