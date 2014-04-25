#ifndef STRING_HPP
#define STRING_HPP

#include "phypp/vec.hpp"
#include "phypp/reflex.hpp"
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <regex.h>

namespace std {
    template<typename O, typename T, std::size_t N>
    O& operator << (O& o, const std::array<T,N>& v) {
        for (std::size_t i = 0; i < N; ++i) {
            if (i != 0) o << ", ";
            o << v[i];
        }

        return o;
    }
}

template<typename T>
struct named_t {
    T&          obj;
    std::string name;
};

template<typename T>
named_t<T> name(T& obj, const std::string& str) {
    return {obj, str};
}

template<typename T>
std::string strn(const T& t) {
    std::ostringstream ss;
    ss << reflex::wrap(t);
    return ss.str();
}

std::string strn(const double& t) {
    std::ostringstream ss;
    ss.precision(12);
    ss << t;
    return ss.str();
}

template<typename T>
std::string strn(const T& t, std::size_t n, char fill = '0') {
    if (n <= 1) {
        return strn(t);
    }
    if (t == 0) {
        return std::string(n-1, fill)+'0';
    } else {
        std::ostringstream ss;
        ss << t;
        std::size_t nz = n-1 - floor(log10(t));
        if (nz > 0 && nz < 6) {
            return std::string(nz, fill) + ss.str();
        } else {
            return ss.str();
        }
    }
}

template<std::size_t Dim, typename Type>
vec_t<Dim,std::string> strna(const vec_t<Dim,Type>& v) {
    vec_t<Dim,std::string> s = strarr(v.dims);
    for (std::size_t i = 0; i < v.data.size(); ++i) {
        s.data[i] = strn(dref<Type>(v.data[i]));
    }

    return s;
}

template<std::size_t Dim, typename Type>
vec_t<Dim,std::string> strna(const vec_t<Dim,Type>& v, std::size_t n, char fill = '0') {
    vec_t<Dim,std::string> s = strarr(v.dims);
    for (std::size_t i = 0; i < v.data.size(); ++i) {
        s.data[i] = strn(dref<Type>(v.data[i]), n, fill);
    }

    return s;
}

template<typename T>
std::string strn_sci(const T& t) {
    std::ostringstream ss;
    ss << std::scientific << reflex::wrap(t);
    return ss.str();
}


template<std::size_t Dim, typename Type>
vec_t<Dim,std::string> strna_sci(const vec_t<Dim,Type>& v) {
    vec_t<Dim,std::string> s = strarr(v.dims);
    for (std::size_t i = 0; i < v.data.size(); ++i) {
        s.data[i] = strn_sci(dref<Type>(v.data[i]));
    }

    return s;
}

template<typename T>
bool from_string(const std::string& s, T& t) {
    std::istringstream ss(s);
    return (ss >> t);
}

template<std::size_t Dim, typename T>
vec_t<Dim,bool> from_string(const vec_t<Dim,std::string>& s, vec_t<Dim,T>& t) {
    vec_t<Dim,bool> res(s.dims);
    t.resize(s.dims);
    for (uint_t i : range(s)) {
        res[i] = from_string(s[i], t[i]);
    }

    return res;
}

template<std::size_t Dim, typename T>
vec_t<Dim,bool> from_string(const vec_t<Dim,std::string*>& s, vec_t<Dim,T>& t) {
    vec_t<Dim,bool> res(s.dims);
    t.resize(s.dims);
    for (uint_t i : range(s)) {
        res[i] = from_string(s[i], t[i]);
    }

    return res;
}

template<typename T>
T nstr(const std::string& s) {
    std::istringstream ss(s);
    T t;
    assert(ss >> t);
    return t;
}

std::string trim(std::string s, const std::string& chars = " \t") {
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

std::string toupper(std::string s) {
    for (auto& c : s) {
        c = ::toupper(c);
    }
    return s;
}

std::string tolower(std::string s) {
    for (auto& c : s) {
        c = ::tolower(c);
    }
    return s;
}

std::string replace(std::string s, const std::string& pattern, const std::string& rep) {
    auto p = s.find(pattern);
    while (p != s.npos) {
        s.replace(p, pattern.size(), rep);
        p = s.find(pattern, p+rep.size());
    }

    return s;
}

bool empty(const std::string& s) {
    return s.empty();
}

uint_t distance(const std::string& t, const std::string& u) {
    uint_t n = std::min(t.size(), u.size());
    uint_t d = abs(t.size() - u.size());
    for (uint_t i = 0; i < n; ++i) {
        if (t[i] != u[i]) ++d;
    }

    return d;
}

uint_t find(const std::string& ts, const std::string& pattern) {
    auto p = ts.find(pattern);
    if (p != ts.npos) {
        return p;
    } else {
        return npos;
    }
}

void build_regex_(const std::string& regex, regex_t& re) {
    int status = regcomp(&re, regex.c_str(), REG_EXTENDED | REG_NOSUB);
    if (status != 0) {
        error("match: parsing regex '"+regex+"'");
        switch (status) {
        case REG_NOMATCH :  print("    regexec() failed to match. "); break;
        case REG_BADPAT :   print("    Invalid regular expression. "); break;
        case REG_ECOLLATE : print("    Invalid collating element referenced. "); break;
        case REG_ECTYPE :   print("    Invalid character class type referenced. "); break;
        case REG_EESCAPE :  print("    Trailing \\ in pattern. "); break;
        case REG_ESUBREG :  print("    Number in \\digit invalid or in error. "); break;
        case REG_EBRACK :   print("    [ ] imbalance. "); break;
        case REG_EPAREN :   print("    \\( \\) or ( ) imbalance. "); break;
        case REG_EBRACE :   print("    \\{ \\} imbalance. "); break;
        case REG_BADBR :    print("    Content of \\{ \\} invalid: not a number, number too large, "
            "more than two numbers, first larger than second. "); break;
        case REG_ERANGE :   print("    Invalid endpoint in range expression. "); break;
        case REG_ESPACE :   print("    Out of memory. "); break;
        case REG_BADRPT :   print("    ?, * or + not preceded by valid regular expression. "); break;
        case REG_ENOSYS :   print("    The implementation does not support the function. "); break;
        default :           print("    Unknown error."); break;
        }
        assert(false);
    }
}

bool match_(const std::string& ts, const std::string& regex, regex_t& re) {
    int status = regexec(&re, ts.c_str(), (size_t)0, nullptr, 0);

    if (status == 0) return true;
    if (status == REG_NOMATCH) return false;

    error("match: parsing regex '"+regex+"'");
    switch (status) {
    case REG_NOMATCH :  return false;
    case REG_BADPAT :   print("    Invalid regular expression. "); break;
    case REG_ECOLLATE : print("    Invalid collating element referenced. "); break;
    case REG_ECTYPE :   print("    Invalid character class type referenced. "); break;
    case REG_EESCAPE :  print("    Trailing \\ in pattern. "); break;
    case REG_ESUBREG :  print("    Number in \\digit invalid or in error. "); break;
    case REG_EBRACK :   print("    [ ] imbalance. "); break;
    case REG_EPAREN :   print("    \\( \\) or ( ) imbalance. "); break;
    case REG_EBRACE :   print("    \\{ \\} imbalance. "); break;
    case REG_BADBR :    print("    Content of \\{ \\} invalid: not a number, number too large, "
        "more than two numbers, first larger than second. "); break;
    case REG_ERANGE :   print("    Invalid endpoint in range expression. "); break;
    case REG_ESPACE :   print("    Out of memory. "); break;
    case REG_BADRPT :   print("    ?, * or + not preceded by valid regular expression. "); break;
    case REG_ENOSYS :   print("    The implementation does not support the function. "); break;
    default :           print("    Unknown error."); break;
    }
    assert(false);
    return false;
}

bool match(const std::string& ts, const std::string& regex) {
    regex_t re;
    build_regex_(regex, re);
    bool ret = match_(ts, regex, re);
    regfree(&re);
    return ret;
}

template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
    std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
vec_t<Dim,bool> match(const vec_t<Dim,Type>& v, const std::string& regex) {
    regex_t re;
    build_regex_(regex, re);
    vec_t<Dim,bool> r = boolarr(v.dims);
    for (uint_t i = 0; i < v.size(); ++i) {
        r.data[i] = match_(dref<Type>(v.data[i]), regex, re);
    }
    regfree(&re);
    return r;
}


uint_t length(const std::string& s) {
    return s.size();
}

template<typename T, typename enable = typename std::enable_if<std::is_arithmetic<T>::value>::type>
uint_t length(T c) {
    return 1u;
}

std::string align_left(std::string s, uint_t width, char fill = ' ') {
    if (s.size() < width) {
        s += std::string(width-s.size(), fill);
    }

    return s;
}

std::string align_right(std::string s, uint_t width, char fill = ' ') {
    if (s.size() < width) {
        s.insert(0, std::string(width-s.size(), fill));
    }

    return s;
}

std::string align_center(std::string s, uint_t width, char fill = ' ') {
    if (s.size() < width) {
        uint_t n1 = (width-s.size())/2, n2 = width-s.size() - n1;
        s.insert(0, std::string(n1, fill));
        s += std::string(n2, fill);
    }

    return s;
}

bool start_with(const std::string& s, const std::string& pattern) {
    if (s.size() < pattern.size()) return false;
    for (uint_t i = 0; i < pattern.size(); ++i) {
        if (s[i] != pattern[i]) return false;
    }

    return true;
}

bool end_with(const std::string& s, const std::string& pattern) {
    if (s.size() < pattern.size()) return false;
    for (uint_t i = 1; i <= pattern.size(); ++i) {
        if (s[s.size()-i] != pattern[pattern.size()-i]) return false;
    }

    return true;
}

std::string erase_begin(std::string s, uint_t n) {
    if (n >= s.size()) {
        s.clear();
    } else {
        s.erase(0, n);
    }

    return s;
}

std::string erase_end(std::string s, uint_t n) {
    if (n >= s.size()) {
        s.clear();
    } else {
        s.erase(s.size()-n, n);
    }

    return s;
}

std::string erase_begin(std::string s, const std::string& pattern) {
    phypp_check(start_with(s, pattern), "unexpected string content: '"+s+"', "
        "should start with '"+pattern+"'");
    s.erase(0, pattern.size());
    return s;
}

std::string erase_end(std::string s, const std::string& pattern) {
    phypp_check(end_with(s, pattern), "unexpected string content: '"+s+"', "
        "should end with '"+pattern+"'");
    s.erase(s.size()-pattern.size(), pattern.size());
    return s;
}

std::string keep_start(std::string s, uint_t n = 1) {
    if (s.size() > n) {
        s.erase(n);
    }

    return s;
}

std::string keep_end(std::string s, uint_t n = 1) {
    if (s.size() > n) {
        s.erase(0, s.size()-n);
    }

    return s;
}

std::string remove_extension(std::string s) {
    auto p = s.find_last_of('.');
    if (p == s.npos) return s;
    return s.substr(0u, p);
}

#define VECTORIZE(name) \
    template<std::size_t Dim, typename Type, typename ... Args, \
        typename enable = typename std::enable_if< \
            std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type> \
    auto name(const vec_t<Dim,Type>& v, const Args& ... args) -> \
        vec_t<Dim,decltype(name(v[0], args...))> { \
        using ntype = decltype(name(v[0], args...)); \
        vec_t<Dim,ntype> r = arr<ntype>(v.dims); \
        for (std::size_t i = 0; i < v.data.size(); ++i) { \
            r.data[i] = name(dref<Type>(v.data[i]), args...); \
        } \
        return r; \
    }

VECTORIZE(trim)
VECTORIZE(toupper)
VECTORIZE(tolower)
VECTORIZE(erase_begin)
VECTORIZE(erase_end)
VECTORIZE(empty)
VECTORIZE(distance)
VECTORIZE(find)
VECTORIZE(replace)
VECTORIZE(length)
VECTORIZE(align_left)
VECTORIZE(align_right)
VECTORIZE(align_center)
VECTORIZE(start_with)
VECTORIZE(end_with)
VECTORIZE(keep_start)
VECTORIZE(keep_end)

#undef VECTORIZE

template<typename T>
vec1s split(const std::string& ts, const T& pattern) {
    vec1s ret;
    std::size_t p = 0, op = 0;
    while ((p = ts.find(pattern, op)) != ts.npos) {
        ret.data.push_back(ts.substr(op, p - op));
        op = p + length(pattern);
    }

    ret.data.push_back(ts.substr(op));
    ret.dims[0] = ret.data.size();

    return ret;
}

vec1s cut(const std::string& ts, uint_t size) {
    uint_t ncut = floor(ts.size()/float(size));
    vec1s res(ncut);
    for (uint_t i = 0; i < ncut; ++i) {
        res[i] = ts.substr(i*size, std::min(size, ts.size() - i*size));
    }

    return res;
}

vec1s wrap(const std::string& ts, uint_t width, const std::string& indent = "", bool ellipse = false) {
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

namespace format {
    void header(const std::string& msg) {
        vec1s w = wrap("  "+msg, 80, "  ");
        for (auto& s : w) {
            print(s);
        }
    }

    void paragraph(const std::string& msg) {
        vec1s w = wrap("  "+msg, 80, "  ");
        for (auto& s : w) {
            print(s);
        }
        print("");
    }

    void bullet(const std::string& name, const std::string& desc) {
        std::string header = "    "+name+": ";
        vec1s w = wrap(header+desc, 80, std::string(header.size(), ' '));
        for (auto& s : w) {
            print(s);
        }
    }
}

template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
    std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
std::string collapse(const vec_t<Dim,Type>& v) {
    std::string r;
    for (auto& s : v) {
        r += s;
    }

    return r;
}

template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
    std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
std::string collapse(const vec_t<Dim,Type>& v, const std::string& sep) {
    std::string r;
    bool first = true;
    for (auto& s : v) {
        if (!first) r += sep;
        r += s;
        first = false;
    }

    return r;
}

std::string system_var(const std::string& name, const std::string& def = "") {
    char* v = getenv(name.c_str());
    if (!v) return def;
    return v;
}

#endif
