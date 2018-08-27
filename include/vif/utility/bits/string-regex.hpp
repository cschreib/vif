#ifndef VIF_INCLUDING_STRING_BITS
#error this file is not meant to be included separately, include "vif/utilty/string.hpp" instead
#endif

#include <regex.h>

namespace vif {
    namespace impl {
    namespace regex {
        inline std::string regex_get_error_(int status) {
            switch (status) {
            case REG_NOMATCH :  return "no match";
            case REG_BADPAT :   return "invalid regular expression";
            case REG_ECOLLATE : return "invalid collating element referenced";
            case REG_ECTYPE :   return "invalid character class type referenced";
            case REG_EESCAPE :  return "trailing \\ in pattern";
            case REG_ESUBREG :  return "number in \\digit invalid or in error";
            case REG_EBRACK :   return "[ ] imbalance";
            case REG_EPAREN :   return "\\( \\) or ( ) imbalance";
            case REG_EBRACE :   return "\\{ \\} imbalance";
            case REG_BADBR :    return "content of \\{ \\} invalid: not a number, number too large, "
                "more than two numbers, first larger than second";
            case REG_ERANGE :   return "invalid endpoint in range expression";
            case REG_ESPACE :   return "out of memory";
            case REG_BADRPT :   return "?, * or + not preceded by valid regular expression";
            case REG_ENOSYS :   return "the implementation does not support the function";
            default :           return "unknown error";
            }
        }

        inline void build_regex_(const std::string& regex, regex_t& re, int flags) {
            int status = regcomp(&re, regex.c_str(), flags);
            vif_check(status == 0, "parsing regex '", regex, "': ", regex_get_error_(status));
        }

        inline bool regex_match_(const std::string& ts, regex_t& re) {
            return regexec(&re, ts.c_str(), std::size_t(0), nullptr, 0) == 0;
        }

        inline uint_t regex_find_nmatch_(const std::string& regex) {
            uint_t nmatch = 0;
            uint_t p = regex.find_first_of('(');
            while (p != npos) {
                if (p == 0 || regex[p-1] != '\\') {
                    ++nmatch;
                }

                p = regex.find_first_of('(', p+1);
            }

            return nmatch;
        }
    }
    }

    inline bool regex_match(const std::string& ts, const std::string& regex) {
        regex_t re;
        impl::regex::build_regex_(regex, re, REG_EXTENDED | REG_NOSUB);
        bool ret = impl::regex::regex_match_(ts, re);
        regfree(&re);
        return ret;
    }

    template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
        std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
    vec<Dim,bool> regex_match(const vec<Dim,Type>& v, const std::string& regex) {
        vec<Dim,bool> r(v.dims);
        regex_t re;
        impl::regex::build_regex_(regex, re, REG_EXTENDED | REG_NOSUB);
        for (uint_t i : range(v)) {
            r.safe[i] = impl::regex::regex_match_(v.safe[i], re);
        }
        regfree(&re);
        return r;
    }

    template<std::size_t Dim, typename Type, typename enable = typename std::enable_if<
        std::is_same<typename std::remove_pointer<Type>::type, std::string>::value>::type>
    vec<Dim,bool> regex_match(const std::string& v, const vec<Dim,Type>& regex) {
        vec<Dim,bool> r(regex.dims);
        for (uint_t i : range(regex)) {
            regex_t re;
            impl::regex::build_regex_(regex.safe[i], re, REG_EXTENDED | REG_NOSUB);
            r.safe[i] = impl::regex::regex_match_(v, re);
            regfree(&re);
        }
        return r;
    }

    inline vec2s regex_extract(const std::string& ts, const std::string& regex) {
        vec2s ret;

        regex_t re;
        impl::regex::build_regex_(regex, re, REG_EXTENDED);

        uint_t nmatch = impl::regex::regex_find_nmatch_(regex);
        if (nmatch == 0) return ret;

        std::vector<regmatch_t> m(nmatch+1);
        uint_t offset = 0;
        int status = regexec(&re, ts.c_str(), nmatch+1, m.data(), 0);

        while (status == 0) {
            vec1s tret;
            tret.reserve(nmatch);
            for (uint_t i : range(nmatch)) {
                if (m[i+1].rm_eo > m[i+1].rm_so) {
                    tret.push_back(ts.substr(offset+m[i+1].rm_so, m[i+1].rm_eo - m[i+1].rm_so));
                } else {
                    tret.push_back("");
                }
            }

            append<0>(ret, reform(std::move(tret), 1, nmatch));

            offset += m[0].rm_eo;
            status = regexec(&re, ts.c_str() + offset, nmatch+1, m.data(), 0);
        }

        return ret;
    }

    template<typename F>
    std::string regex_replace(const std::string& ts, const std::string& regex, F&& func) {
        std::string s;

        regex_t re;
        impl::regex::build_regex_(regex, re, REG_EXTENDED);

        uint_t nmatch = impl::regex::regex_find_nmatch_(regex);

        std::vector<regmatch_t> m(nmatch+1);
        uint_t offset = 0;
        int status = regexec(&re, ts.c_str(), nmatch+1, m.data(), 0);

        while (status == 0) {
            s += ts.substr(offset, m[0].rm_so);

            vec1s ext;
            ext.reserve(nmatch);
            for (uint_t i : range(nmatch)) {
                ext.push_back(ts.substr(offset+m[i+1].rm_so, m[i+1].rm_eo - m[i+1].rm_so));
            }

            s += func(std::move(ext));

            offset += m[0].rm_eo;
            status = regexec(&re, ts.c_str() + offset, nmatch+1, m.data(), 0);
        }

        s += ts.substr(offset);

        return s;
    }
}
