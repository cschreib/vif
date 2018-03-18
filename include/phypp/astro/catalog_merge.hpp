#ifndef ASTRO_CATALOG_MERGE_HPP
#define ASTRO_CATALOG_MERGE_HPP

#include <list>
#include "phypp/astro/astro.hpp"
#include "phypp/astro/qxmatch.hpp"

namespace phypp {
namespace astro {
    struct comment_pool_t {
        std::map<std::string, std::string> var_pool;
        std::vector<std::string> com_list;
    };

    comment_pool_t create_comment_pool() {
        return comment_pool_t();
    }

    void add_comment(comment_pool_t& pool, const std::string& com) {
        pool.com_list.push_back(com);
    }

    void add_comment(comment_pool_t& pool, const std::string& var, const std::string& com) {
        pool.var_pool.insert(std::make_pair(var, com));
    }

    std::string build_comments(const comment_pool_t& pool, uint_t width = 80) {
        std::string str;

        auto add_cmt = [&str, width] (const std::string& c, const std::string& tab = "") {
            vec1s spl = wrap(c, width, tab);
            for (uint_t i = 0; i < spl.size(); ++i) {
                str += spl[i] + '\n';
            }
        };

        for (auto& s : pool.com_list) {
            add_cmt(s);
        }

        vec1s tree;
        for (auto& p : pool.var_pool) {
            vec1s ltree = trim(split(p.first, "."));
            uint_t nmax = std::min(tree.size(), ltree.size());
            uint_t icom;
            std::string tab = "";
            for (icom = 0; icom < nmax; ++icom) {
                if (tree[icom] != ltree[icom]) break;
                if (icom % 2 == 0) {
                    tab += "| ";
                } else {
                    tab += ": ";
                }
            }

            uint_t i;
            for (i = icom; i < ltree.size()-1; ++i) {
                str += tab+ltree[i]+'\n';
                if (i % 2 == 0) {
                    tab += "| ";
                } else {
                    tab += ": ";
                }
            }

            std::string header = ltree[i]+": ";
            add_cmt(tab+header+p.second, tab+std::string(header.size(), ' '));

            tree = ltree;
        }

        return str;
    }

    struct catalog_pool;

    struct catalog_t {
        catalog_pool& pool;

        vec1u       idi;
        vec1u       idm;
        vec1u       idb;
        vec2d       d;

        std::string name;
        vec1s       sources;
        vec1s       files;
        std::string comment;
        std::string ref;

        template<std::size_t Dim, typename T, typename U, typename V,
            typename enable = typename std::enable_if<Dim != 1>::type>
        void merge(vec<Dim,T>& in, const vec<Dim,U>& out, const V& def);

        template<typename T, typename U, typename V>
        void merge(vec<1,T>& in, const U& out, const V& def);

        template<typename T, typename U, typename V>
        void merge(vec<1,T>& in, const vec<1,U>& out, const V& def);

        template<std::size_t Dim, typename T, typename U, typename V,
            typename enable = typename std::enable_if<Dim != 1>::type>
        void merge(vec<Dim,T>& in, const vec<Dim,U>& out, const V& def, const std::string& com);

        template<typename T, typename U, typename V>
        void merge(vec<1,T>& in, const U& out, const V& def, const std::string& com);

        template<typename T, typename U, typename V>
        void merge(T& in, const U& out, const V& def);

        template<typename T, typename U, typename V>
        void merge(T& in, const U& out, const V& def, const std::string& com);

        template<typename T, typename U>
        void assign(T& in, const U& out, const std::string& com);

        void merge_flux(const vec2f& flux, const vec2f& err, const vec1s& bands, const vec1s& notes);
    };

    struct catalog_pool {
        comment_pool_t coms;
        uint_t ngal;
        std::list<catalog_t> pool;

        vec1u origin;
        vec1d ra, dec;
        vec2f flux, flux_err;
        vec1s bands, notes;

        bool xmatch;
        std::string xmatch_file;

        catalog_pool(bool xmatch_, const std::string& xmatch_file_) :
            xmatch(xmatch_), xmatch_file(xmatch_file_) {
            file::mkdir(file::get_directory(xmatch_file));
        }

        catalog_t& add_catalog(const std::string& name, const vec1s& sources, const vec1s& files,
            const std::string& comment = "") {

            std::string ref = "["+to_string(pool.size()+1)+"]";
            pool.push_back({*this, uindgen(ngal), uindgen(ngal), {}, {}, name, sources, files, comment, ref});
            return pool.back();
        }

    private :
        catalog_t& add_catalog_(const vec1d& cra, const vec1d& cdec, bool no_new, vec1u sel,
            const std::string& name, const vec1s& sources, const vec1s& files,
            const std::string& comment = "") {

            phypp_check(cra.size() == cdec.size(), "need ra.size() == dec.size()");

            if (sel.empty()) {
                sel = uindgen(cra.size());
            }

            if (pool.empty()) {
                ngal = sel.size();
                origin = replicate(1, ngal);
                ra = cra[sel];
                dec = cdec[sel];
                vec1u idm = uindgen(ngal);
                std::string ref = "[1]";
                pool.push_back({*this, sel, idm, {}, {}, name, sources, files, comment, ref});
                return pool.back();
            } else {
                print("cross-matching "+name+"...");

                std::string hashed = hash(name, cra.size(), sel.size(), ngal, sources, files);
                std::string file = xmatch_file+"_"+hashed+".fits";

                bool rematch = false;
                if (!xmatch && file::exists(file)) {
                    struct {
                        vec1u sel;
                        uint_t ngal;
                    } tmp;

                    fits::read_table(file, ftable(tmp.sel, tmp.ngal));
                    if (tmp.sel.size() != sel.size() || tmp.ngal != ngal) {
                        warning("incompatible cross-match data, re-matching");
                        rematch = true;
                    }
                } else {
                    rematch = true;
                }

                qxmatch_res xm;

                if (rematch) {
                    qxmatch_params p; p.nth = 2; p.thread = 4; p.verbose = true;
                    xm = qxmatch(cra[sel], cdec[sel], ra, dec, p);

                    fits::write_table(file, ftable(
                        sel, ngal, xm.id, xm.d, xm.rid, xm.rd
                    ));
                } else {
                    print("reading data from "+file);
                    fits::read_table(file, ftable(xm.id, xm.d, xm.rid, xm.rd));
                }

                const uint_t n = sel.size();
                vec1u idm; idm.reserve(n);
                vec1u idn;
                vec1u tsel;
                if (no_new) {
                    tsel.reserve(n);
                } else {
                    tsel = uindgen(sel.size());
                }

                for (uint_t i = 0; i < n; ++i) {
                    if (xm.rid[xm.id(0,i)] == i) {
                        idm.push_back(xm.id(0,i));
                        if (no_new) {
                            tsel.push_back(i);
                        }
                    } else {
                        if (!no_new) {
                            idn.push_back(sel[i]);
                            idm.push_back(ra.size());
                            ra.push_back(cra[sel[i]]);
                            dec.push_back(cdec[sel[i]]);
                            origin.push_back(pool.size()+1);
                        }
                    }
                }

                vec2d td = replicate(dnan, xm.d.dims[0], cra.size());
                td(_,sel[tsel]) = xm.d(_,tsel);

                sel = sel[tsel];

                if (no_new) {
                    print("> ", idm.size(), " matched sources, ", n - idm.size(), " lost");
                } else {
                    print("> ", n - idn.size(), " matched sources, ", idn.size(), " new sources");
                    ngal += idn.size();
                }

                std::string ref = "["+to_string(pool.size()+1)+"]";
                pool.push_back({*this, sel, idm, {}, std::move(td), name, sources, files, comment, ref});
                return pool.back();
            }
        }

    public :
        catalog_t& add_catalog(const vec1d& cra, const vec1d& cdec, vec1u sel,
            const std::string& name, const vec1s& sources, const vec1s& files,
            const std::string& comment = "") {
            return add_catalog_(cra, cdec, false, sel, name, sources, files, comment);
        }

        catalog_t& match_catalog(const vec1d& cra, const vec1d& cdec, vec1u sel,
            const std::string& name, const vec1s& sources, const vec1s& files,
            const std::string& comment = "") {
            return add_catalog_(cra, cdec, true, sel, name, sources, files, comment);
        }

        std::string build_catalog_list(uint_t width = 80) {
            std::string str;

            for (auto& c : pool) {
                str += c.ref+": "+c.name+"\n";
                if (!c.sources.empty()) {
                    if (c.sources.size() == 1) {
                        str += "  source: "+c.sources[0]+"\n";
                    } else {
                        std::string header = "  sources: ";
                        str += header+c.sources[0]+"\n";
                        for (uint_t i = 1; i < c.sources.size(); ++i) {
                            str += std::string(header.size(), ' ')+c.sources[i]+"\n";
                        }
                    }
                }
                if (!c.files.empty()) {
                    if (c.files.size() == 1) {
                        str += "  file: "+c.files[0]+"\n";
                    } else {
                        std::string header = "  files: ";
                        str += header+c.files[0]+"\n";
                        for (uint_t i = 1; i < c.files.size(); ++i) {
                            str += std::string(header.size(), ' ')+c.files[i]+"\n";
                        }
                    }
                }
                if (!c.comment.empty()) {
                    vec1s spl = wrap("  comment: "+c.comment, width, std::string(4, ' '));
                    for (uint_t i = 0; i < spl.size(); ++i) {
                        str += spl[i] + '\n';
                    }
                }
                str += "\n";
            }

            return str;
        }
    };

    template<std::size_t Dim, typename T, typename U, typename V, typename enable>
    void catalog_t::merge(vec<Dim,T>& in, const vec<Dim,U>& out, const V& def) {
        phypp_check(in.empty(), name+": merging twice into the same variable");
        in.dims = out.dims;
        in.dims[0] = pool.ngal;
        in.resize();
        in[_] = def;
        in(idm, repeat<Dim-1>(_)) = out(idi, repeat<Dim-1>(_));
    }

    template<typename T, typename U, typename V>
    void catalog_t::merge(vec<1,T>& in, const U& out, const V& def) {
        phypp_check(in.empty(), name+": merging twice into the same variable");
        in = replicate(def, pool.ngal);
        in[idm] = out;
    }

    template<typename T, typename U, typename V>
    void catalog_t::merge(vec<1,T>& in, const vec<1,U>& out, const V& def) {
        phypp_check(in.empty(), name+": merging twice into the same variable");
        in = replicate(def, pool.ngal);
        in[idm] = out[idi];
    }

    template<std::size_t Dim, typename T, typename U, typename V, typename enable>
    void catalog_t::merge(vec<Dim,T>& in, const vec<Dim,U>& out, const V& def,
        const std::string& com) {
        merge(in, out, def);

        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
    }

    template<typename T, typename U, typename V>
    void catalog_t::merge(vec<1,T>& in, const U& out, const V& def, const std::string& com) {
        merge(in, out, def);

        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
    }
}

namespace impl {
    namespace astro_impl {
        template<typename M, typename V>
        struct do_catalog_merge_run {
            const reflex::member_t& m;
            const M& v;
            const V& def;
            catalog_t& cat;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    this->cat.merge(p, this->v, this->def);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename M, typename V>
        struct do_catalog_merge_run<reflex::struct_t<M>,V> {
            const reflex::member_t& m;
            const M& v;
            const V& def;
            catalog_t& cat;

            void operator () (reflex::member_t& n, M& p) {
                if (this->m.name == n.name) {
                    this->cat.merge(p, this->v, this->def);
                }
            }

            template<typename P>
            void operator () (reflex::member_t& n, reflex::struct_t<P> p) {
                this->cat.merge(p, this->v, this->def);
            }

            template<typename P>
            void operator () (reflex::member_t& n, P&& p) {
                phypp_check(this->m.name != n.name, "incompatible types in merging '",
                    this->m.full_name(), "' into '", n.full_name(), "'"
                );
            }
        };

        template<typename T, typename V>
        struct do_catalog_merge {
            reflex::struct_t<T> t;
            const V& def;
            catalog_t& cat;

            template<typename M>
            void operator () (const reflex::member_t& m, const M& v) {
                do_catalog_merge_run<M,V> run{m, v, this->def, this->cat};
                reflex::foreach_member(this->t, run);
            }
        };
    }
}

namespace astro {
    template<typename T, typename U, typename V>
    void catalog_t::merge(T& in, const U& out, const V& def) {
        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        impl::astro_impl::do_catalog_merge<T,V> run{reflex::wrap(in), def, *this};
        reflex::foreach_member(reflex::wrap(out), run);
    }

    template<typename T, typename U, typename V>
    void catalog_t::merge(T& in, const U& out, const V& def, const std::string& com) {
        merge(in, out, def);

        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
    }

    template<typename T, typename U>
    void catalog_t::assign(T& in, const U& out, const std::string& com) {
        in = out;

        #ifdef NO_REFLECTION
        static_assert(!std::is_same<T,T>::value, "this function requires reflection capabilities (NO_REFLECTION=0)");
        #endif

        add_comment(pool.coms, reflex::seek_name(in), com+" (from "+ref+")");
    }

    void catalog_t::merge_flux(const vec2f& flux, const vec2f& err, const vec1s& bands,
        const vec1s& notes) {

        phypp_check(flux.dims[1] == err.dims[1], name+": flux and error do not match");
        phypp_check(flux.dims[1] == bands.size(), name+": flux and band do not match");

        vec1s tnotes = notes;
        if (tnotes.empty()) {
            tnotes = strarr(bands.size());
        }

        vec1u idne = where(!empty(tnotes));
        tnotes[idne] = tnotes[idne]+" ";
        tnotes += ref;

        vec2f tflux, terr;
        merge(tflux, flux, fnan);
        merge(terr, err, fnan);

        idb = uindgen(bands.size()) + pool.bands.size();

        append(pool.bands, bands);
        append(pool.notes, tnotes);
        append<1>(pool.flux, tflux);
        append<1>(pool.flux_err, terr);
    }
}
}

#endif
