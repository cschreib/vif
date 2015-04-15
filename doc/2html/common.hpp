struct item_t {
    std::string key;
    std::string value;
};

using sorted_cache_t = std::map<std::string, std::string>;
using cache_t = std::vector<item_t>;

void insert_cache(sorted_cache_t& cache, std::string k, std::string v) {
    cache.insert(std::make_pair(k, v));
}

void insert_cache(cache_t& cache, std::string k, std::string v) {
    item_t i;
    i.key = k;
    i.value = v;
    cache.push_back(i);
}

template<typename T>
void read_cache(T& cache, std::string tfile) {
    std::ifstream in(tfile);
    std::string tline;
    while (std::getline(in, tline)) {
        tline = trim(tline);
        uint_t p = tline.find_first_of('=');
        std::string key = tline.substr(0, p);
        std::string value = replace(tline.substr(p+1), "__[\\n]__", "\n");
        insert_cache(cache, key, value);
    }
}

void write_cache(const sorted_cache_t& cache, std::string tfile) {
    std::ofstream tout(tfile);
    for (auto& kv : cache) {
        tout << kv.first << "=" << replace(kv.second, "\n", "__[\\n]__") << "\n";
    }
}

void write_cache(const cache_t& cache, std::string tfile) {
    std::ofstream tout(tfile);
    for (auto& kv : cache) {
        tout << kv.key << "=" << replace(kv.value, "\n", "__[\\n]__") << "\n";
    }
}

struct context_t {
    sorted_cache_t pygment_cache;
    sorted_cache_t latex_cache;
    sorted_cache_t section_index;
    sorted_cache_t footnotes;

    std::string base;

    std::string out_dir;
    std::string file_name;
    uint_t l = 0;
};

std::string pygmentize(context_t& con, std::string lang, std::string s) {
    std::string ret;
    std::string hsh = hash(s);
    auto iter = con.pygment_cache.find(hsh);
    if (iter == con.pygment_cache.end()) {
        std::ofstream out(temporary_dir+"tmp.cpp");
        out << s;
        out.close();

        spawn("pygmentize -f html -l "+lang+" -o "+temporary_dir+"tmp.html "+
                                              temporary_dir+"tmp.cpp");

        s = file::to_string(temporary_dir+"tmp.html");
        auto p = s.find_last_of('\n');
        if (trim(s.substr(p+1)).empty()) {
            s = s.substr(0, p);
        }

        con.pygment_cache.insert(std::make_pair(hsh, s));
    } else {
        s = iter->second;
    }

    return s;
}

std::string pygmentize_inline(context_t& con, std::string lang, std::string s) {
    s = pygmentize(con, lang, s);
    s = replace(s, "class=\"highlight\"", "class=\"highlight-inline\"");
    s = replace(s, "\n", "");
    return s;
}

std::string latex_math(context_t& con, std::string s) {
    static const std::string latex_header = file::to_string("header.tex");
    static const std::string latex_footer = file::to_string("footer.tex");

    std::string ret;
    std::string hsh = hash(s);
    auto iter = con.latex_cache.find(hsh);
    if (iter == con.latex_cache.end()) {
        static uint_t count = 0;
        std::ofstream out(temporary_dir+"tmp.tex");
        out << latex_header;
        out << "\\newcommand\\formula{$" << s << "$}\n";
        out << latex_footer;
        out.close();

        std::string out_file = "latex/"+con.base+"_"+strn(count)+".svg";
        file::remove(temporary_dir+"tmp.dvi");
        spawn("cd "+temporary_dir+"; latex -interaction batchmode tmp.tex");
        if (!file::exists(temporary_dir+"tmp.dvi")) {
            error("reading ", con.file_name);
            error("parsing $"+s+"$ on l.", con.l);
            return "$"+s+"$";
        }

        spawn("dvisvgm -v 0 -S -c 1.3 -n -o "+con.out_dir+out_file+" "+
            temporary_dir+"tmp.dvi");

        std::string info = file::to_string(temporary_dir+"tmp.dat");
        vec1s spl = split(split(info, "\n")[0], ",");

        double bh = 0, bd = 0, bw = 0;
        if (!from_string(erase_end(split(spl[0], "=")[1], "pt"), bh) ||
            !from_string(erase_end(split(spl[1], "=")[1], "pt"), bd) ||
            !from_string(erase_end(split(spl[2], "=")[1], "pt"), bw)) {
            warning("reading ", con.file_name);
            warning("parsing $"+s+"$ on l.", con.l);
            warning("could not read LaTeX info from tmp.dat");
        }

        double h = 0; {
            std::ifstream svg(con.out_dir+out_file);
            std::string tline;

            bool found = false;
            while (std::getline(svg, tline)) {
                tline = trim(tline);

                replace_block(tline, "<svg height='", "'", [&](std::string ts) {
                    ts = erase_end(ts, "pt");
                    if (!found && from_string(ts, h)) {
                        found = true;
                    }

                    return "";
                });

                if (found) break;
            }

            if (!found) {
                warning("reading ", con.file_name);
                warning("parsing $"+s+"$ on l.", con.l);
                warning("could not read SVG file height");
            }
        }

        double d = round(h*bd/(bh+bd));

        ret = "<img class=\"latex\" style=\"vertical-align: -"+strn(d)+"pt\" "
            "src="+out_file+">";

        con.latex_cache.insert(std::make_pair(hsh, ret));
    } else {
        ret = iter->second;
    }

    return ret;
}

template<typename F>
std::string parse_latex_command(context_t& con, const std::string& os, const std::string& cmd,
    uint_t narg, F&& parse) {

    std::string s;

    uint_t p = os.find("\\"+cmd);
    uint_t p0 = 0;
    while (p != npos) {
        s += os.substr(p0, p-p0);

        p += cmd.size()+1;
        p0 = p;

        vec1s args(narg);
        bool found = true;
        uint_t op0 = p0;
        for (uint_t i : range(narg)) {
            if (os[p] != '{') {
                found = false;
                break;
            }

            uint_t nopen = 1;

            ++p;
            ++p0;

            while (p < os.size() && nopen != 0) {
                if (os[p] == '\\') ++p;
                if (os[p] == '{') ++nopen;
                if (os[p] == '}') --nopen;
                ++p;
            }

            phypp_check(nopen == 0, "wrong balance of '{' and '}' in \\", cmd, " command "
                "parsing ", con.file_name, " l.", con.l);

            args[i] = os.substr(p0, p-1-p0);
            p0 = p;
        }

        if (found) {
            s += parse(std::move(args));
        } else {
            p0 = op0;
        }

        p = os.find("\\"+cmd, p0);
    }

    s += os.substr(p0);

    return s;
}

std::string latexify(context_t& con, std::string s) {
    static const std::string strue = pygmentize_inline(con, "c++", "true");
    static const std::string sfalse = pygmentize_inline(con, "c++", "false");

    s = parse_latex_command(con, s, "cpptrue", 0, [&](vec1s v) {
        return strue;
    });

    s = parse_latex_command(con, s, "cppfalse", 0, [&](vec1s v) {
        return sfalse;
    });

    s = parse_latex_command(con, s, "phypp", 0, [&](vec1s v) {
        return "<i>phy</i><sub>++</sub>";
    });

    s = parse_latex_command(con, s, "cppinline", 1, [&](vec1s v) {
        return pygmentize_inline(con, "c++", v[0]);
    });

    s = replace(s, "``", "&ldquo;");
    s = replace(s, "''", "&rdquo;");

    s = replace_block(s, "$", "$", [&](std::string ts) {
        return latex_math(con, ts);
    });

    s = parse_latex_command(con, s, "href", 2, [&](vec1s v) {
        return "<a href="+v[0]+">"+v[1]+"</a>";
    });

    s = parse_latex_command(con, s, "emph", 1, [&](vec1s v) {
        return "<span class=\"emphasis\">"+v[0]+"</span>";
    });

    s = parse_latex_command(con, s, "texttt", 1, [&](vec1s v) {
        return "<span class=\"fixed-width\">"+v[0]+"</span>";
    });
    s = parse_latex_command(con, s, "itt", 1, [&](vec1s v) {
        return "";
    }) ;

    s = parse_latex_command(con, s, "rsec", 1, [&](vec1s v) {
        auto iter = con.section_index.find(v[0]);
        if (iter != con.section_index.end()) {
            std::string link = replace_block(iter->second, "__(", ")__", [](std::string) {
                return "";
            });

            return "<a class=\"section-link\" href=\""+link+"\"><img src=\"section-link.png\"></a>";
        } else {
            return v[0];
        }
    }) ;

    s = parse_latex_command(con, s, "footnote", 1, [&](vec1s v) {
        std::string key = std::string(1, con.footnotes.size()+'a');
        insert_cache(con.footnotes, key, latexify(con, v[0]));
        return "<a class=\"footnote\" href=\"#footnote-"+key+"\"><sup>"+key+"</sup></a>";
    });

    s = replace(s, ".~", ". ");

    return s;
}

struct text_line {
    enum type {
        text,
        code,
        paragraph,
        begin,
        end,
        item
    };

    text_line(type t, std::string s = "") : ty(t), content(std::move(s)) {}

    type ty;
    std::string content;
};

using latex_lines = std::vector<text_line>;

void latexify(context_t& con, latex_lines& lines) {
    for (auto& d : lines) {
        if (d.ty == text_line::text) {
            std::string tline = trim(d.content);
            if (start_with(tline, "\\begin{itemize}")) {
                d.content = "<ul>";
                d.ty = text_line::begin;
            } else if (start_with(tline, "\\end{itemize}")) {
                d.content = "</ul>";
                d.ty = text_line::end;
            } else if (start_with(tline, "\\item")) {
                d.content = replace(d.content, "\\item", "<li>");
                d.content += "</li>";
                d.ty = text_line::item;
            }
        }
    }

    for (uint_t i : range(lines)) {
        auto& d = lines[i];

        if (d.ty == text_line::paragraph) {
            if (i != 0 && i != lines.size()-1 &&
                lines[i-1].ty == text_line::text &&
                lines[i+1].ty == text_line::text) {
                d.content = "<br><br>\n";
            }
        } else {
            if (d.ty != text_line::code) {
                d.content = latexify(con, d.content);
            }

            d.content = d.content + "\n";
        }
    }
}
