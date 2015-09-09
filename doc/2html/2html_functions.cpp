#include <phypp.hpp>
#include "common.hpp"

struct function_t {
    std::string name;
    std::string fname;
    std::string signature;
    vec1s requires;
    bool vectorized = false;
};

int main(int argc, char* argv[]) {
    if (argc < 2) {
        return 1;
    }

    context_t con;
    con.file_name = argv[1];
    if (!file::exists(con.file_name)) {
        error("could not find '", con.file_name, "'");
        return 1;
    }

    read_args(argc-1, argv+1, arg_list(con.out_dir));

    if (!con.out_dir.empty()) {
        con.out_dir = file::directorize(con.out_dir);
        file::mkdir(con.out_dir);
    }

    file::mkdir(con.out_dir+"latex/");

    con.base = file::get_basename(file::remove_extension(con.file_name));
    std::ifstream file(con.file_name);

    cache_t index;

    std::vector<function_t> funcs;
    bool open = true;
    bool incode = false;

    latex_lines desc;
    std::string code;

    std::string header = file::to_string("header.html");
    std::string footer = file::to_string("footer.html");

    // Read cache, if any
    std::string cache_dir = con.out_dir+"cache/";
    file::mkdir(cache_dir);

    read_cache(con.pygment_cache, cache_dir+con.base+"_pygment.cache");
    read_cache(con.latex_cache, cache_dir+con.base+"_latex.cache");

    for (auto ifile : file::list_files(cache_dir+"category*_index.cache")) {
        read_cache(con.section_index, cache_dir+ifile);
    }

    // Find the categories associated with this file
    std::vector<std::pair<std::string, std::string>> categories; {
        sorted_cache_t catfun_cache;
        for (auto ifile : file::list_files(cache_dir+"category*_catfun.cache")) {
            read_cache(catfun_cache, cache_dir+ifile);
        }

        auto iter = catfun_cache.find(file::get_basename(con.file_name));
        if (iter != catfun_cache.end()) {
            vec1s spl = split(iter->second, "|||");
            for (auto& ts : spl) {
                std::string tit, link;
                link = replace_block(ts, "__(", ")__", [&](std::string s) {
                    tit = s;
                    return "";
                });

                categories.push_back(std::make_pair(tit, link));
            }
        }
    }


    std::string line;

    auto find_function = [&]() {
        function_t f;

        // Locate the name of the function
        vec2s p = regex_extract(line, R"(\\itt\{([a-zA-Z0-9_:]+)\})");
        if (!p.empty()) {
            if (p.dims[0] != 1) {
                error("reading ", con.file_name);
                error("too many \\itt{...} on l.", con.l);
                return false;
            }

            f.name = p[0];
        } else {
            // Not found, use that of the previous function in the list
            if (funcs.empty()) {
                error("reading ", con.file_name);
                error("missing function name l.", con.l);
                return false;
            }

            f.name = funcs.back().name;
        }

        // Create a "file friendly" name
        f.fname = replace(f.name, "::", "_");

        // Locate the signature
        p = regex_extract(line, R"(\\cppinline\|([^|]+)\|)");
        if (!p.empty()) {
            if (p.dims[0] != 1) {
                error("reading ", con.file_name);
                error("multiple \\cppinline on l.", con.l);
                return false;
            }

            f.signature = p[0];
        } else {
            // This must be a multi-line signature
            uint_t p0 = line.find("\\begin{cppcode}");
            if (p0 == npos) {
                error("reading ", con.file_name);
                error("could not extract signature from l.", con.l);
                return false;
            }

            f.signature = line.substr(p0+std::string("\\begin{cppcode}").size());

            bool closed = false;
            while (std::getline(file, line)) {
                if (start_with(trim(line), "\\end{cppcode}")) {
                    closed = true;
                    break;
                } else {
                    f.signature += "\n"+line;
                }
            }

            if (!closed) {
                error("reading ", con.file_name);
                error("reached end of file looking for \\end{cppcode}");
                return false;
            }
        }

        // Check if the function is vectorized
        f.vectorized = find(line, "\\vectorfunc") != npos;

        // Find dependencies
        vec2s ext = regex_extract(line, R"(\\requirelib\{([^}]+)\})");
        if (!ext.empty()) {
            f.requires = ext(0,_);
        }

        funcs.push_back(f);

        insert_cache(index, f.name, con.base+"_"+funcs.front().fname+".html");

        return true;
    };

    auto write_function = [&]() {
        if (!funcs.empty()) {
            std::ofstream out(con.out_dir+con.base+"_"+funcs[0].fname+".html");
            out << header << std::endl;

            out << "<div class=\"content\">\n";

            if (!categories.empty()) {
                out << "<span class=\"section\">Categories</span><br>\n"
                    << "<div class=\"catlist\">\n";

                vec1s scats;

                for (auto& c : categories) {
                    scats.push_back("<a href=\""+c.second+"\">"+c.first+"</a>");
                }

                out << collapse(scats, "<span class=\"catsep\">&gt;</span>") << "\n";

                out << "</div><br>\n";
            }

            out << "<span class=\"section\">Signature"
                << (funcs.size() > 1 ? "s" : "")
                << "</span>\n";
            out << "<table style=\"signatures\">\n";

            bool hasdep = false;
            bool hasvec = false;
            for (auto& f : funcs) {
                if (!f.requires.empty()) {
                    hasdep = true;
                }
                if (f.vectorized) {
                    hasvec = true;
                }
            }

            for (auto& f : funcs) {
                out << "<tr>";
                if (hasdep) {
                    out << "<td>";
                    out << collapse("<span class=\"libsymbol\">["+f.requires+"]</span>", " ");
                    out << "</td>";
                }
                if (hasvec) {
                    out << "<td>";
                    if (f.vectorized) {
                        out << "<img src=\"vectorized.png\" alt=\"vectorized\">";
                    } else {
                        out << "<img src=\"vectorized.png\" alt=\"vectorized\" class=\"hidden\">";
                    }
                    out << "</td>";
                }
                std::string sig = pygmentize(con, "c++", f.signature);
                out << "<td>" << replace(sig, "class=\"highlight\"", "class=\"highlight-sig\"") << "</td>";
                out << "</tr>\n";
            }
            out << "</table>\n";

            out << "<br><span class=\"section\">Description</span><br>\n";
            out << "<div class=\"description\">\n";

            latexify(con, desc);
            for (auto& d : desc) {
                out << d.content;
            }

            out << "</div>\n";

            if (!con.footnotes.empty()) {
                out << "<hr>\n";
                out << "<ul class=\"footnotes\">\n";
                for (auto& f : con.footnotes) {
                    out << "<li><a name=\"footnote-" << f.first << "\"></a>"
                        << "<span class=\"footnote\"><sup>" << f.first << "</sup></span> "
                        << f.second << "</li>\n";
                }
                out << "</ul>\n";
            }

            out << "</div>\n";

            out << footer << std::endl;
            out.close();
        }

        con.footnotes.clear();
        funcs.clear();
        desc.clear();
        open = true;

        return true;
    };

    while (std::getline(file, line)) {
        ++con.l;

        std::string tline = trim(line);
        if (start_with(tline, "\\funcitem")) {
            if (!write_function()) {
                return 1;
            }

            if (!find_function()) {
                return 1;
            }
        } else if (open && (start_with(tline, "\\vectorfunc") ||
                            start_with(tline, "\\requirelib") ||
                            start_with(tline, "\\cppinline"))) {
            if (!find_function()) {
                return 1;
            }
        } else if (start_with(tline, "\\begin{cppcode}")) {
            incode = true;
        } else if (start_with(tline, "\\end{cppcode}")) {
            desc.push_back(text_line(text_line::code, pygmentize(con, "c++", code)));
            code = "";
            incode = false;
        } else if (start_with(tline, "\\begin{bashcode}")) {
            incode = true;
        } else if (start_with(tline, "\\end{bashcode}")) {
            desc.push_back(text_line(text_line::code, pygmentize(con, "bash", code)));
            code = "";
            incode = false;
        } else if (start_with(tline, "\\begin{example}")) {
            desc.push_back(text_line(text_line::begin, "<br><br><span class=\"example\">Example</span><br><hr>"));
        } else if (start_with(tline, "\\end{example}")) {
            desc.push_back(text_line(text_line::end, "<hr><br>"));
        } else if (start_with(tline, "\\begin{advanced}")) {
            desc.push_back(text_line(text_line::begin, "<br><br><span class=\"advanced-title\">Advanced</span>"
                "<div class=\"advanced\"><hr>"));
        } else if (start_with(tline, "\\end{advanced}")) {
            desc.push_back(text_line(text_line::end, "<hr></div><br>"));
        } else {
            if (incode) {
                if (!code.empty()) code += "\n";
                code += line;
            } else if (!empty(line)) {
                open = false;
                desc.push_back(text_line(text_line::text, line));
            } else {
                desc.push_back(text_line(text_line::paragraph, line));
            }
        }
    }

    if (incode) {
        error("reading ", con.file_name);
        error("reached end of file looking for \\end{incode}");
        return 1;
    }

    if (!write_function()) {
        return 1;
    }

    // Save cache
    write_cache(con.pygment_cache, cache_dir+con.base+"_pygment.cache");
    write_cache(con.latex_cache, cache_dir+con.base+"_latex.cache");
    write_cache(index, cache_dir+con.base+"_index.cache");

    return 0;
}
