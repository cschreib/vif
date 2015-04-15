#include <phypp.hpp>
#include "common.hpp"

struct section_t {
    std::string title;
    std::string id, fid;
    std::string funcfile;
    std::vector<text_line> desc;
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

    sorted_cache_t index;

    section_t mainsec;
    section_t* sec = &mainsec;
    std::vector<section_t> subsecs;
    bool open = true;
    bool incode = false;

    std::string code;

    std::string header = file::to_string("header.html");
    std::string footer = file::to_string("footer.html");
    std::string menu = file::to_string(con.out_dir+"menu.html");

    // Read cache, if any
    std::string cache_dir = con.out_dir+"cache/";
    file::mkdir(cache_dir);

    read_cache(con.pygment_cache, cache_dir+con.base+"_pygment.cache");
    read_cache(con.latex_cache, cache_dir+con.base+"_latex.cache");

    for (auto ifile : file::list_files(cache_dir+"category*_index.cache")) {
        read_cache(con.section_index, cache_dir+ifile);
    }

    std::string line;
    while (std::getline(file, line)) {
        ++con.l;

        line = parse_latex_command(con, line, "texorhtml", 2, [](vec1s v) {
            return v[1];
        });

        std::string tline = trim(line);

        if (start_with(tline, "\\section")) {
            // Locate the ID of the section
            line = regex_replace(line, R"(\\label\{([^}]+)\})", [&](vec1s s) {
                mainsec.id = s[0];
                return "";
            });

            if (mainsec.id.empty()) {
                // Not found
                error("reading ", con.file_name);
                error("missing section ID l.", con.l);
                return 1;
            }

            // Locate the title
            vec2s p = regex_extract(line, R"(\\section\{([^}]+)\})");
            if (!p.empty()) {
                if (p.dims[0] != 1) {
                    error("reading ", con.file_name);
                    error("multiple \\section on l.", con.l);
                    return 1;
                }

                mainsec.title = trim(p[0]);
            } else {
                // Not found
                error("reading ", con.file_name);
                error("missing section name l.", con.l);
                return 1;
            }

            insert_cache(index, sec->id, "__("+sec->title+")__"+con.base+".html");
        } else if (start_with(tline, "\\subsection")) {
            section_t sub;

            // Locate the ID of the section
            line = regex_replace(line, R"(\\label\{([^}]+)\})", [&](vec1s s) {
                sub.id = s[0];
                sub.fid = replace(sub.id, ":", "_");
                return "";
            });

            if (sub.id.empty()) {
                // Not found
                warning("reading ", con.file_name);
                warning("missing subsection ID l.", con.l);
            }

            // Locate the title
            vec2s p = regex_extract(line, R"(\\subsection\{([^}]+)\})");
            if (!p.empty()) {
                if (p.dims[0] != 1) {
                    error("reading ", con.file_name);
                    error("multiple \\subsection on l.", con.l);
                    return 1;
                }

                sub.title = trim(p[0]);
            } else {
                // Not found
                error("reading ", con.file_name);
                error("missing subection name l.", con.l);
                return 1;
            }

            subsecs.push_back(sub);
            sec = &subsecs.back();

            insert_cache(index, sec->id, "__("+sec->title+")__"+con.base+".html#"+sec->fid);
        } else if (start_with(tline, "\\loadfunctions")) {
            vec2s p = regex_extract(line, R"(\\loadfunctions\{([^}]+)\})");
            if (!p.empty()) {
                if (p.dims[0] != 1) {
                    error("reading ", con.file_name);
                    error("multiple \\loadfunctions on l.", con.l);
                    return 1;
                }

                sec->funcfile = p[0];
            }
        } else if (start_with(tline, "\\begin{cppcode}")) {
            incode = true;
        } else if (start_with(tline, "\\end{cppcode}")) {
            sec->desc.push_back(text_line(text_line::code, pygmentize(con, "c++", code)));
            code = "";
            incode = false;
        } else if (start_with(tline, "\\begin{bashcode}")) {
            incode = true;
        } else if (start_with(tline, "\\end{bashcode}")) {
            sec->desc.push_back(text_line(text_line::code, pygmentize(con, "bash", code)));
            code = "";
            incode = false;
        } else if (start_with(tline, "\\begin{example}")) {
            sec->desc.push_back(text_line(text_line::begin, "<br><span class=\"example\">Example</span><br><hr>"));
        } else if (start_with(tline, "\\end{example}")) {
            sec->desc.push_back(text_line(text_line::end, "<hr><br>"));
        } else {
            if (incode) {
                if (!code.empty()) code += "\n";
                code += line;
            } else if (!empty(line)) {
                open = false;
                sec->desc.push_back(text_line(text_line::text, line));
            } else {
                sec->desc.push_back(text_line(text_line::paragraph));
            }
        }
    }

    if (incode) {
        error("reading ", con.file_name);
        error("reached end of file looking for \\end{incode}");
        return 1;
    }

    auto write_section = [&](std::ofstream& out, section_t& s) {
        out << "<div class=\"section-content\">\n";

        if (!s.desc.empty()) {
            bool empty = true;
            latexify(con, s.desc);

            for (auto& d : s.desc) {
                if (d.ty != text_line::paragraph) {
                    empty = false;
                    break;
                }
            }

            if (!empty) {
                out << "<div class=\"description\">\n";

                for (auto d : s.desc) {
                    out << d.content;
                }

                out << "</div>\n";
            }
        }

        // Load index cache of function list
        std::string ffile = con.out_dir+"cache/"+file::remove_extension(s.funcfile)+"_index.cache";
        if (file::exists(ffile)) {
            out << "<ul class=\"category-funclist\">\n";
            cache_t findex;
            read_cache(findex, ffile);
            std::string last;
            for (auto kv : findex) {
                if (kv.key != last) {
                    out << "<li><a href=\"" << kv.value << "\">" << kv.key << "</a></li>\n";
                    last = kv.key;
                }
            }
            out << "</ul>\n";
        }

        out << "</div>\n";
    };

    std::ofstream out(con.out_dir+con.base+".html");
    out << header << std::endl;

    out << "<div class=\"content\">\n";

    out << "<span class=\"section\">"
        << mainsec.title
        << "</span>\n";

    write_section(out, mainsec);

    for (auto s : subsecs) {
        out << "<span class=\"subsection\"><a name=\"" << s.fid << "\"></a>"
            << s.title
            << "</span>\n";

        write_section(out, s);
    }

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

    // Save cache
    write_cache(con.pygment_cache, cache_dir+con.base+"_pygment.cache");
    write_cache(con.latex_cache, cache_dir+con.base+"_latex.cache");
    write_cache(index, cache_dir+con.base+"_index.cache");

    return 0;
}
