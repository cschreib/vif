#include <phypp.hpp>
#include "common.hpp"

int main(int argc, char* argv[]) {
    std::string out_dir = file::directorize(argv[1]);

    sorted_cache_t funindex;
    for (auto& f : file::list_files(out_dir+"cache/functions*_index.cache")) {
        read_cache(funindex, out_dir+"cache/"+f);
    }

    std::vector<cache_t> catindex;
    for (auto& f : file::list_files(out_dir+"cache/category*_index.cache")) {
        cache_t tmp;
        read_cache(tmp, out_dir+"cache/"+f);
        catindex.push_back(std::move(tmp));
    }

    print("found ", funindex.size(), " functions");
    print("found ", catindex.size(), " categories");

    std::string header = file::to_string("header.html");
    std::string footer = file::to_string("footer.html");

    std::ofstream out(out_dir+"menu-left.js");

    out << "document.write('\\\n";
    out << "<ul class=\"menu-left\"> \\\n";
    out << "<li><span class=\"menu-heading\">Categories</span><br><hr> \\\n";
    out << "<ul class=\"menu-cat\"> \\\n";
    for (auto& c : catindex) {
        if (c.empty()) continue;

        bool first = true;
        for (auto& i : c) {
            std::string tit, link;
            link = replace_block(i.value, "__(", ")__", [&](std::string s) {
                tit = s;
                return "";
            });

            if (c.size() != 1) {
                if (first) {
                    out << "<li><a class=\"menu-cat-sec\" href=\""+link+"\">"+tit+"</a> \\\n";
                    out << "<ul class=\"menu-cat-sub\"> \\\n";
                } else {
                    out << "<li><a href=\""+link+"\">"+tit+"</a></li> \\\n";
                }
            } else {
                out << "<li><a class=\"menu-cat-sec\" href=\""+link+"\">"+tit+"</a></li> \\\n";
            }

            first = false;
        }

        if (c.size() != 1) {
            out << "</ul></li> \\\n";
        }
    }

    out << "</ul></li> \\\n";
    out << "</ul> \\\n";
    out << "');";

    out.close();
    out.open(out_dir+"menu-right.js");

    out << "document.write('\\\n";
    out << "<ul class=\"menu-right\"> \\\n";
    out << "<li><span class=\"menu-heading\">Alphabetical list</span><br><hr> \\\n";
    out << "<ul class=\"menu-alpha\"> \\\n";
    char last = '\0';
    for (auto& kv : funindex) {
        if (kv.first[0] != last) {
            if (last != '\0') {
                out << "</ul><br></li> \\\n";
            }

            out << "<li><span class=\"menu-alpha-heading\">"+toupper(kv.first).substr(0,1)+"</span><br> \\\n";
            out << "<ul> \\\n";

            last = kv.first[0];
        }

        out << "<li><a href=\""+kv.second+"\">"+kv.first+"</a></li> \\\n";
    }

    out << "</ul></li> \\\n";
    out << "</ul></li> \\\n";
    out << "</ul> \\\n";
    out << "');";

    return 0;
}
