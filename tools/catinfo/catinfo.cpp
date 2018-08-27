#include <vif.hpp>

int vif_main(int argc, char* argv[]) {
    if (argc < 2) return 0;

    bool version = false;
    std::string fcat = argv[1];
    read_args(argc-1, argv+1, arg_list(version));

    struct {
        std::string catalogs;
        std::string comments;
        std::string version;
        std::string changelog;
    } cat;
        
    fits::read_table_loose(fcat, ftable(
        cat.catalogs, cat.comments, cat.version, cat.changelog
    ));

    print(fcat+"\n");

    if (version) {
        if (!cat.version.empty()) {
            print("version: "+cat.version);
        }
        if (!cat.changelog.empty()) {
            print(cat.changelog);
        }
    } else {
        if (!cat.comments.empty()) {
            print(cat.comments);
        }
        if (!cat.catalogs.empty()) {
            print(cat.catalogs);
        }
    }

    return 0;
}
