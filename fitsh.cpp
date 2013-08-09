#include <phypp.hpp>

void print_help() {
    print("fitsh v1.0");
    print("  usage: fitsh file [options]\n");
    print("  Available options:");
    print("    edit: start interactive mode, to edit some of the keywords");
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        print_help();
        return 0;
    }

    bool edit = false;
    read_args(argc-1, argv+1, arg_list(edit));

    if (edit) {
        fitsfile* fptr;
        int status = 0;

        fits_open_image(&fptr, argv[1], READWRITE, &status);
        fits::phypp_check_cfitsio(status, "cannot open file '"+std::string(argv[1])+"'");

        while (true) {
            print("please type the name of the keyword you want to edit (empty to exit):");
            std::string name;
            std::cout << "> ";
            std::getline(std::cin, name);
            name = toupper(trim(name));
            if (name.empty()) break;

            char value[80];
            char comment[80];
            fits_read_keyword(fptr, const_cast<char*>(name.c_str()), value, comment, &status);

            print(name, " = ", value, " (", comment, ")");
            print("please type the new value for this keyword (empty to abort):");
            std::string new_value;
            std::cout << "> ";
            std::getline(std::cin, new_value);
            new_value = trim(new_value);
            if (new_value.empty()) break;

            if (new_value[0] == '\'') {
                new_value.erase(0,1);
                for (uint_t i = new_value.size()-1; i != npos; --i) {
                    if (new_value[i] == '\'') {
                        new_value.resize(i);
                        break;
                    }
                }

                fits_update_key(fptr, TSTRING, const_cast<char*>(name.c_str()),
                    const_cast<char*>(new_value.c_str()), comment, &status);
            } else {
                double d;
                from_string(new_value, d);
                fits_update_key(fptr, TDOUBLE, const_cast<char*>(name.c_str()), &d, comment,
                    &status);
            }
        }

        fits_close_file(fptr, &status);
    } else {
        std::string header = fits::read_header(argv[1]);
        for (uint_t i = 0; i < header.size(); i += 80) {
            print(header.substr(i, std::min(header.size()-i, std::size_t(80))));
        }
    }

    return 0;
}
