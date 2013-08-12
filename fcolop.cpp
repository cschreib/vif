#include <phypp.hpp>

void print_help() {
    using namespace format;

    print("fcolop v1.0");
}

int main(int argc, char* argv[]) {
    if (argc < 4) {
        print_help();
        return 0;
    }

    std::string file = argv[1];
    std::string op = tolower(argv[2]);
    vec1s cols;
    bool force = false;

    read_args(argc-2, argv+2, arg_list(cols, force));

    if (cols.empty()) {
        error("missing column(s) name(s)");
        print_help();
        return 0;
    }

    cols = toupper(cols);

    fitsfile* fptr;
    int status = 0;
    fits_open_table(&fptr, file.c_str(), READWRITE, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+file+"'");

    vec1s fcols;
    for (auto& col : cols) {
        int id = -1;
        fits_get_colnum(fptr, CASESEN, const_cast<char*>(col.c_str()), &id, &status);
        if (status != 0) {
            status = 0;
            if (!force) {
                std::cout << "warning: no column named '" << col << "' in this file, continue? (y/n) "
                    << std::flush;

                std::string line;
                while (line.empty()) {
                    std::getline(std::cin, line);
                    line = tolower(line);
                    if (line == "y" || line == "yes") break;
                    if (line == "n" || line == "no") {
                        fits_close_file(fptr, &status);
                        return 0;
                    }
                    std::cout << "please answer either 'yes' (y) or 'no' (n): " << std::flush;
                }
            }

            continue;
        }

        fcols.push_back(col);
    }

    if (op == "remove") {
        if (!force) {
            std::cout << "warning: this operation will permanently remove information from this file\n";
            std::cout << "warning: are you sure ? (y/n) " << std::flush;
            std::string line;
            while (line.empty()) {
                std::getline(std::cin, line);
                line = tolower(line);
                if (line == "y" || line == "yes") break;
                if (line == "n" || line == "no") {
                    fits_close_file(fptr, &status);
                    return 0;
                }
                std::cout << "please answer either 'yes' (y) or 'no' (n): " << std::flush;
            }
        }

        for (auto& col : fcols) {
            int id;
            fits_get_colnum(fptr, CASESEN, const_cast<char*>(col.c_str()), &id, &status);

            if (status != 0) {
                status = 0;
                continue;
            }

            fits_delete_col(fptr, id, &status);
            fits::phypp_check_cfitsio(status, "cannot remove column '"+col+"'");
        }
    } else if (op == "transpose") {
        for (auto& col : fcols) {
            vec2d v;

            int id;
            fits_get_colnum(fptr, CASESEN, const_cast<char*>(col.c_str()), &id, &status);

            char otform[80], comment[80];
            std::string tformn = "TFORM"+strn(id);
            fits_read_key(fptr, TSTRING, const_cast<char*>(tformn.c_str()), otform, comment, &status);
            fits::phypp_check_cfitsio(status, "cannot read TFORM for column '"+col+"'");
            print(otform);

            int type, naxis;
            long repeat, width;
            std::array<long,2> naxes;
            fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
            fits_read_tdim(fptr, id, 2, &naxis, naxes.data(), &status);
            phypp_check(naxis == 2, "wrong dimension for column '"+col+"' (", naxis, "), only "
                "2D columns can be transposed");
            phypp_check(type != TSTRING, "cannot transpose a string column");

            status = 0;

            v.dims[0] = naxes[1];
            v.dims[1] = naxes[0];

            v.resize();

            double def = dnan;
            int null;
            fits_read_col(
                fptr, fits::traits<double>::ttype, id, 1, 1, repeat, &def, v.data.data(),
                &null, &status
            );

            fits_delete_col(fptr, id, &status);
            fits::phypp_check_cfitsio(status, "cannot temporarily delete column '"+col+"'");

            v = transpose(v);

            fits_insert_col(fptr, id, const_cast<char*>(col.c_str()), otform, &status);
            fits::phypp_check_cfitsio(status, "cannot re-insert column '"+col+"'");

            std::swap(naxes[0], naxes[1]);
            fits_write_tdim(fptr, id, 2, naxes.data(), &status);
            fits::phypp_check_cfitsio(status, "cannot re-write TDIM for column '"+col+"'");

            fits_write_col(
                fptr, fits::traits<double>::ttype, id, 1, 1, n_elements(v),
                const_cast<vec2d::dtype*>(v.data.data()), &status
            );
            fits::phypp_check_cfitsio(status, "cannot write transposed data for column '"+col+"'");
        }
    } else {
        error("unknown operation '", op, "'");
    }

    fits_close_file(fptr, &status);

    return 0;
}
