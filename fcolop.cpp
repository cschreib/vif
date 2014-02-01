#include <phypp.hpp>

void print_help();

bool get_columns(const vec1s& cols, fitsfile* fptr, vec1s& fcols, bool force = false);

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string file = argv[1];
    if (!file::exists(file)) {
        error("cannot open file '"+file+"'");
        return 1;
    }

    std::string op = tolower(argv[2]);

    fitsfile* fptr;
    int status = 0;
    fits_open_table(&fptr, file.c_str(), READWRITE, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+file+"'");

    if (op == "remove") {
        vec1s cols;
        bool force = false;
        read_args(argc-2, argv+2, arg_list(cols, force));

        if (cols.empty()) {
            error("missing column(s) name(s)");
            print_help();
            return 0;
        }

        vec1s fcols;
        cols = toupper(cols);
        if (!get_columns(cols, fptr, fcols, force)) {
            return 0;
        }

        if (!force) {
            std::cout << "warning: this operation will permanently remove information from this file\n";
            std::cout << "warning: are you sure? (y/n) " << std::flush;
            std::string line;
            while (line.empty()) {
                std::getline(std::cin, line);
                line = tolower(line);
                if (line == "y" || line == "yes") break;
                if (line == "n" || line == "no") {
                    fits_close_file(fptr, &status);
                    return 0;
                }
                line.clear();
                std::cout << "error: please answer either 'yes' (y) or 'no' (n): " << std::flush;
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
        vec1s cols;
        bool force = false;
        read_args(argc-2, argv+2, arg_list(cols, force));

        if (cols.empty()) {
            error("missing column(s) name(s)");
            print_help();
            return 0;
        }

        vec1s fcols;
        cols = toupper(cols);
        if (!get_columns(cols, fptr, fcols, force)) {
            return 0;
        }

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
    } else if (op == "meta") {
        bool force = false;
        read_args(argc-2, argv+2, arg_list(force));

        vec1u dims;
        int ncols;
        fits_get_num_cols(fptr, &ncols, &status);

        if (ncols == 0) {
            print("no columns in this file");
            fits_close_file(fptr, &status);
            return 0;
        }

        for (int i = 1; i <= ncols; ++i) {
            int naxis;
            std::array<long,100> naxes;
            fits_read_tdim(fptr, i, 100, &naxis, naxes.data(), &status);
            for (int j = 0; j < naxis; ++j) {
                dims.push_back(naxes[j]);
            }
        }

        vec1u udim = dims[uniq(dims, sort(dims))];
        uint_t dim;
        if (udim.size() == 1) {
            dim = udim[0];
        } else {
            vec1u count(udim.size());
            for (uint_t i = 0; i < udim.size(); ++i) {
                count[i] = total(dims == udim[i]);
            }

            vec1u weight = uindgen(udim.size()) + sort(count);
            dim = udim[max_id(weight)];

            if (!force) {
                std::cout << "identified 'row' dimension: " << dim << "\n";
                std::cout << "is it correct? (y/n) " << std::flush;
                std::string line;
                while (line.empty()) {
                    std::getline(std::cin, line);
                    line = tolower(line);
                    if (line == "y" || line == "yes") break;
                    if (line == "n" || line == "no") {
                        print("available dimensions:");
                        print(udim);
                        std::cout << "please enter the correct 'row' dimension: ";
                        line.clear();
                        while (line.empty()) {
                            std::getline(std::cin, line);
                            if (from_string(line, dim) && total(dim == udim) != 0) {
                                break;
                            }

                            error("please enter a number that matches one of the following "
                                "dimensions:");
                            print(udim);
                            line.clear();
                            std::cout << "please enter the correct 'row' dimension: ";
                        }
                        break;
                    }
                    line.clear();
                    std::cout << "error: please answer either 'yes' (y) or 'no' (n): " << std::flush;
                }
            }
        }

        vec1u mcols;
        for (int i = 1; i <= ncols; ++i) {
            int naxis;
            std::array<long,100> naxes;
            fits_read_tdim(fptr, i, 100, &naxis, naxes.data(), &status);
            bool found = false;
            for (int j = 0; j < naxis; ++j) {
                if (uint_t(naxes[j]) == dim) {
                    found = true;
                    break;
                }
            }

            if (!found) {
                mcols.push_back(i);
            }
        }

        if (mcols.empty()) {
            print("no meta data column");
            fits_close_file(fptr, &status);
            return 0;
        }

        if (!force) {
            vec1s names;
            for (auto i : mcols) {
                std::string id = strn(i);
                char name[80];
                int num;
                fits_get_colname(fptr, CASESEN, const_cast<char*>(id.c_str()), name, &num, &status);
                names.push_back(name);
            }

            warning("this operation will move the following columns in a new table:");
            print(names);
            std::cout << "warning: are you sure? (y/n) " << std::flush;
            std::string line;
            while (line.empty()) {
                std::getline(std::cin, line);
                line = tolower(line);
                if (line == "y" || line == "yes") break;
                if (line == "n" || line == "no") {
                    fits_close_file(fptr, &status);
                    return 0;
                }
                line.clear();
                std::cout << "error: please answer either 'yes' (y) or 'no' (n): " << std::flush;
            }
        }

        fitsfile* ofptr;
        fits_open_table(&ofptr, file.c_str(), READWRITE, &status);

        fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

        for (auto i : mcols) {
            fits_copy_col(ofptr, fptr, i, mcols.size(), 1, &status);
            fits::phypp_check_cfitsio(status, "cannot copy column "+strn(i));
        }

        for (auto i : reverse(mcols)) {
            fits_delete_col(ofptr, i, &status);
        }
    } else {
        error("unknown operation '", op, "'");
    }

    fits_close_file(fptr, &status);

    return 0;
}

bool get_columns(const vec1s& cols, fitsfile* fptr, vec1s& fcols, bool force) {
    int status = 0;
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
                        return false;
                    }
                    line.clear();
                    std::cout << "error: please answer either 'yes' (y) or 'no' (n): " << std::flush;
                }
            }

            continue;
        }

        fcols.push_back(col);
    }

    return true;
}

void print_help() {
    using namespace format;

    print("fcolop v1.0");
    header("Usage: fcolop op cols=[...] [force]");
    header("Available commands:");
    bullet("remove", "remove the specified columns from this FITS file");
    bullet("transpose", "transpose the specified columns within this FITS file");
    bullet("meta", "shift all 'meta' columns into a new extension so that the file can be read in "
        "programs like TOPCAT that expect a single invariant dimension for all columns");
}
