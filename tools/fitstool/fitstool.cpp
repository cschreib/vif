#include <vif.hpp>

using namespace vif;

void print_help();

bool get_columns(const vec1s& cols, fitsfile* fptr, vec1s& fcols, bool force = false);

bool remove_columns(int argc, char* argv[], const std::string& file);
bool transpose_columns(int argc, char* argv[], const std::string& file);
bool rows_to_columns(int argc, char* argv[], const std::string& file);
bool copy_wcs_header(int argc, char* argv[], const std::string& file);
bool show_header(int argc, char* argv[], const std::string& file);
bool edit_keyword(int argc, char* argv[], const std::string& file);
bool read_keyword(int argc, char* argv[], const std::string& file);
bool remove_keyword(int argc, char* argv[], const std::string& file);
bool make_2d(int argc, char* argv[], const std::string& file);
bool move_meta_columns(int argc, char* argv[], const std::string& file);
bool extract_extension(int argc, char* argv[], const std::string& file);

void print_remove_help();
void print_transpose_help();
void print_r2c_help();
void print_cpwcs_help();
void print_hdr_help();
void print_editkwd_help();
void print_rmkwd_help();
void print_readkwd_help();
void print_make2d_help();
void print_meta_help();
void print_extract_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string file = argv[1];
    std::string op = to_lower(argv[2]);

    if (op == "help") {
        op = file;
        if (op == "remove") {
            print_remove_help();
        } else if (op == "transpose") {
            print_transpose_help();
        } else if (op == "r2c") {
            print_r2c_help();
        } else if (op == "cpwcs") {
            print_cpwcs_help();
        } else if (op == "hdr") {
            print_hdr_help();
        } else if (op == "editkwd") {
            print_editkwd_help();
        } else if (op == "rmkwd") {
            print_rmkwd_help();
        } else if (op == "readkwd") {
            print_readkwd_help();
        } else if (op == "make2d") {
            print_make2d_help();
        } else if (op == "meta") {
            print_meta_help();
        } else if (op == "extract") {
            print_extract_help();
        } else {
            error("unknown operation '", op, "'");
        }
    } else {
        if (!file::exists(file)) {
            error("cannot open file '"+file+"'");
            return 1;
        }

        if (op == "remove") {
            remove_columns(argc-2, argv+2, file);
        } else if (op == "transpose") {
            transpose_columns(argc-2, argv+2, file);
        } else if (op == "r2c") {
            rows_to_columns(argc-2, argv+2, file);
        } else if (op == "cpwcs") {
            copy_wcs_header(argc-2, argv+2, file);
        } else if (op == "hdr") {
            show_header(argc-2, argv+2, file);
        } else if (op == "editkwd") {
            edit_keyword(argc-2, argv+2, file);
        } else if (op == "rmkwd") {
            remove_keyword(argc-2, argv+2, file);
        } else if (op == "readkwd") {
            read_keyword(argc-2, argv+2, file);
        } else if (op == "make2d") {
            make_2d(argc-2, argv+2, file);
        } else if (op == "meta") {
            move_meta_columns(argc-2, argv+2, file);
        } else if (op == "extract") {
            extract_extension(argc-2, argv+2, file);
        } else {
            error("unknown operation '", op, "'");
        }
    }

    return 0;
}

void print_remove_help() {
    using namespace terminal_format;

    paragraph("The program will remove the provided columns from this FITS file. Be warned that "
        "the contained data is permantently lost.");

    header("List of available command line options:");
    bullet("cols", "[strings] list of columns to remove");
    bullet("force", "[flag] do not ask for confirmation at each step");
    bullet("help", "[flag] print this text");
    print("");
}

bool remove_columns(int argc, char* argv[], const std::string& file) {
    vec1s cols;
    bool force = false;
    bool help = false;
    read_args(argc, argv, arg_list(cols, force, help));

    if (help) {
        print_remove_help();
        return true;
    }

    if (cols.empty()) {
        error("missing column(s) name(s)");
        print_remove_help();
        return false;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_table(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        vec1s fcols;
        cols = to_upper(cols);
        if (!get_columns(cols, fptr, fcols, force)) {
            return true;
        }

        if (!force) {
            std::cout << "warning: this operation will permanently remove information from this file\n";
            std::cout << "warning: are you sure? (y/n) " << std::flush;
            std::string line;
            while (line.empty()) {
                std::getline(std::cin, line);
                line = to_lower(line);
                if (line == "y" || line == "yes") break;
                if (line == "n" || line == "no") {
                    fits_close_file(fptr, &status);
                    return true;
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
            fits::vif_check_cfitsio(status, "cannot remove column '"+col+"'");
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_transpose_help() {
    using namespace terminal_format;

    paragraph("The program will transpose the provided 2D columns, effectively transforming a NxM "
        "column into an MxN one. No data is lost in the process.");

    header("List of available command line options:");
    bullet("cols", "[strings] list of columns to transpose");
    bullet("force", "[flag] do not ask for confirmation at each step");
    bullet("help", "[flag] print this text");
    print("");
}

bool transpose_columns(int argc, char* argv[], const std::string& file) {
    vec1s cols;
    bool force = false;
    bool help = false;
    read_args(argc, argv, arg_list(cols, force, help));

    if (help) {
        print_transpose_help();
        return true;
    }

    if (cols.empty()) {
        error("missing column(s) name(s)");
        print_transpose_help();
        return false;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_table(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        vec1s fcols;
        cols = to_upper(cols);
        if (!get_columns(cols, fptr, fcols, force)) {
            return true;
        }

        for (auto& col : fcols) {
            vec2d v;

            int id;
            fits_get_colnum(fptr, CASESEN, const_cast<char*>(col.c_str()), &id, &status);

            char otform[80], comment[80];
            std::string tformn = "TFORM"+to_string(id);
            fits_read_key(fptr, TSTRING, const_cast<char*>(tformn.c_str()), otform, comment, &status);
            fits::vif_check_cfitsio(status, "cannot read TFORM for column '"+col+"'");
            print(otform);

            int type, naxis;
            long repeat, width;
            std::array<long,2> naxes;
            fits_get_coltype(fptr, id, &type, &repeat, &width, &status);
            fits_read_tdim(fptr, id, 2, &naxis, naxes.data(), &status);
            vif_check(naxis == 2, "wrong dimension for column '"+col+"' (", naxis, "), only "
                "2D columns can be transposed");
            vif_check(type != TSTRING, "cannot transpose a string column");

            status = 0;

            v.dims[0] = naxes[1];
            v.dims[1] = naxes[0];

            v.resize();

            double def = dnan;
            int null;
            fits_read_col(
                fptr, impl::fits_impl::traits<double>::ttype, id, 1, 1, repeat, &def, v.data.data(),
                &null, &status
            );

            fits_delete_col(fptr, id, &status);
            fits::vif_check_cfitsio(status, "cannot temporarily delete column '"+col+"'");

            v = transpose(v);

            fits_insert_col(fptr, id, const_cast<char*>(col.c_str()), otform, &status);
            fits::vif_check_cfitsio(status, "cannot re-insert column '"+col+"'");

            std::swap(naxes[0], naxes[1]);
            fits_write_tdim(fptr, id, 2, naxes.data(), &status);
            fits::vif_check_cfitsio(status, "cannot re-write TDIM for column '"+col+"'");

            fits_write_col(
                fptr, impl::fits_impl::traits<double>::ttype, id, 1, 1, n_elements(v),
                const_cast<vec2d::dtype*>(v.data.data()), &status
            );
            fits::vif_check_cfitsio(status, "cannot write transposed data for column '"+col+"'");
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_r2c_help() {
    using namespace terminal_format;

    paragraph("The program will convert a row-oriented FITS file to a column-oriented FITS file. "
        "No data is lost in the process.");

    header("List of available command line options:");
    bullet("out", "[string] path to the output file (has to be different than the input file)");
    bullet("silent", "[flag] do not print information on the conversion process");
    bullet("help", "[flag] print this text");
    print("");
}

bool rows_to_columns(int argc, char* argv[], const std::string& file) {
    std::string out;
    bool silent = false;
    bool help = false;
    read_args(argc, argv, arg_list(out, silent, help));

    if (help) {
        print_r2c_help();
        return true;
    }

    if (out == file) {
        error("output file cannot be equal to input file");
        return false;
    }

    file::mkdir(file::get_directory(out));

    if (!begins_with(out, "!")) {
        out = "!"+out;
    }

    fitsfile* ifptr = nullptr;
    fitsfile* ofptr = nullptr;
    int status = 0;

    try {
        fits_open_table(&ifptr, file.c_str(), READONLY, &status);
        fits::vif_check_cfitsio(status, "cannot open '"+file+"'");

        long nrow;
        fits_get_num_rows(ifptr, &nrow, &status);
        if (nrow == 1) {
            error("this FITS file is already column-oriented");
            fits_close_file(ifptr, &status);
            return 1;
        }

        int ncol;
        fits_get_num_cols(ifptr, &ncol, &status);
        if (!silent) print("col=", ncol, ", row=", nrow);

        fits_create_file(&ofptr, out.c_str(), &status);
        fits::vif_check_cfitsio(status, "cannot open '"+out+"'");
        fits_create_tbl(ofptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

        int oc = 1;
        for (int c = 1; c <= ncol; ++c) {
            // Read column info
            int type;
            long repeat, width;
            fits_get_eqcoltype(ifptr, c, &type, &repeat, &width, &status);

            // TODO: understand why TINT32BIT needs more than this function suggests
            if (type == TINT32BIT) {
                // 'fits_get_eqcoltype' suggests width=4, but actually width=8 is read
                // in 'fits_read_col'
                width *= 2;
            }

            char name[80] = {0};
            char comment[80];
            fits_read_key(ifptr, TSTRING, const_cast<char*>(("TTYPE"+to_string(c)).c_str()),
                name, comment, &status);
            char ittform[80] = {0};
            fits_read_key(ifptr, TSTRING, const_cast<char*>(("TFORM"+to_string(c)).c_str()),
                ittform, comment, &status);
            std::string itform = ittform;
            itform = trim(ittform);

            if (!silent) print(name);

            // Convert it to a column oriented format
            std::string tform, tdim;
            long nelem = 1;
            bool string = false;
            if (itform.size() == 1) {
                tform = to_string(nrow)+itform;
                tdim = "("+to_string(nrow)+")";
            } else {
                {
                    std::string tmp = itform; erase_end(tmp, 1);
                    from_string(tmp, nelem);
                }

                if (itform.back() == 'A') {
                    tform = to_string(nrow*nelem)+'A'+to_string(nelem);
                    string = true;
                } else {
                    tform = to_string(nrow*nelem)+itform.back();
                }

                tdim = "("+to_string(nelem)+","+to_string(nrow)+")";
            }

            // Create the new column
            fits_insert_col(ofptr, oc, name, const_cast<char*>(tform.c_str()), &status);
            fits::vif_check_cfitsio(status, "error writing '"+std::string(name)+"'");
            char tdim_com[] = "size of the multidimensional array";
            fits_write_key(ofptr, TSTRING, const_cast<char*>(("TDIM"+to_string(c)).c_str()),
                const_cast<char*>(tdim.c_str()), tdim_com, &status);

            if (string) {
                // Read the input data
                char nul[] = "?";
                int nnul = 0;

                std::vector<char> tdata(nrow*(nelem+1));
                char** data = new char*[nrow];
                for (long i : range(nrow)) {
                    data[i] = tdata.data() + i*(nelem+1);
                }

                fits_read_col_str(ifptr, c, 1, 1, nrow, nul, data, &nnul, &status);
                fits::vif_check_cfitsio(status, "error reading '"+std::string(name)+"'");

                // Write the output data
                fits_write_col(ofptr, impl::fits_impl::traits<std::string>::ttype, oc, 1, 1, nrow, data, &status);
                fits::vif_check_cfitsio(status, "error writing '"+std::string(name)+"'");

                delete[] data;
            } else {
                // Read the input data
                long nul = 0;
                int nnul = 0;
                std::vector<char> data(width*nelem*nrow);
                fits_read_col(ifptr, type, c, 1, 1, nelem*nrow, &nul, data.data(), &nnul, &status);
                fits::vif_check_cfitsio(status, "error reading '"+std::string(name)+"'");

                // Write the output data
                fits_write_col(ofptr, type, oc, 1, 1, nelem*nrow, data.data(), &status);
                fits::vif_check_cfitsio(status, "error writing '"+std::string(name)+"'");
            }

            ++oc;
        }

        fits_close_file(ifptr, &status);
        fits_close_file(ofptr, &status);
    } catch (fits::exception& e) {
        if (ifptr) fits_close_file(ifptr, &status);
        if (ofptr) fits_close_file(ofptr, &status);
        print(e.msg);
        return false;
    }

    return true;
}

void print_cpwcs_help() {
    using namespace terminal_format;

    paragraph("The program will copy all the WCS related header keywords from the provided file to "
        "another file. Previous WCS data in the other file will be lost.");

    header("List of available command line options:");
    bullet("out", "[string] path to the file to copy the keywords to");
    bullet("help", "[flag] print this text");
    print("");
}

bool copy_wcs_header(int argc, char* argv[], const std::string& file) {
    std::string out;
    bool help = false;
    read_args(argc, argv, arg_list(out, help));

    if (help) {
        print_cpwcs_help();
        return true;
    }

    if (out.empty()) {
        error("missing output file, please provide 'out=...'");
        return false;
    }

    // List of keywords taken from 'cphead' (WCSTools)
    vec1s keywords = {"RA", "DEC", "EPOCH", "EQUINOX", "RADECSYS", "SECPIX", "IMWCS",
        "PC001001", "PC001002", "PC002001", "PC002002", "LATPOLE", "LONPOLE"};

    for (uint_t i = 1; i < 13; ++i) {
        keywords.push_back("CO1_"+to_string(i));
        keywords.push_back("CO2_"+to_string(i));
    }
    for (uint_t i = 0; i < 10; ++i) {
        keywords.push_back("PROJP"+to_string(i));
        keywords.push_back("PV1_"+to_string(i));
        keywords.push_back("PV2_"+to_string(i));
        keywords.push_back("CTYPE"+to_string(i));
        keywords.push_back("CUNIT"+to_string(i));
        keywords.push_back("CRVAL"+to_string(i));
        keywords.push_back("CRPIX"+to_string(i));
        keywords.push_back("CDELT"+to_string(i));
        keywords.push_back("CROTA"+to_string(i));
        keywords.push_back("SECPIX"+to_string(i));
        keywords.push_back("CD"+to_string(i));
        for (uint_t j = 0; j < 10; ++j) {
            keywords.push_back("CD"+to_string(i)+"_"+to_string(j));
            keywords.push_back("PC"+to_string(i)+"_"+to_string(j));
        }
    }

    fitsfile* fptr1;
    fitsfile* fptr2;
    int status = 0;

    try {
        fits_open_image(&fptr1, file.c_str(), READONLY, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");
        fits_open_image(&fptr2, out.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+out+"'");

        for (auto& k : keywords) {
            char value[80] = {0};
            char comment[80] = {0};
            fits_read_keyword(fptr1, const_cast<char*>(k.c_str()), value, comment, &status);
            if (status != 0) {
                status = 0;
                continue;
            }

            if (value[0] == '\'') {
                for (uint_t i = 79; i != npos; --i) {
                    if (value[i] == '\'') {
                        value[i] = '\0';
                        break;
                    }
                }

                fits_update_key(fptr2, TSTRING, const_cast<char*>(k.c_str()), value+1, comment,
                    &status);
            } else {
                double d;
                from_string(value, d);
                fits_update_key(fptr2, TDOUBLE, const_cast<char*>(k.c_str()), &d, comment, &status);
            }
        }

        fits_close_file(fptr1, &status);
        fits_close_file(fptr2, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr1) fits_close_file(fptr1, &status);
        if (fptr2) fits_close_file(fptr2, &status);
        return false;
    }

    return true;
}

void print_hdr_help() {
    using namespace terminal_format;

    paragraph("The program will just display the content of the provided FITS file's header.");

    header("List of available command line options:");
    bullet("raw", "[flag] do not attempt to pretty print the header, just dump it as it is");
    bullet("verbose", "[flag] print additional information");
    bullet("help", "[flag] print this text");
    print("");
}

bool show_header(int argc, char* argv[], const std::string& file) {
    bool pretty = false;
    bool help = false;
    bool verbose = false;
    uint_t hdu = npos;
    read_args(argc, argv, arg_list(pretty, hdu, help, verbose));

    if (help) {
        print_remove_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_file(&fptr, file.c_str(), READONLY, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        if (hdu != npos) {
            int nhdu = 0;
            fits_get_num_hdus(fptr, &nhdu, &status);
            vif_check(hdu < uint_t(nhdu), "requested HDU does not exists in this "
                "FITS file (", hdu, " vs. ", nhdu, ")");

            fits_movabs_hdu(fptr, hdu+1, nullptr, &status);
        }

        if (!pretty) {
            // Read the header as a string
            char* hstr = nullptr;
            int nkeys  = 0;
            fits_hdr2str(fptr, 0, nullptr, 0, &hstr, &nkeys, &status);
            std::string header = hstr;
            free(hstr);

            for (uint_t i = 0; i < header.size(); i += 80) {
                print(header.substr(i, std::min(header.size()-i, std::size_t(80))));
            }
        } else {
            print("WIP");
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_editkwd_help() {
    using namespace terminal_format;

    paragraph("The program will allow you to edit manually any keyword within this FITS file. It "
        "is an interactive tool, that first asks your for the name of the keyword, then for the "
        "new value");

    header("List of available command line options:");
    bullet("verbose", "[flag] print additional information");
    bullet("help", "[flag] print this text");
    print("");
}

bool write_keyword_value(fitsfile* fptr, bool newk, std::string name, std::string new_value,
    char* comment) {

    int status = 0;
    double d;
    bool is_number = from_string(new_value, d);

    if (!is_number) {
        if (new_value[0] == '\'') {
            new_value.erase(0,1);
            for (uint_t i = new_value.size()-1; i != npos; --i) {
                if (new_value[i] == '\'') {
                    new_value.resize(i);
                    return false;
                }
            }
        }

        if (newk) {
            fits_write_key(fptr, TSTRING, const_cast<char*>(name.c_str()),
                const_cast<char*>(new_value.c_str()), comment, &status);
            fits::vif_check_cfitsio(status, "cannot write keyword '"+name+"'");
        } else {
            fits_update_key(fptr, TSTRING, const_cast<char*>(name.c_str()),
                const_cast<char*>(new_value.c_str()), comment, &status);
            fits::vif_check_cfitsio(status, "cannot update keyword '"+name+"'");
        }
    } else {
        if (ends_with(new_value, ".") || find(new_value, ".") == npos) {
            int di = d;
            if (newk) {
                fits_write_key(fptr, TINT, const_cast<char*>(name.c_str()), &di, comment,
                    &status);
                fits::vif_check_cfitsio(status, "cannot write keyword '"+name+"'");
            } else {
                fits_update_key(fptr, TINT, const_cast<char*>(name.c_str()), &di, comment,
                    &status);
                fits::vif_check_cfitsio(status, "cannot update keyword '"+name+"'");
            }
        } else {
            if (newk) {
                fits_write_key(fptr, TDOUBLE, const_cast<char*>(name.c_str()), &d, comment,
                    &status);
                fits::vif_check_cfitsio(status, "cannot write keyword '"+name+"'");
            } else {
                fits_update_key(fptr, TDOUBLE, const_cast<char*>(name.c_str()), &d, comment,
                    &status);
                fits::vif_check_cfitsio(status, "cannot update keyword '"+name+"'");
            }
        }
    }

    return true;
}

bool edit_keyword(int argc, char* argv[], const std::string& file) {
    bool help = false;
    bool verbose = false;
    uint_t hdu = npos;
    std::string key, val;
    read_args(argc, argv, arg_list(help, verbose, hdu, key, val));

    if (help) {
        print_remove_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_image(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        int naxis = 0;
        fits_get_img_dim(fptr, &naxis, &status);
        if (naxis == 0) {
            fits_close_file(fptr, &status);
            fits_open_table(&fptr, file.c_str(), READWRITE, &status);
            fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");
            if (verbose) print("loaded table file");
        } else {
            if (verbose) print("loaded image file");
        }

        if (hdu != npos) {
            int nhdu = 0;
            fits_get_num_hdus(fptr, &nhdu, &status);
            vif_check(hdu < uint_t(nhdu), "requested HDU does not exists in this "
                "FITS file (", hdu, " vs. ", nhdu, ")");

            fits_movabs_hdu(fptr, hdu+1, nullptr, &status);
        }

        if (key.empty() && val.empty()) {
            while (true) {
                print("please type the name of the keyword you want to edit (empty to exit):");
                std::string name;
                std::cout << "> ";
                std::getline(std::cin, name);
                name = to_upper(trim(name));
                if (name.empty()) break;

                char value[80];
                char comment[80];
                bool newk = false;
                fits_read_keyword(fptr, const_cast<char*>(name.c_str()), value, comment, &status);
                if (status != 0) {
                    status = 0;
                    newk = true;
                    std::cout << "note: this keyword does not exist, do you want to create it? (y/n) "
                        << std::flush;

                    std::string line;
                    bool stop = false;
                    while (line.empty()) {
                        std::getline(std::cin, line);
                        line = to_lower(line);
                        if (line == "y" || line == "yes") break;
                        if (line == "n" || line == "no") {
                            stop = true;
                            break;
                        }

                        line.clear();
                        std::cout << "error: please answer either 'yes' (y) or 'no' (n): "
                            << std::flush;
                    }

                    if (stop) {
                        continue;
                    }

                    print("");
                }

                if (!newk) {
                    print(name, " = ", value, " (", comment, ")");
                }

                print("please type the new value for this keyword (empty to abort):");
                std::string new_value;
                std::cout << "> ";
                std::getline(std::cin, new_value);
                new_value = trim(new_value);
                if (new_value.empty()) break;

                if (!write_keyword_value(fptr, newk, name, new_value, comment)) {
                    error("could not write keyword");
                    break;
                }

                fits_read_keyword(fptr, const_cast<char*>(name.c_str()), value, comment, &status);
                print(name, " = ", value, " (", comment, ")");
            }
        } else {
            char value[80];
            char comment[80];
            bool newk = false;
            fits_read_keyword(fptr, const_cast<char*>(key.c_str()), value, comment, &status);
            if (status != 0) {
                status = 0;
                newk = true;
            }

            if (!write_keyword_value(fptr, newk, key, val, comment)) {
                error("could not write keyword");
            }
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_rmkwd_help() {
    using namespace terminal_format;

    paragraph("The program will remove the listed keywords (if found) from the provided FITS "
        "file. This is a destructive procedure: the FITS file will be modified and data will "
        "be lost.");

    header("List of available command line options:");
    bullet("key", "[string array] list of keywords to remove");
    bullet("verbose", "[flag] print additional information");
    bullet("help", "[flag] print this text");
    print("");
}

bool remove_keyword(int argc, char* argv[], const std::string& file) {
    vec1s key;
    bool help = false;
    bool verbose = false;
    read_args(argc, argv, arg_list(key, help, verbose));

    if (help) {
        print_rmkwd_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_image(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        int naxis = 0;
        fits_get_img_dim(fptr, &naxis, &status);
        if (naxis == 0) {
            fits_close_file(fptr, &status);
            fits_open_table(&fptr, file.c_str(), READWRITE, &status);
            fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");
            if (verbose) print("loaded table file");
        } else {
            if (verbose) print("loaded image file");
        }

        for (auto& k : key) {
            while (status == 0) {
                fits_delete_key(fptr, k.c_str(), &status);
            }
            status = 0;
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_readkwd_help() {
    using namespace terminal_format;

    paragraph("The program will print the values of the provided keywords from the provided FITS "
        "file. Keywords that do not exist will print empty values.");

    header("List of available command line options:");
    bullet("key", "[string array] list of keywords to remove");
    bullet("verbose", "[flag] print additional information");
    bullet("help", "[flag] print this text");
    print("");
}

bool read_keyword(int argc, char* argv[], const std::string& file) {
    vec1s key;
    uint_t hdu = npos;
    bool help = false;
    bool verbose = false;
    read_args(argc, argv, arg_list(key, help, verbose, hdu));

    if (help) {
        print_readkwd_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_image(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        int naxis = 0;
        fits_get_img_dim(fptr, &naxis, &status);
        if (naxis == 0) {
            fits_close_file(fptr, &status);
            fits_open_table(&fptr, file.c_str(), READWRITE, &status);
            fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");
            if (verbose) print("loaded table file");
        } else {
            if (verbose) print("loaded image file");
        }

        if (hdu != npos) {
            int nhdu = 0;
            fits_get_num_hdus(fptr, &nhdu, &status);
            vif_check(hdu < uint_t(nhdu), "requested HDU does not exists in this "
                "FITS file (", hdu, " vs. ", nhdu, ")");

            fits_movabs_hdu(fptr, hdu+1, nullptr, &status);
        }

        for (auto& k : key) {
            status = 0;
            char val[80];
            fits_read_keyword(fptr, k.c_str(), val, nullptr, &status);
            if (status == 0) {
                if (val[0] == '\'') {
                    // String, remove quotes and trim spaces
                    std::string sval = val;
                    auto p0 = sval.find_first_of('\'');
                    auto p1 = sval.find_last_of('\'');
                    if (p0 == p1) {
                        print("");
                    } else {
                        print(trim(sval.substr(p0+1, p1-p0-1)));
                    }
                } else {
                    print(val);
                }
            } else {
                print("");
            }
        }

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_make2d_help() {
    using namespace terminal_format;

    paragraph("The program will remove and modify specific keywords to transform a "
        "multidimensional image into a simple 2D image. Note that the pixel content "
        "will not be changed! This is just keyword manipulation. This procedure can "
        "only run on images for which the extra dimensions have size 1.");

    header("List of available command line options:");
    bullet("help", "[flag] print this text");
    print("");
}

bool make_2d(int argc, char* argv[], const std::string& file) {
    bool help = false;
    read_args(argc, argv, arg_list(help));

    if (help) {
        print_make2d_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    int status = 0;

    try {
        fits_open_image(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        int naxis = 0;
        fits_get_img_dim(fptr, &naxis, &status);
        if (naxis < 2) {
            error("'make2d' requires an image FITS file");
            return false;
        } else if (naxis == 2) {
            // Nothing to do
            return true;
        }

        // Check the dimensions are of size 1
        std::vector<long> naxes(naxis);
        fits_get_img_size(fptr, naxis, naxes.data(), &status);
        for (int i = 2; i < naxis; ++i) {
            if (naxes[i] != 1) {
                error("cannot remove extra dimensions from '"+file+"'");
                note("these dimensions contain data");
                return false;
            }
        }

        // First remove the extra keywords
        vec1s keybase = {"NAXIS", "CTYPE", "CUNIT", "CRVAL", "CRPIX", "CROTA", "CDELT"};

        for (int i = 3; i <= naxis; ++i) {
            for (auto& k : keybase) {
                fits_delete_key(fptr, (k+to_string(i)).c_str(), &status);
                if (status != 0) {
                    status = 0;
                }
            }

            for (int j = 1; j <= naxis; ++j) {
                fits_delete_key(fptr, ("PC"+align_right(to_string(i),2,'0')+"_"+
                                            align_right(to_string(j),2,'0')).c_str(), &status);
                if (status != 0) {
                    status = 0;
                }
                fits_delete_key(fptr, ("PC"+align_right(to_string(j),2,'0')+"_"+
                                            align_right(to_string(i),2,'0')).c_str(), &status);
                if (status != 0) {
                    status = 0;
                }
            }
        }

        // Then modify NAXIS
        naxis = 2;
        char comment[] = "";
        fits_update_key(fptr, TINT, "NAXIS", &naxis, comment, &status);

        fits_close_file(fptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        return false;
    }

    return true;
}

void print_meta_help() {
    using namespace terminal_format;

    paragraph("The program will identify 'meta' columns in a column oriented FITS table. It will "
        "then move these columns to a new table within the FITS file, so that the file can be "
        "open in programs that do not support meta columns (such as Topcat). No data is lost in "
        "the process.");

    header("List of available command line options:");
    bullet("force", "[flag] do not ask for confirmation at each step");
    bullet("help", "[flag] print this text");
    print("");
}

bool move_meta_columns(int argc, char* argv[], const std::string& file) {
    bool force = false;
    bool help = false;
    read_args(argc, argv, arg_list(force, help));

    if (help) {
        print_meta_help();
        return true;
    }

    fitsfile* fptr = nullptr;
    fitsfile* ofptr = nullptr;
    int status = 0;

    try {
        fits_open_table(&fptr, file.c_str(), READWRITE, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

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

        vec1u udim = unique_values(dims);
        uint_t dim;
        if (udim.size() == 1) {
            dim = udim[0];
        } else {
            vec1u count(udim.size());
            for (uint_t i = 0; i < udim.size(); ++i) {
                count[i] = total(dims == udim[i]);
            }

            vec1u weight = indgen(udim.size()) + sort(count);
            dim = udim[max_id(weight)];

            if (!force) {
                std::cout << "identified 'row' dimension: " << dim << "\n";
                std::cout << "is it correct? (y/n) " << std::flush;
                std::string line;
                while (line.empty()) {
                    std::getline(std::cin, line);
                    line = to_lower(line);
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
            return true;
        }

        if (!force) {
            vec1s names;
            for (auto i : mcols) {
                std::string id = to_string(i);
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
                line = to_lower(line);
                if (line == "y" || line == "yes") break;
                if (line == "n" || line == "no") {
                    fits_close_file(fptr, &status);
                    return true;
                }
                line.clear();
                std::cout << "error: please answer either 'yes' (y) or 'no' (n): " << std::flush;
            }
        }

        fits_open_table(&ofptr, file.c_str(), READWRITE, &status);

        fits_create_tbl(fptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

        for (auto i : mcols) {
            fits_copy_col(ofptr, fptr, i, mcols.size(), 1, &status);
            fits::vif_check_cfitsio(status, "cannot copy column "+to_string(i));
        }

        for (auto i : reverse(mcols)) {
            fits_delete_col(ofptr, i, &status);
        }

        fits_close_file(fptr, &status);
        fits_close_file(ofptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (fptr) fits_close_file(fptr, &status);
        if (ofptr) fits_close_file(ofptr, &status);
        return false;
    }

    return true;
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
                    line = to_lower(line);
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

void print_extract_help() {
    using namespace terminal_format;

    paragraph("The program will extract one extension from the provided FITS file and save it "
        "(with its header and content) in a new FITS file.");

    header("List of available command line options:");
    bullet("out", "[string] path to the output file, must be provided");
    bullet("hdu", "[integer] HDU ID of the extension to copy (default: 0, the primary)");
    bullet("help", "[flag] print this text");
    print("");
}

bool extract_extension(int argc, char* argv[], const std::string& file) {
    std::string out;
    bool help = false;
    uint_t hdu = 0;
    read_args(argc, argv, arg_list(out, hdu, help));

    if (help) {
        print_extract_help();
        return true;
    }

    if (out.empty()) {
        error("missing output file name (out=...)");
        print_extract_help();
        return false;
    }

    file::mkdir(file::get_directory(out));

    if (!begins_with(out, "!")) {
        out = "!"+out;
    }

    fitsfile* ifptr = nullptr;
    fitsfile* ofptr = nullptr;
    int status = 0;

    try {
        fits_open_file(&ifptr, file.c_str(), READONLY, &status);
        fits::vif_check_cfitsio(status, "cannot open file '"+file+"'");

        int hdut = 0;
        fits_movabs_hdu(ifptr, hdu+1, &hdut, &status);
        fits::vif_check_cfitsio(status, "cannot reach HDU "+to_string(hdu));

        fits_create_file(&ofptr, out.c_str(), &status);
        fits::vif_check_cfitsio(status, "cannot create file '"+out+"'");

        fits_copy_hdu(ifptr, ofptr, 0, &status);

        fits_close_file(ifptr, &status);
        fits_close_file(ofptr, &status);
    } catch (fits::exception& e) {
        print(e.msg);
        if (ofptr) fits_close_file(ofptr, &status);
        if (ifptr) fits_close_file(ifptr, &status);
        return false;
    }

    return true;
}

void print_help() {
    using namespace terminal_format;

    print("fitstool v1.0");
    header("Usage: fitstool operation [options]");
    header("Available operations:");
    bullet("remove", "remove the specified columns from this FITS file");
    bullet("transpose", "transpose the specified columns within this FITS file");
    bullet("r2c", "convert a row oriented FITS table to a column oriented one");
    bullet("cpwcs", "copy the WCS keywords from this FIST file into another");
    bullet("hdr", "print the header of the provided FITS file");
    bullet("editkwd", "modify or add keywords to the provided FIST file interactively");
    bullet("rmkwd", "remove keywords from the provided FIST file");
    bullet("readkwd", "print keywords from the provided FIST file");
    bullet("make2d", "remove extra dimensions from FITS image");
    bullet("meta", "shift all 'meta' columns into a new extension so that the file can be read in "
        "programs like TOPCAT that expect a single invariant dimension for all columns");
    bullet("extract", "extracts one extension (or HDU) from a FITS file and save it in another");
    print("");
    header("To learn more about each operations and see the list of avilable options, run "
        "'fitstool operation help'.");

    print("");
    paragraph("Copyright (c) 2014 C. Schreiber (corentin.schreiber@cea.fr)");

    paragraph("This software is provided 'as-is', without any express or implied warranty. In no "
        "event will the authors be held liable for any damages arising from the use of this "
        "software.");

    paragraph("Permission is granted to anyone to use this software for any purpose, including "
        "commercial applications, and to alter it and redistribute it freely, subject to the "
        "following restrictions:");

    bullet("1", "The origin of this software must not be misrepresented; you must not claim that "
        "you wrote the original software. If you use this software in a product, an acknowledgment "
        "in the product documentation would be appreciated but is not required.");
    bullet("2", "Altered source versions must be plainly marked as such, and must not be "
        "misrepresented as being the original software.");
    bullet("3", "This notice may not be removed or altered from any source distribution.");

    print("");
}
