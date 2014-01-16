#include <phypp.hpp>

void print_help();

int main(int argc, char* argv[]) {
    std::string in, out;
    bool silent = false;

    if (argc < 3) {
        print_help();
        return 1;
    }

    in = argv[1];
    out = argv[2];

    if (in == out) {
        error("output file cannot be equal to input file");
        return 1;
    }

    read_args(argc-2, argv+2, arg_list(silent));

    file::mkdir(file::get_directory(out));

    if (!start_with(out, "!")) {
        out = "!"+out;
    }

    fitsfile* ifptr = nullptr;
    fitsfile* ofptr = nullptr;
    int status = 0;

    try {
        fits_open_table(&ifptr, in.c_str(), READONLY, &status);
        fits::phypp_check_cfitsio(status, "cannot open '"+in+"'");

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
        fits::phypp_check_cfitsio(status, "cannot open '"+out+"'");
        fits_create_tbl(ofptr, BINARY_TBL, 1, 0, 0, 0, 0, nullptr, &status);

        int oc = 1;
        for (int c = 1; c <= ncol; ++c) {
            // Read column info
            int type;
            long repeat, width;
            fits_get_coltype(ifptr, c, &type, &repeat, &width, &status);
            char name[80] = {0};
            char comment[80];
            fits_read_key(ifptr, TSTRING, const_cast<char*>(("TTYPE"+strn(c)).c_str()),
                name, comment, &status);
            char ittform[80] = {0};
            fits_read_key(ifptr, TSTRING, const_cast<char*>(("TFORM"+strn(c)).c_str()),
                ittform, comment, &status);
            std::string itform = ittform;
            itform = trim(ittform);

            if (!silent) print(name);

            // Convert it to a column oriented format
            std::string tform, tdim;
            long nelem = 1;
            int otype = type;
            bool string = false;
            if (itform.size() == 1) {
                tform = strn(nrow)+itform;
                tdim = "("+strn(nrow)+")";
            } else {
                {
                    std::string tmp = itform; erase_end(tmp, 1);
                    from_string(tmp, nelem);
                }

                if (itform.back() == 'A') {
                    tform = strn(nrow*nelem)+'B';
                    string = true;
                } else {
                    tform = strn(nrow*nelem)+itform.back();
                }

                tdim = "("+strn(nelem)+","+strn(nrow)+")";
            }

            if (string) {
                // Read the input data
                char nul[] = "?";
                int nnul = 0;

                std::vector<char> tdata(nrow*(nelem+1));
                char** idata = new char*[nrow];
                for (long i = 0; i < nrow; ++i) {
                    idata[i] = tdata.data() + i*(nelem+1);
                }

                fits_read_col_str(ifptr, c, 1, 1, nrow, nul, idata, &nnul, &status);
                fits::phypp_check_cfitsio(status, "error reading '"+std::string(name)+"'");

                delete[] idata;

                std::vector<char> data(nrow*nelem);
                for (long i = 0; i < nrow; ++i) {
                    for (long j = 0; j < nelem; ++j) {
                        data[i*nelem+j] = tdata[i*(nelem+1)+j];
                    }
                }

                // Write the output data
                fits_insert_col(ofptr, oc, name, const_cast<char*>(tform.c_str()), &status);
                fits::phypp_check_cfitsio(status, "error writing '"+std::string(name)+"'");
                char tdim_com[] = "size of the multidimensional array";
                fits_write_key(ofptr, TSTRING, const_cast<char*>(("TDIM"+strn(c)).c_str()),
                    const_cast<char*>(tdim.c_str()), tdim_com, &status);
                fits_write_col(ofptr, TBYTE, oc, 1, 1, nelem*nrow, data.data(), &status);
                fits::phypp_check_cfitsio(status, "error writing '"+std::string(name)+"'");
            } else {
                // Read the input data
                long nul = 0;
                int nnul = 0;
                std::vector<char> data(width*nelem*nrow);
                fits_read_col(ifptr, type, c, 1, 1, nelem*nrow, &nul, data.data(), &nnul, &status);
                fits::phypp_check_cfitsio(status, "error reading '"+std::string(name)+"'");

                // Write the output data
                fits_insert_col(ofptr, oc, name, const_cast<char*>(tform.c_str()), &status);
                fits::phypp_check_cfitsio(status, "error writing '"+std::string(name)+"'");
                char tdim_com[] = "size of the multidimensional array";
                fits_write_key(ofptr, TSTRING, const_cast<char*>(("TDIM"+strn(c)).c_str()),
                    const_cast<char*>(tdim.c_str()), tdim_com, &status);
                fits_write_col(ofptr, type, oc, 1, 1, nelem*nrow, data.data(), &status);
                fits::phypp_check_cfitsio(status, "error writing '"+std::string(name)+"'");
            }

            ++oc;
        }

        fits_close_file(ifptr, &status);
        fits_close_file(ofptr, &status);
    } catch (fits::exception& e) {
        if (ifptr) fits_close_file(ifptr, &status);
        if (ofptr) fits_close_file(ofptr, &status);
        error(e.msg);
        throw;
    }

    return 0;
}

void print_help() {
    using namespace format;

    print("fitsr2c v1.0");
    paragraph("usage: fitsr2c input output [options=...]");

    paragraph(
        "The program will convert a row-oriented FITS file to a column-oriented FITS file. No data "
        "is lost in the process."
    );

    header("List of available command line options:");
    bullet("silent", "[flag] do not print information on the conversion process");
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
