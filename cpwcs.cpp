#include <phypp.hpp>

int main(int argc, char* argv[]) {
    if (argc < 3) return 0;

    std::string file1 = argv[1];
    std::string file2 = argv[2];

    // List of keywords taken from 'cphead' from WCSTools
    vec1s keywords = {"RA", "DEC", "EPOCH", "EQUINOX", "RADECSYS", "SECPIX", "SECPIX1", "SECPIX2",
        "CTYPE1", "CTYPE2", "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CRPIX1", "CRPIX2", "CROTA1",
        "CROTA2", "IMWCS", "CD1_1", "CD1_2", "CD2_1", "CD2_2", "PC1_1", "PC1_2", "PC2_1", "PC2_2",
        "PC001001", "PC001002", "PC002001", "PC002002", "LATPOLE", "LONPOLE"};

    for (uint_t i = 1; i < 13; ++i) {
        keywords.push_back("CO1_"+strn(i));
    }
    for (uint_t i = 1; i < 13; ++i) {
        keywords.push_back("CO2_"+strn(i));
    }
    for (uint_t i = 0; i < 10; ++i) {
        keywords.push_back("PROJP"+strn(i));
    }
    for (uint_t i = 0; i < 10; ++i) {
        keywords.push_back("PV1_"+strn(i));
    }
    for (uint_t i = 0; i < 10; ++i) {
        keywords.push_back("PV2_"+strn(i));
    }

    fitsfile* fptr1;
    fitsfile* fptr2;
    int status = 0;

    fits_open_image(&fptr1, file1.c_str(), READONLY, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+file1+"'");
    fits_open_image(&fptr2, file2.c_str(), READWRITE, &status);
    fits::phypp_check_cfitsio(status, "cannot open file '"+file2+"'");

    for (auto& s : keywords) {
        char value[80] = {0};
        char comment[80] = {0};
        fits_read_keyword(fptr1, const_cast<char*>(s.c_str()), value, comment, &status);
        if (status != 0) {
            status = 0;
            continue;
        }

        if (value[0] == '\'') {
            for (uint_t i = 80; i != npos; --i) {
                if (value[i] == '\'') {
                    value[i] = '\0';
                    break;
                }
            }

            fits_update_key(fptr2, TSTRING, const_cast<char*>(s.c_str()), value+1, comment, &status);
        } else {
            double d;
            from_string(value, d);
            fits_update_key(fptr2, TDOUBLE, const_cast<char*>(s.c_str()), &d, comment, &status);
        }
    }

    fits_close_file(fptr1, &status);
    fits_close_file(fptr2, &status);

    return 0;
}

