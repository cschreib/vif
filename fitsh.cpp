#include <CCfits/CCfits>
#include <iostream>
#include <algorithm>

namespace fits {
    using namespace CCfits;
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "fitsh v1.0" << std::endl;
        std::cout << "usage: fitsh [file]" << std::endl;
        return 0;
    }
    
    // Open the file
    fits::FITS f(argv[1], fits::Read);
    
    // Get the first HDU (extension)
    fits::PHDU& hdu = f.pHDU();
    
    fitsfile* fp = hdu.fitsPointer();
    
    // Read the header as a string
    char* hstr = nullptr;
    int nkeys  = 0;
    int status = 0;
    fits_hdr2str(fp, 0, nullptr, 0, &hstr, &nkeys, &status);
    
    std::string header = hstr;
    free(hstr);
    
    for (std::size_t i = 0; i < header.size(); i += 80) {
        std::cout << header.substr(i, std::min(header.size()-i, std::size_t(80))) << std::endl;
    }
    
    return 0;
}
