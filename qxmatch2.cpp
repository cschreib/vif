#include <phypp.hpp>

void print_help() {
    std::cout << "qxmatch2 v2.0\n";
    std::cout << "  usage: qxmatch2 cats={file1,file2} output=ofile [options=...]\n\n";
    
    std::cout << "  'file1' and 'file2' must be FITS files, expected to contain an extension with\n";
    std::cout << "  at least two vector columns named 'RA' and 'DEC', each containing a single row.\n";
    std::cout << "  Both are assumed to contain galactic coordinates in degrees.\n";
    std::cout << "  They can be created in IDL using:\n";
    std::cout << "          mwrfits, {ra:[...], dec:[...]}, \"filename.fits\", /create\n";
    std::cout << "  The output file is written in a similar way: it contains one extension named\n";
    std::cout << "  'RESULT_BINARY', which contains a single column named 'ID', itself being a two\n";
    std::cout << "  dimensional vector containing the ID of the n'th closest sources. It can be \n";
    std::cout << "  read in IDL using:\n";
    std::cout << "          res = mrdfits(\"output.fits\", 1)\n";
    std::cout << "          best = res.id[*,0]\n";
    std::cout << "          second_best = res.id[*,1]\n";
    std::cout << "  The program also outputs the crossmatch distance in a column called 'D', with a\n";
    std::cout << "  format similar to that of 'ID', in arcseconds.\n\n";
    
    std::cout << "  It is also possible to compute an 'auto crossmatch' of a single catalog by only\n";
    std::cout << "  providing one file in the 'cats' parameter. The program will then compute the\n";
    std::cout << "  n'th nearest neighbors for each source within this catalog.\n";

    std::cout << "  List of available command line options:\n";
    std::cout << "    verbose: set this flag to print additional information in the standard output\n";
    std::cout << "    nth=[number]: set this value to the number of closest neighbors you want to\n";
    std::cout << "                  retrieve (default: 1).\n";
    std::cout << "    thread=[number]: set this value to the number of concurrent threads you want\n";
    std::cout << "                     to run (default: 1).\n";
    
    std::cout << "  Copyright (c) 2013 C. Schreiber (corentin.schreiber@cea.fr)\n\n";

    std::cout << "  This software is provided 'as-is', without any express or implied warranty.\n";
    std::cout << "  In no event will the authors be held liable for any damages arising from the\n";
    std::cout << "  use of this software.\n";
      
    std::cout << "  Permission is granted to anyone to use this software for any purpose,\n";
    std::cout << "  including commercial applications, and to alter it and redistribute it\n";
    std::cout << "  freely, subject to the following restrictions:\n\n";
      
    std::cout << "    1. The origin of this software must not be misrepresented; you must not\n";
    std::cout << "       claim that you wrote the original software. If you use this software in\n";
    std::cout << "       a product, an acknowledgment in the product documentation would be\n";
    std::cout << "       appreciated but is not required.\n\n";
                  
    std::cout << "    2. Altered source versions must be plainly marked as such, and must not be\n";
    std::cout << "       misrepresented as being the original software.\n\n";
                  
    std::cout << "    3. This notice may not be removed or altered from any source distribution.\n\n";
    
    std::cout << std::flush;
}

int main(int argc, char* argv[]) {
    vec1s cats;
    std::string output; 

    uint_t nth = 1;
    uint_t thread = 1;
    bool   verbose = false;
    bool   quiet = false;

    read_args(argc, argv, arg_list(cats, output, nth, thread, verbose, quiet));

    if (quiet) verbose = false;

    if (output.empty()) {
        if (!quiet) print_help();
        return 0;
    }

    struct coord_t {
        vec1d ra, dec;
    };

    qxmatch_res res;

    if (n_elements(cats) == 2) {
        coord_t cat1, cat2;
        fits::read_table(cats[0], cat1);
        fits::read_table(cats[1], cat2);

        if (verbose) {
            print("qxmatch: crossmatching ", n_elements(cat1.ra), " sources from '", cats[0], 
                "' with ", n_elements(cat2.ra), " sources from '", cats[1], "'...");
        }

        res = qxmatch(cat1, cat2, keywords(_nth(nth), _thread(thread), _verbose(verbose)));
    } else if (n_elements(cats) == 1) {
        coord_t cat;
        fits::read_table(cats[0], cat);

        if (verbose) {
            print("qxmatch: self matching ", n_elements(cat.ra), " sources from '", cats[0], "'...");
        }

        res = qxmatch(cat, cat, keywords(_nth(nth), _thread(thread), _verbose(verbose), _self(true)));
    } else {
        if (!quiet) print_help();
        return 0;
    }

    qxmatch_save(output, res);

    return 0;
}
