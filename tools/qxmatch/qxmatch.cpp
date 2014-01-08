//  Copyright (c) 2013 C. Schreiber
//  
//  This software is provided 'as-is', without any express or implied warranty. In no event will the
//  authors be held liable for any damages arising from the use of this software.
//  
//  Permission is granted to anyone to use this software for any purpose, including commercial
//  applications, and to alter it and redistribute it freely, subject to the following restrictions:
//  
//      1. The origin of this software must not be misrepresented; you must not claim that you wrote the
//         original software. If you use this software in a product, an acknowledgment in the product
//         documentation would be appreciated but is not required.
//  
//      2. Altered source versions must be plainly marked as such, and must not be misrepresented as
//         being the original software.
//  
//      3. This notice may not be removed or altered from any source distribution.
//

#include <CCfits/CCfits>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <limits>
#include <algorithm>
#include <chrono>
#include <thread>
#include <atomic>

std::string strtrim(const std::string& str) {
    if (str.empty()) return str;
    
    std::size_t is = 0;
    while (is != str.size() && (str[is] == ' ' || str[is] == '\t')) ++is;
    std::size_t ie = str.size()-1;
    while (ie != 0          && (str[ie] == ' ' || str[ie] == '\t')) --ie;
    return str.substr(is, ie-is+1);
}

template<typename T>
T from_string(const std::string& str) {
    std::istringstream ss(str);
    T t;
    ss >> t;
    return t;
}

template<typename T>
bool from_string(T& t, const std::string& str) {
    std::istringstream ss(str);
    return ss >> t;
}

template<typename T>
std::string to_string(const T& t) {
    std::ostringstream ss;
    ss << t;
    return ss.str();
}

template<typename T>
std::string to_string(const T& t, std::size_t n) {
    if (n <= 1) {
        return to_string(t);
    }
    if (t == 0) {
        return std::string(floor(log10(n-1))+1, '0');
    } else {
        std::ostringstream ss;
        ss << t;
        std::size_t nz = floor(log10(n-1)) - floor(log10(t));
        if (nz > 0 && nz < 6) {
            return std::string(nz, '0') + ss.str();
        } else {
            return ss.str();    
        }
    }
}

std::string date_str(double t) {
    std::string date;
    
    std::size_t day = floor(t/(24*60*60));
    std::size_t hour = floor(t/(60*60)) - day*24;
    std::size_t min = floor(t/60) - day*24*60 - hour*60;
    std::size_t sec = floor(t) - day*24*60*60 - hour*60*60 - min*60;
    
    if (day != 0) date += to_string(day)+'d';
    if (hour != 0) date += to_string(hour,2)+'h';
    if (min != 0) date += to_string(min,2)+'m';
    date += to_string(sec,2)+'s';
    
    if (date[0] == '0' && date.size() != 2) date.erase(0,1);
    
    return date;
}

using arguments_t = std::map<std::string, std::vector<std::string>>;
arguments_t read_arguments(int argc, char** argv) {
    arguments_t map;
    for (int i = 1; i < argc; ++i) {
        std::string s = strtrim(argv[i]);
        if (s.empty()) continue;
        if (s[0] != '-') {
            map[""].push_back(s);
            continue;
        }
        
        s.erase(0,1);
        std::size_t eqpos = s.find_first_of('=');
        if (eqpos != s.npos) {
            map[strtrim(s.substr(0, eqpos))].push_back(strtrim(s.substr(eqpos+1)));
        } else {
            map[s].push_back("");
        }
    }
    
    map[""];
    
    return map;
}

template<class T>
bool set_value(T& arg, const std::string& name, const arguments_t& args) {
    auto iter = args.find(name);
    if (iter != args.end()) {
        return from_string(arg, iter->second.back());
    }
    
    return false;
}

bool set_value(bool& arg, const std::string& name, const arguments_t& args) {
    auto iter = args.find(name);
    if (iter != args.end()) {
        if (iter->second.empty() || iter->second.back().empty()) {
            arg = true;
            return true;
        } else {
            return from_string(arg, iter->second.back());
        }
    }
    
    return false;
}

bool argument_is_provided(const std::string& name, const arguments_t& args) {
    return args.find(name) != args.end();
}

bool argument_is_provided(const arguments_t& args) {
    return args.find("") != args.end();
}

const std::vector<std::string>& argument_list(const arguments_t& args) {
    return args.find("")->second;
}

struct list_t {
    std::valarray<double> x;
    std::valarray<double> y;
    std::valarray<double> cy;
};

namespace fits {
    using namespace CCfits;
}

bool readlist(list_t& data, const std::string& file, bool isfits) {
    if (isfits) {
        try {
            fits::FITS in(file);
            fits::ExtHDU& table = in.currentExtension();
            fits::Column& ra   = table.column("RA");
            fits::Column& dec  = table.column("DEC");
            
            ra.read(data.x, 1);
            dec.read(data.y, 1);
            
            return true;
        } catch (fits::FitsException& e) {
            std::cout << "error: qxmatch: FITS exception occured when loading source list:" <<
                e.message() << std::endl;
            return false;
        }
    } else {
        std::ifstream in(file);
        std::size_t nline = 0;
        std::string line;
        while (std::getline(in, line)) ++nline;
        
        data.x.resize(nline);
        data.y.resize(nline);
        in.seekg(0);
        
        std::size_t i = 0;
        while (!in.eof()) {
            in >> data.x[i] >> data.y[i];
            ++i;
        }
        
        return true;
    }
}

void make_radian(list_t& data) {
    const double conv = 3.14159265359/180.0;
    const std::size_t n = data.x.size();
    for (std::size_t i = 0; i < n; ++i) {
        data.x[i] *= conv;
        data.y[i] *= conv;
    }
}

void make_cos(list_t& data) {
    data.cy.resize(data.y.size());
    data.cy = cos(data.y);
}

struct res_t {
    res_t() : id(-1), d(std::numeric_limits<double>::infinity()) {}

    int    id;
    double d;
    
    bool operator < (const res_t& t) const {
        return d < t.d;
    }
};

void thread_job(const list_t& list1, const list_t& list2, std::vector<res_t>& res, const std::size_t nth, 
    const std::size_t i0, const std::size_t i1, std::atomic<std::size_t>& cnt) {
    const std::size_t nobj2 = list2.x.size();
    for (std::size_t i = i0; i < i1; ++i) {
        for (std::size_t j = 0; j < nobj2; ++j) {
            // For each pair of source, compute a distance indicator.
            // Note that this is not the 'true' distance in arseconds, but this is sufficient
            // to find the nearest neighbors (the true distance is obtained by computing 
            // 2*asin(sqrt(sd))), but all these functions are monotonous, hence not applying
            // them does not change the relative distances).
            double sra = sin(0.5*(list2.x[j] - list1.x[i]));
            double sde = sin(0.5*(list2.y[j] - list1.y[i]));
            double sd = sde*sde + sra*sra*list2.cy[j]*list1.cy[i];
            
            // We then compare this new distance to the largest one that is in the Nth nearest
            // neighbor list. If it is lower than that, we insert it in the list, removing the
            // old one, and sort the whole thing so that the largest distance is as the end of
            // the list.
            auto& wr = res[(i+1)*nth-1]; 
            if (sd < wr.d) {
                wr.id = j;
                wr.d = sd;
                if (nth != 1) {
                    std::sort(res.begin() + i*nth, res.begin() + (i+1)*nth);
                }
            }
        }
        
        ++cnt;
    }
}

template<typename T>
void print_progress(std::size_t i, std::size_t nobj, T start) {
    auto now = std::chrono::high_resolution_clock::now();
    auto total = std::chrono::duration_cast<std::chrono::milliseconds>(now - start).count();
    double remaining = total*double(nobj)/i - total;

    const std::size_t ndash = 100;
    std::cout << "\r[" << std::string(floor(ndash*i/double(nobj)),'-') << std::string(ndash - floor(ndash*i/double(nobj)),' ') << "]";
    std::cout << " " << to_string(i, floor(log10(double(nobj))) + 1);
    std::cout << " " << to_string(std::size_t(floor(100.0*i/double(nobj))), 3) << "%, " << date_str(remaining/1000.0) << " left" << std::flush;
}

void print_help() {
    std::cout << "qxmatch v1.2\n";
    std::cout << "  usage: qxmatch [-options=...] file1 file2 output\n\n";
    
    std::cout << "  'file1' and 'file2' must be ASCII files containing only right ascension and\n";
    std::cout << "  declination coordinates of objects (both in degrees), organized in two columns.\n";
    std::cout << "  'output' is the file in which the result will be stored also as an ASCII file,\n";
    std::cout << "  containing a line per object in 'file1'. Each line contains 'N' columns,\n";
    std::cout << "  displaying the ID of the Nth closest object in 'file2'. 'N' is configurable\n";
    std::cout << "  using the command line arguments (see below) and defaults to 1.\n\n";
    
    std::cout << "  If using very high number of objects, it will be preferable to work with\n";
    std::cout << "  binary files instead of ASCII files (use the 'fits' command line flag). Such\n";
    std::cout << "  FITS files are expected to point to an extension with two vector columns named\n";
    std::cout << "  'RA' and 'DEC', each containing a single row. They can be created in IDL using:\n";
    std::cout << "          mwrfits, {ra:[...], dec:[...]}, \"filename.fits\", /create\n";
    std::cout << "  The output file is written in a similar way: it contains one extension named\n";
    std::cout << "  'RESULT_BINARY', which contains 'N' columns named 'ID'+i (where 'i' stands for\n"; 
    std::cout << "  the 'i'th nearest neighbor). It can be read in IDL using:\n";
    std::cout << "          res = mrdfits(\"output.fits\", 1)\n";
    std::cout << "          best = res.id1[...]\n\n";
    
    std::cout << "  List of available command line options:\n";
    std::cout << "   -verbose: set this flag to print additional information in the standard output\n";
    std::cout << "   -quiet: set this flag to prevent writing informations to the standard output\n";
    std::cout << "   -fits: set this flag to read and write files as FITS binary tables (faster)\n";
    std::cout << "   -nth=[number]: set this value to the number of closest neighbors you want to\n";
    std::cout << "                  retrieve (default: 1).\n";
    std::cout << "   -j=[number]: set this value to the number of concurrent threads you want to\n";
    std::cout << "                run (default: 1).\n";
    std::cout << "   -radian: set this flag if coordinates are provided in radian (faster).\n\n";
    std::cout << "   -dist: set this flag to also output the crossmatch distance. If reading ASCII\n";
    std::cout << "          data, the distance will be saved as an aditionnal column along each ID\n";
    std::cout << "          (ID1, D1, ID2, D2, ...). If reading FITS data, 'N' new columns are added\n";
    std::cout << "          and called 'D'+i (where 'i' stands for the 'i'th neighbor). In both\n";
    std::cout << "          cases, if 'radian' is set, the distance will be saved in radians, else\n";
    std::cout << "          it will be saved in arcseconds.\n";
    
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
    if (argc < 2) {
        print_help();
        return 0;
    }

    // Mandatory parameters
    std::string sfile1 = "";
    std::string sfile2 = "";
    std::string soutput = "";
    
    // Optional parameters
    std::size_t nth = 1;
    std::size_t nthread = 1;
    bool        radians = false;
    bool        verbose = false;
    bool        quiet = false;
    bool        out_dist = false;
    bool        isfits = false;
    
    // Read configuration from the command line arguments
    auto args = read_arguments(argc, argv);
    
    // Get the file names
    auto arglist = argument_list(args);
    if (arglist.size() != 3) {
        print_help();
        return 1;
    } else {
        sfile1  = arglist[0];
        sfile2  = arglist[1];
        soutput = arglist[2];
    }
    
    // Fetch optional keywords and parameters
    set_value(nth,      "nth",     args);
    set_value(nthread,  "j",       args);
    set_value(radians,  "radians", args);
    set_value(verbose,  "verbose", args);
    set_value(quiet,    "quiet",   args);
    set_value(isfits,   "fits",    args);
    set_value(out_dist, "dist",    args);

    // Read the two lists
    if (verbose) {
        std::cout << "reading object lists..." << std::endl;
    }
    
    list_t list1, list2; 
    if (!readlist(list1, sfile1, isfits) || !readlist(list2, sfile2, isfits)) {
        return 1;
    }
    
    if (verbose) {
        std::cout << "cross-matching " << list1.x.size() << " objects ('" << sfile1 << "') with "
            << list2.x.size() << " objects ('" << sfile2 << "')" << std::endl;
    }
    
    // By default, coordinates are expected in degrees, but we need radians for the computation
    if (!radians) {
        if (verbose) {
            std::cout << "converting to radians..." << std::endl;
        }
    
        make_radian(list1);
        make_radian(list2);
    }
    
    // To speed up the computation, we can precompute a few values
    if (verbose) {
        std::cout << "computing cosines..." << std::endl;
    }
    
    make_cos(list1);
    make_cos(list2);
    
    std::vector<res_t> res(list1.x.size()*nth);
    
    // And now do the computation itself
    if (verbose) {
        std::cout << "find best matches..." << std::endl;
    }
    
    if (nthread == 1) {
        // When using a single thread, all the work is done in the main thread
        const std::size_t nobj1 = list1.x.size();
        const std::size_t nobj2 = list2.x.size();
        
        auto start = std::chrono::high_resolution_clock::now();
        auto last = start;
        
        for (std::size_t i = 0; i < nobj1; ++i) {
            for (std::size_t j = 0; j < nobj2; ++j) {
                // For each pair of source, compute a distance indicator.
                // Note that this is not the 'true' distance in arseconds, but this is sufficient
                // to find the nearest neighbors (the true distance is obtained by computing 
                // 2*asin(sqrt(sd))), but all these functions are monotonous, hence not applying
                // them does not change the relative distances).
                double sra = sin(0.5*(list2.x[j] - list1.x[i]));
                double sde = sin(0.5*(list2.y[j] - list1.y[i]));
                double sd = sde*sde + sra*sra*list2.cy[j]*list1.cy[i];
                
                // We then compare this new distance to the largest one that is in the Nth nearest
                // neighbor list. If it is lower than that, we insert it in the list, removing the
                // old one, and sort the whole thing so that the largest distance goes as the end of
                // the list.
                auto& wr = res[(i+1)*nth-1];
                if (sd < wr.d) {
                    wr.id = j;
                    wr.d = sd; 
                    if (nth != 1) {
                        std::sort(res.begin() + i*nth, res.begin() + (i+1)*nth);
                    }
                }
            }
            
            if (!quiet) {
                auto now = std::chrono::high_resolution_clock::now();
                auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(now - last).count();
                
                if (i == 0 || elapsed > 1000.0) {
                    print_progress(i, nobj1, start);
                    last = now;
                }
            }
        }
    } else {
        // When using more than one thread, the work load is evenly separated between all the
        // available threads, such that they should more or less all end at the same time.
        std::atomic<std::size_t> iter(0);
    
        // Create the thread pool and launch the threads
        std::vector<std::unique_ptr<std::thread>> pool(nthread);
        const std::size_t nobj = list1.x.size();
        std::size_t assigned_total = 0;
        for (std::size_t t = 0; t < nthread; ++t) {
            std::size_t assigned = floor(nobj/float(nthread));
            if (t == nthread-1) {
                assigned = nobj - assigned_total;
            }
            
            pool[t] = std::unique_ptr<std::thread>(new std::thread(&thread_job, list1, list2, std::ref(res), nth, assigned_total, assigned_total + assigned, std::ref(iter)));
            assigned_total += assigned;
        }
        
        // Wait for the computation to finish
        // Here the main thread does nothing except sleeping, occasionally waking up every second to
        // update the progress bar if any.
        auto start = std::chrono::high_resolution_clock::now();
        while (iter < nobj) {
            std::this_thread::sleep_for(std::chrono::seconds(1));
            if (!quiet) {
                print_progress(iter, nobj, start);
            }
        }
        
        // By now, all thread should have ended their tasks.
        // We must ask them to terminate nicely.
        for (std::size_t t = 0; t < nthread; ++t) {
            pool[t]->join();
        }
    }
    
    if (!quiet) {
        std::cout << std::endl;
    }
    
    // The computation is finished, now we only need to store the result somewhere
    if (verbose) {
        std::cout << "saving results..." << std::endl;
    }
    
    if (isfits) {
        try {
            fits::FITS out("!"+soutput, fits::Write);
            auto* tab = out.addTable("RESULT_BINARY", 1);
            
            const std::size_t nobj = list1.x.size();
            for (std::size_t n = 0; n < nth; ++n) {
                tab->addColumn(fits::Tint, "ID"+to_string(n+1), nobj);
                std::valarray<std::size_t> tres(nobj);
                for (std::size_t i = 0; i < nobj; ++i) {
                    tres[i] = res[i*nth+n].id;
                }
                
                tab->column("ID"+to_string(n+1)).write(tres, 1, 1);
                
                if (out_dist) {
                    tab->addColumn(fits::Tdouble, "D"+to_string(n+1), nobj);
                    std::valarray<double> tresd(nobj);
                    if (radians) {
                        for (std::size_t i = 0; i < nobj; ++i) {
                            tresd[i] = 2.0*asin(sqrt(res[i*nth+n].d));
                        }
                    } else {
                        for (std::size_t i = 0; i < nobj; ++i) {
                            tresd[i] = 3600.0*(180.0/3.14159265359)*2.0*asin(sqrt(res[i*nth+n].d));
                        }
                    }
                    
                    tab->column("D"+to_string(n+1)).write(tresd, 1, 1);
                }
            }
    
            if (verbose) {
                std::cout << "done." << std::endl;
            }
            
            return 0;
        } catch (fits::FitsException& e) {
            std::cout << "warning: qxmatch: FITS exception occured when saving result:\n" << e.message()
                << "note: qxmatch: using ASCII output as a fall back solution" << std::endl;
        }
    } 
    
    std::ofstream out(soutput);
    const std::size_t nobj = list1.x.size();
    for (std::size_t i = 0; i < nobj; ++i) {
        for (std::size_t n = 0; n < nth; ++n) {
            out << '\t' << res[i*nth + n].id;
            if (out_dist) {
                if (radians) {
                    out << '\t' << 2.0*asin(sqrt(res[i*nth + n].d));
                } else {
                    out << '\t' << 3600.0*(180.0/3.14159265359)*2.0*asin(sqrt(res[i*nth + n].d));
                }
            }
        }
        out << '\n';
    }
    
    if (verbose) {
        std::cout << "done." << std::endl;
    }
    
    return 0;
}

