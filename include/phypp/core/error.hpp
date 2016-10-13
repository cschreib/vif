#ifndef PHYPP_CORE_ERROR_HPP
#define PHYPP_CORE_ERROR_HPP

#include "phypp/core/print.hpp"
#include "phypp/core/typedefs.hpp"

#ifndef NO_LIBUNWIND
 #include <vector>
 #include <cxxabi.h>
 #define UNW_LOCAL_ONLY
 #include <libunwind.h>
 #ifndef NO_LIBDWARF
  #include <dwarf.h>
  #include <libdwarf.h>
  #include <sys/types.h>
  #include <fcntl.h>
  #include <unistd.h>
 #endif
#endif

// Store the path to the executable, to get debug information in case of error
extern const char* executable_path;

#ifndef NO_LIBUNWIND
#ifndef NO_LIBDWARF
// Inspired from:
// http://eli.thegreenplace.net/2011/02/07/how-debuggers-work-part-3-debugging-information
// ... and the code of add2line from:
// https://sourceforge.net/projects/elftoolchain/
struct dwarf_helper {
    struct source_info {
        Dwarf_Addr addr;
        Dwarf_Unsigned line;
        std::string filename;

        bool operator< (const source_info& s) const {
            return addr > s.addr; // reverse order is intentionnal
        }
        bool operator< (Dwarf_Addr a) const {
            return addr > a; // reverse order is intentionnal
        }
    };

    std::vector<source_info> db;

    bool get_info(Dwarf_Addr addr, source_info& info) const {
        auto iter = std::lower_bound(db.begin(), db.end(), addr);
        if (iter != db.end()) {
            info = *iter;
            return true;
        } else {
            return false;
        }
    }

    explicit dwarf_helper() {
        // Initialize dwarf
        int fd = -1;
        if ((fd = open(executable_path, O_RDONLY)) < 0) {
            return;
        }

        Dwarf_Debug dbg = 0;
        Dwarf_Error err;
        if (dwarf_init(fd, DW_DLC_READ, 0, 0, &dbg, &err) != DW_DLV_OK) {
            close(fd);
            return;
        }

        // Iterate through compilation unit (CU) headers
        while (true) {
            // Find next CU
            Dwarf_Unsigned cu_header_length, abbrev_offset, next_cu_header;
            Dwarf_Half version_stamp, offset_size, extension_size, address_size;
            int curc = dwarf_next_cu_header_b(dbg, &cu_header_length, &version_stamp, &abbrev_offset,
                &address_size, &offset_size, &extension_size, &next_cu_header, &err);
            if (curc == DW_DLV_ERROR || curc == DW_DLV_NO_ENTRY) {
                break;
            }

            // Get the CU DIE
            Dwarf_Die cu_die = 0, tmp_die;
            while (dwarf_siblingof(dbg, cu_die, &tmp_die, &err) == DW_DLV_OK) {
                cu_die = tmp_die;

                Dwarf_Half tag;
                if (dwarf_tag(cu_die, &tag, &err) != DW_DLV_OK) {
                    cu_die = 0;
                    break;
                }

                if (tag == DW_TAG_compile_unit) {
                    break;
                }
            }

            if (cu_die == 0) continue;

            // Grab the mappings between address and source lines in this CU
            Dwarf_Line* lines;
            Dwarf_Signed nlines = 0;
            if (dwarf_srclines(cu_die, &lines, &nlines, &err) == DW_DLV_ERROR) {
                continue;
            }

            for (Dwarf_Signed i = 0; i < nlines; ++i) {
                char* filename;
                if (dwarf_linesrc(lines[i], &filename, &err) != DW_DLV_OK) {
                    continue;
                }

                Dwarf_Unsigned lineno;
                if (dwarf_lineno(lines[i], &lineno, &err) != DW_DLV_OK) {
                    continue;
                }

                Dwarf_Addr lineaddr;
                if (dwarf_lineaddr(lines[i], &lineaddr, &err) != DW_DLV_OK) {
                    continue;
                }

                db.push_back({lineaddr, lineno, std::string(filename)});
            }

            dwarf_srclines_dealloc(dbg, lines, nlines);
        }

        // Close dwarf
        dwarf_finish(dbg, &err);
        close(fd);

        // Sort line database for lookup later on
        std::sort(db.begin(), db.end());
    }
};
#endif

namespace impl {
    inline std::string dbg_pretty_demangled(std::string name) {
        // Functions from "string.hpp"
        auto local_replace = [](std::string s, std::string pattern, std::string rep) {
            auto p = s.find(pattern);
            while (p != s.npos) {
                s.replace(p, pattern.size(), rep);
                p = s.find(pattern, p+rep.size());
            }

            return s;
        };

        auto local_split = [](std::string s, std::string pattern) {
            std::vector<std::string> ret;
            std::size_t p = 0, op = 0;
            while ((p = s.find(pattern, op)) != s.npos) {
                ret.push_back(s.substr(op, p - op));
                op = p + pattern.size();
            }

            ret.push_back(s.substr(op));
            return ret;
        };

        auto local_trim = [](std::string s, std::string chars = " \t") {
            std::size_t spos = s.find_first_of(chars);
            if (spos == 0) {
                std::size_t epos = s.find_first_not_of(chars);
                if (epos == s.npos) return std::string{};
                s = s.substr(epos);
            }

            spos = s.find_last_of(chars);
            if (spos == s.size()-1) {
                std::size_t epos = s.find_last_not_of(chars);
                s = s.erase(epos+1, s.size() - epos+1);
            }

            return s;
        };

        // Use type shortcuts
        name = local_replace(name, "unsigned long", "uint_t");

        // Remove spaces between template arguments and remove "enable_if" void parameter
        std::string oname = name;
        name.clear();
        auto p = oname.find('<');
        auto p0 = p*0;
        while (p != oname.npos) {
            name += oname.substr(p0, p-p0+1);
            uint_t open = 1;
            ++p;

            auto p1 = p;
            for (; p < oname.size() && open != 0; ++p) {
                if (oname[p] == '<') ++open;
                if (oname[p] == '>') --open;
            }

            auto params = local_split(oname.substr(p1, p-p1-1), ",");
            for (uint_t i = 0; i != params.size(); ++i) {
                params[i] = local_trim(params[i]);
                if (i == params.size()-1 && params[i] == "void") continue;
                if (i != 0) name += ",";
                name += params[i];
            }

            name += '>';

            p0 = p;
            p = oname.find('<', p0);
        }

        name += oname.substr(p0);

        // Simplify integral constants ("1ul" -> "1")
        oname = name;
        name.clear();
        p = oname.find("ul");
        p0 = p*0;
        while (p != oname.npos) {
            if (p == 0 || oname[p-1] < '0' || oname[p-1] > '9') {
                p += 2;
                p = oname.find("ul", p);
                continue;
            }

            name += oname.substr(p0, p-p0);
            p0 = p+2;
            p = oname.find("ul", p0);
        }

        name += oname.substr(p0);

        return name;
    }

    inline std::string dbg_pretty_filename(std::string filename) {
        auto p = filename.find_last_of('/');
        if (p != filename.npos) {
            std::string reflex_header = "._reflex_";
            if (filename.substr(p+1, reflex_header.size()) == reflex_header) {
                filename = filename.substr(p+1+reflex_header.size());
            }
        }

        return filename;
    }
}

// Inspired from:
// http://eli.thegreenplace.net/2015/programmatic-access-to-the-call-stack-in-c/
inline void backtrace() {
    unw_cursor_t cursor;
    unw_context_t context;

    // Initialize cursor to current frame for local unwinding.
    unw_getcontext(&context);
    unw_init_local(&cursor, &context);

    #ifndef NO_LIBDWARF
    // Initialize DWARF debug information
    dwarf_helper debug_info;
    #endif

    // Unwind frames one by one, going up the frame stack.
    while (unw_step(&cursor) > 0) {
        unw_word_t offset, pc;
        unw_get_reg(&cursor, UNW_REG_IP, &pc);
        if (pc == 0) {
            break;
        }

        char sym[2048];
        if (unw_get_proc_name(&cursor, sym, sizeof(sym), &offset) == 0) {
            std::string fname = sym;

            // Try to get demangled name from C++ ABI
            int status;
            char* demangled = abi::__cxa_demangle(sym, nullptr, nullptr, &status);
            if (status == 0) {
                fname = impl::dbg_pretty_demangled(demangled);
                free(demangled);
            }

            // Filter out system calls
            if (fname == "__libc_start_main" || fname == "_start") continue;

            #ifndef NO_LIBDWARF
            // auto iter = debug_info.db.find(sym);
            dwarf_helper::source_info info;
            if (debug_info.get_info(pc, info)) {
                // Print function name and pointer offset
                print(" - ", fname);
                // Print file and line number
                std::string filename = impl::dbg_pretty_filename(info.filename);
                print("   at ", filename, ":", info.line);
            } else {
                // Print function name and pointer offset
                print(" - ", fname, " (+", offset, ")");
            }
            #else
            // Print function name
            print(" - ", fname, " (+", offset, ")");
            #endif
        } else {
            print(" - unknown function (missing debug information?)");
        }
    }
}

#define phypp_check(value, ...) \
    if (!(value)) { \
        print("\n"); \
        error(__VA_ARGS__, "\n"); \
        print("backtrace:"); \
        backtrace(); \
        print("\n"); \
        exit(EXIT_FAILURE); \
    }
#else
#define phypp_check(value, ...) \
    if (!(value)) { \
        print("\n"); \
        error(__FILE__, ":", __LINE__, ": ", __VA_ARGS__, "\n"); \
        exit(EXIT_FAILURE); \
    }
#endif

#endif
