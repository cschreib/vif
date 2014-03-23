#include <clang-c/Index.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>
#include <unistd.h>
#include <sys/ioctl.h>

std::ostream& out = std::cerr;

// Note: code for coloring the terminal is inspired from LLVM
#define COLOR(FGBG, CODE, BOLD) "\033[0;" BOLD FGBG CODE "m"

#define ALLCOLORS(FGBG,BOLD) {\
    COLOR(FGBG, "0", BOLD),\
    COLOR(FGBG, "1", BOLD),\
    COLOR(FGBG, "2", BOLD),\
    COLOR(FGBG, "3", BOLD),\
    COLOR(FGBG, "4", BOLD),\
    COLOR(FGBG, "5", BOLD),\
    COLOR(FGBG, "6", BOLD),\
    COLOR(FGBG, "7", BOLD)\
  }

static const char color_codes[2][8][10] = {
    ALLCOLORS("3",""), ALLCOLORS("3","1;")
};

#undef COLOR
#undef ALLCOLORS

namespace color {
    enum color_value : char {
        black = 0, red, green, yellow, blue, magenta, cyan, white, normal = -1
    };

    std::ostream& reset(std::ostream& o) {
        return o << "\033[0m";
    }

    std::ostream& bold(std::ostream& o) {
        return o << "\033[1m";
    }

    struct set {
        color_value col_;
        bool bold_;

        explicit set(color_value col, bool bold = false) : col_(col), bold_(bold) {}
    };

    std::ostream& operator << (std::ostream& o, set s) {
        if (s.col_ == color::normal) {
            reset(o);
            if (s.bold_) bold(o);
            return o;
        } else {
            return o << color_codes[s.bold_ ? 1 : 0][(char)s.col_ & 7];
        }
    }
}

template<typename T>
void get_location(T& t, CXSourceRange s) {
    CXFile f;
    unsigned int offset;
    clang_getSpellingLocation(clang_getRangeStart(s), &f, &t.lstart, &t.cend, &offset);
    clang_getSpellingLocation(clang_getRangeEnd(s), &f, &t.lend, &t.cstart, &offset);

    // For some reason, one line declarations have start after end...
    if (t.lstart >= t.lend && t.cstart > t.cend) {
        std::swap(t.lstart, t.lend);
        std::swap(t.cstart, t.cend);
    } else {
        // For some reason, the starting character is always one too much
        t.cstart -= 1;
        // and the ending one is always one too few
        t.cend += 1;
    }
}

std::string get_file_name(CXCursor c) {
    CXSourceRange r = clang_getCursorExtent(c);

    CXFile cf;
    unsigned int line;
    unsigned int column;
    unsigned int offset;
    clang_getSpellingLocation(clang_getRangeStart(r), &cf, &line, &column, &offset);
    CXString tsf = clang_getFileName(cf);
    const char* tt = clang_getCString(tsf);
    if (!tt) return "";

    std::string sf = tt;
    clang_disposeString(tsf);

    auto pos = sf.find_last_of('/');
    if (pos != sf.npos) {
        sf = sf.substr(pos+1);
    }

    return sf;
}

bool file_is_same(CXFile f1, CXFile f2) {
    CXFileUniqueID id1, id2;
    clang_getFileUniqueID(f1, &id1);
    clang_getFileUniqueID(f2, &id2);

    for (std::size_t i = 0; i < 3; ++i) {
        if (id1.data[i] != id2.data[i]) return false;
    }

    return true;
}

bool is_in_file(CXCursor c, CXFileUniqueID id) {
    CXSourceRange r = clang_getCursorExtent(c);

    CXFile cf;
    unsigned int line;
    unsigned int column;
    unsigned int offset;
    clang_getSpellingLocation(clang_getRangeStart(r), &cf, &line, &column, &offset);

    CXFileUniqueID idc;
    clang_getFileUniqueID(cf, &idc);

    for (std::size_t i = 0; i < 3; ++i) {
        if (id.data[i] != idc.data[i]) return false;
    }

    return true;
}

bool is_in_parent(CXCursor child, CXCursor parent) {
    CXSourceRange r1 = clang_getCursorExtent(parent);
    CXSourceRange r2 = clang_getCursorExtent(child);

    unsigned int s1, s2, e1, e2;

    CXFile f;
    unsigned int line;
    unsigned int column;
    clang_getSpellingLocation(clang_getRangeStart(r1), &f, &line, &column, &s1);
    clang_getSpellingLocation(clang_getRangeStart(r2), &f, &line, &column, &s2);
    clang_getSpellingLocation(clang_getRangeEnd(r1), &f, &line, &column, &e1);
    clang_getSpellingLocation(clang_getRangeEnd(r2), &f, &line, &column, &e2);

    return (s2 >= s1 && s2 < e1) && (e2 >= s1 && e2 < e1);
}

struct struct_t {
    std::string usr;
    std::string name;
    unsigned int lstart, lend;
    unsigned int cstart, cend;
    std::vector<std::string> members;

    bool operator < (const struct_t& s) const {
        return usr < s.usr;
    }
    bool operator < (const std::string& s) const {
        return usr < s;
    }
    bool operator == (const struct_t& s) const {
        return usr == s.usr;
    }
    bool operator == (const std::string& s) const {
        return usr == s;
    }
};

struct cpos_t {
    CXCursor cur;
    struct_t& str;
};

std::deque<struct_t> db;
std::vector<cpos_t> cstack;
CXFileUniqueID cpp;

CXChildVisitResult visitor(CXCursor cursor, CXCursor parent, CXClientData client_data) {
    CXCursorKind cKind = clang_getCursorKind(cursor);
    if (cKind == CXCursor_StructDecl || cKind == CXCursor_ClassDecl) {
        if (is_in_file(cursor, cpp)) {
            while (!cstack.empty() && !is_in_parent(cursor, cstack.back().cur)) {
                cstack.pop_back();
            }

            CXString nameString = clang_getCursorDisplayName(cursor);
            CXString usr = clang_getCursorUSR(cursor);
            if (std::find(db.begin(), db.end(), clang_getCString(usr)) == db.end()) {
                struct_t s;
                s.usr = clang_getCString(usr);
                s.name = clang_getCString(nameString);
                get_location(s, clang_getCursorExtent(cursor));
                db.push_back(s);

                cstack.push_back({cursor, db.back()});

                clang_disposeString(usr);
                clang_disposeString(nameString);
            } else {
                clang_disposeString(usr);
                clang_disposeString(nameString);

                return CXChildVisit_Continue;
            }
        } else {
            return CXChildVisit_Continue;
        }
    } else if (cKind == CXCursor_FieldDecl) {
        if (is_in_file(cursor, cpp)) {
            while (!cstack.empty() && !is_in_parent(cursor, cstack.back().cur)) {
                cstack.pop_back();
            }

            if (!cstack.empty()) {
                CXString nameString = clang_getCursorDisplayName(cursor);
                cstack.back().str.members.push_back(clang_getCString(nameString));
                clang_disposeString(nameString);
            }
        }
    } else if (cKind == CXCursor_VarDecl) {
        if (!db.empty() && is_in_file(cursor, cpp)) {
            struct_t s;
            get_location(s, clang_getCursorExtent(cursor));
            for (auto& st : db) {
                if (st.name.empty() && s.lstart == st.lend) {
                    CXString nameString = clang_getCursorDisplayName(cursor);
                    st.name = clang_getCString(nameString);
                    clang_disposeString(nameString);
                    break;
                }
            }
        }
    }

    return CXChildVisit_Recurse;
}

std::string location_str(CXSourceLocation s) {
    std::ostringstream ss;
    CXFile f;
    unsigned int line;
    unsigned int column;
    unsigned int offset;
    clang_getSpellingLocation(s, &f, &line, &column, &offset);
    CXString fname = clang_getFileName(f);
    ss << clang_getCString(fname) << ":" << line << ":" << column;
    clang_disposeString(fname);

    return ss.str();
}

std::size_t terminal_width() {
    struct winsize w;
    ioctl(STDERR_FILENO, TIOCGWINSZ, &w);
    return w.ws_col;
}

template<typename CharType>
bool is_any_of(CharType c, std::basic_string<CharType> chars) {
    return chars.find(c) != std::basic_string<CharType>::npos;
}

template<typename CharType>
std::basic_string<CharType> string_range(std::basic_string<CharType> str, std::size_t b, std::size_t e) {
    using string = std::basic_string<CharType>;
    if (e == string::npos) {
        return str.substr(b);
    } else {
        return str.substr(b, e-b);
    }
}

template<typename CharType>
std::basic_string<CharType> wrap(std::basic_string<CharType> msg,
    std::size_t head, std::size_t indent, std::size_t maxwidth) {
    using string = std::basic_string<CharType>;

    if (indent >= maxwidth) return "...";

    string res;
    bool first = true;

    if (head >= maxwidth) {
        // The header is too large already, no text will be written on the first line.
        first = false;
    }

    const string spaces = " \t";

    std::size_t line_begin = msg.find_first_not_of(spaces);
    std::size_t line_end = line_begin;

    while (line_begin + maxwidth - (first ? head : indent) < msg.size()) {
        std::size_t pos = line_begin + maxwidth - (first ? head : indent);
        auto c = msg[pos];
        std::size_t new_begin;
        if (is_any_of(c, spaces)) {
            // Clipping occurs in the middle of a empty region.
            // Look for the end of the previous word.
            line_end = msg.find_last_not_of(spaces, pos);
            if (line_end == string::npos) {
                // No word found, discard this line.
                line_begin = line_end;
            } else {
                // Word found, keep it on this line.
                ++line_end;
            }

            // Set the begining of the new line to the next word.
            new_begin = msg.find_first_not_of(spaces, pos);
        } else {
            // Clipping occurs in the middle of a word.
            // Look for the begining of this word.
            line_end = msg.find_last_of(spaces, pos);
            if (line_end == string::npos) {
                if (first) {
                    // The header is too large, just start the message on the next line
                    new_begin = line_begin;
                    line_begin = string::npos;
                } else {
                    // This is the only word for this line, we have no choice but to keep it,
                    // even if it is too long.
                    line_end = msg.find_first_of(spaces, pos);
                    new_begin = msg.find_first_not_of(spaces, line_end);
                }
            } else {
                // Keep this word as the beginning of the next line.
                ++line_end;
                new_begin = line_end;
            }
        }

        if (line_begin != string::npos) {
            if (!first) res += '\n'+string(indent, ' ');
            res += string_range(msg, line_begin, line_end);
        }

        line_begin = new_begin;
        first = false;

        if (line_begin == string::npos) {
            // No word remaining.
            break;
        }
    }

    if (line_begin != string::npos) {
        // There are some words remaining, put them on the last line.
        line_end = msg.find_last_not_of(spaces);
        if (line_end != string::npos) {
            ++line_end;

            if (!first) res += '\n'+string(indent, ' ');
            res += string_range(msg, line_begin, line_end);
        }
    }

    return res;
}

void format_diagnostic(CXDiagnostic d);

void format_diagnostics(CXDiagnosticSet ds) {
    for (std::size_t i = 0; i < clang_getNumDiagnosticsInSet(ds); ++i) {
        CXDiagnostic d = clang_getDiagnosticInSet(ds, i);
        format_diagnostic(d);
        clang_disposeDiagnostic(d);
    }
}

void format_range(CXDiagnostic d, CXSourceLocation sl) {
    CXFile floc;
    unsigned int lloc, cloc; {
        unsigned int offset;
        clang_getSpellingLocation(sl, &floc, &lloc, &cloc, &offset);
    }

    std::string filename; {
        CXString tmp = clang_getFileName(floc);
        filename = clang_getCString(tmp);
        clang_disposeString(tmp);
    }

    std::ifstream fs(filename);
    if (!fs.is_open()) return;

    std::string line;
    for (std::size_t i = 0; i < lloc; ++i) {
        std::getline(fs, line);
    }

    std::size_t coffset = 0;
    std::size_t width = line.size();
    std::size_t max_width = terminal_width();
    if (width > max_width) {
        // The line is too long to fit on the terminal.
        // Just truncate it for now (TODO: improve that?)
        if (max_width >= 3) {
            line.erase(max_width-3);
            line += "...";
        } else {
            line.erase(max_width);
        }

        width = max_width;
    }

    out << line << '\n';

    std::string highlight = std::string(width, ' ');

    std::size_t nrange = clang_getDiagnosticNumRanges(d);
    for (std::size_t i = 0; i < nrange; ++i) {
        CXSourceRange r = clang_getDiagnosticRange(d, i);

        CXFile fstart, fend;
        unsigned int lstart, lend, cstart, cend; {
            unsigned int offset;
            clang_getSpellingLocation(clang_getRangeStart(r), &fstart, &lstart, &cstart, &offset);
            clang_getSpellingLocation(clang_getRangeEnd(r), &fend, &lend, &cend, &offset);
        }

        if (cend < cstart) std::swap(cend, cstart);

        if (file_is_same(floc, fstart) && lloc == lstart &&
            cstart >= coffset && cend <= coffset+width) {
            for (std::size_t j = cstart; j < cend; ++j) {
                highlight[j-1 - coffset] = '~';
            }
        }
    }

    if (cloc-1 >= coffset && cloc-1 < coffset + width) {
        highlight[cloc-1 - coffset] = '^';
    }

    out << color::set(color::green, true) << highlight << color::reset << '\n';
}

void format_diagnostic(CXDiagnostic d) {
    auto col = color::normal;
    std::string kind = "";
    bool bold_message = true;

    switch (clang_getDiagnosticSeverity(d)) {
    case CXDiagnostic_Note :
        col = color::black;
        kind = "note: ";
        bold_message = false;
        break;
    case CXDiagnostic_Warning :
        col = color::magenta;
        kind = "warning: ";
        break;
    case CXDiagnostic_Error :
    case CXDiagnostic_Fatal :
        col = color::red;
        kind = "error: ";
        break;
    case CXDiagnostic_Ignored :
    default : return;
    }

    CXSourceLocation sl = clang_getDiagnosticLocation(d);

    std::string loc;
    if (clang_equalLocations(sl, clang_getNullLocation()) == 0) {
        loc = location_str(sl)+": ";
    }

    std::string message; {
        CXString tmp = clang_getDiagnosticSpelling(d);
        message = wrap(std::string(clang_getCString(tmp)),
            loc.size()+kind.size(), 6, terminal_width()
        );
        clang_disposeString(tmp);
    }

    out << color::set(color::normal, true) << loc
        << color::set(col, true) << kind
        << color::set(color::normal, bold_message) << message
        << color::reset << "\n";

    if (!loc.empty()) {
        format_range(d, sl);
    }

    format_diagnostics(clang_getChildDiagnostics(d));
    out << std::flush;
}

bool file_exist(const std::string& file) {
    std::fstream f(file);
    return f.is_open();
}

int main(int argc, char* argv[]) {
    bool verbose = false;

    std::string file = std::string(argv[1]);

    if (!file_exist(file)) {
        std::cout << color::set(color::red, true) << "error: " << color::reset <<
            "cannot find '" << file << "'" << std::endl;
        return 1;
    }

    std::string rname = argv[2];

    // Parse the file with clang
    CXIndex cidx = clang_createIndex(0, 0);
    CXTranslationUnit ctu = clang_parseTranslationUnit( // TODO: v optimize
        cidx, file.c_str(), argv+3, argc-3, 0, 0, CXTranslationUnit_None
    );

    CXFile main_file = clang_getFile(ctu, file.c_str());
    clang_getFileUniqueID(main_file, &cpp);

    std::size_t ndiag = clang_getNumDiagnostics(ctu);
    for (std::size_t i = 0; i < ndiag; ++i) {
        CXDiagnostic d = clang_getDiagnostic(ctu, i);
        format_diagnostic(d);
        clang_disposeDiagnostic(d);
    }

    if (ndiag == 0) {
        clang_visitChildren(clang_getTranslationUnitCursor(ctu), &visitor, nullptr);
    }

    clang_disposeTranslationUnit(ctu);
    clang_disposeIndex(cidx);

    if (ndiag != 0) return 1;

    // First sort the struct list by order of appearance in the file
    std::sort(db.begin(), db.end(), [](const struct_t& s1, const struct_t& s2) {
        return s1.lend < s2.lend;
    });

    // Output the list if required
    if (verbose) {
        for (auto iter = db.rbegin(); iter != db.rend(); ++iter) {
            std::cout << iter->name << " (" << iter->usr << ") "
                << iter->lstart << ":" << iter->cstart << " > "
                << iter->lend << ":" << iter->cend << std::endl;
            for (auto& m : iter->members) {
                std::cout << m << std::endl;
            }
        }
    }

    // Create a new code with added reflexion data
    std::ifstream code(file);

    std::ofstream enh(rname.c_str());
    std::size_t l = 1;
    while (!code.eof()) {
        if (db.empty()) {
            std::string line;
            std::getline(code, line);
            enh << line << '\n';
            ++l;
        } else {
            while (l < db.front().lend) {
                std::string line;
                std::getline(code, line);
                enh << line << '\n';
                ++l;
            }

            std::string line;
            std::getline(code, line);

            std::size_t rpos = 0;
            while (!db.empty() && l == db.front().lend) {
                enh << line.substr(rpos, db.front().cend-2 - rpos);

                if (db.front().members.empty()) {
                    enh << "NO_MEMBER(\"" << db.front().name << "\"); ";
                } else {
                    enh << "MEMBERS1(";
                    bool first = true;
                    for (auto& s : db.front().members) {
                        if (!first) enh << ", ";
                        enh << s;
                        first = false;
                    }
                    enh << "); MEMBERS2(\"";
                    enh << db.front().name << "\"";
                    for (auto& s : db.front().members) {
                        enh << ", MAKE_MEMBER(" << s << ")";
                    }
                    enh << "); ";
                }

                enh << line.substr(db.front().cend-2, 2);
                rpos = db.front().cend;
                db.pop_front();
            }

            enh << line.substr(rpos) << "\n";
            ++l;
        }
    }

    enh.flush();

    return 0;
}
