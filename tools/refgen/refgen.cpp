#include <clang-c/Index.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>
#include <algorithm>

std::string location_str(CXSourceLocation s) {
    std::ostringstream ss;
    CXFile f;
    unsigned int line;
    unsigned int column;
    unsigned int offset;
    clang_getSpellingLocation(s, &f, &line, &column, &offset);
    ss << line << ":" << column << ":" << offset;

    return ss.str();
}

template<typename T>
void get_location(T& t, CXSourceRange s) {
    std::ostringstream ss;
    CXFile f;
    unsigned int offset;
    clang_getSpellingLocation(clang_getRangeStart(s), &f, &t.lstart, &t.cend, &offset);
    clang_getSpellingLocation(clang_getRangeEnd(s), &f, &t.lend, &t.cstart, &offset);

    // For some reason, onle line declarations have start after end...
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
std::string cpp;

CXChildVisitResult visitor(CXCursor cursor, CXCursor parent, CXClientData client_data) {
    CXCursorKind cKind = clang_getCursorKind(cursor);
    if (cKind == CXCursor_StructDecl || cKind == CXCursor_ClassDecl) {
        if (get_file_name(cursor) == cpp) {
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
        if (get_file_name(cursor) == cpp) {
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
        if (!db.empty() && get_file_name(cursor) == cpp) {
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

int main(int argc, char* argv[]) {
    bool verbose = false;

    // "clang file.cpp -emit-ast -o file.ast"
    std::string ast = std::string(argv[1]) + ".ast";
    cpp = std::string(argv[1]) + ".cpp";

    // Parse the file with clang
    CXIndex cidx = clang_createIndex(0, 0);
    CXTranslationUnit ctu = clang_createTranslationUnit(cidx, ast.c_str());
    clang_visitChildren(clang_getTranslationUnitCursor(ctu), &visitor, nullptr);
    clang_disposeTranslationUnit(ctu);
    clang_disposeIndex(cidx);

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
    std::ifstream code(cpp);
    std::string rname = "._reflex_"+cpp;

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
