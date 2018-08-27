#include <vif.hpp>

using namespace vif;

template<std::size_t D, typename T>
void read_column_(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v, meta::type_list<T>) {

    vec<D,T> tmp;
    fits::read_table(filename, colname, tmp);

    if (std::is_same<T,float>::value) {
        v = to_string_vector(format::scientific(tmp));
    } else if (std::is_same<T,double>::value) {
        v = to_string_vector(format::precision(format::scientific(tmp), 17));
    } else {
        v = to_string_vector(tmp);
    }
}

template<std::size_t D>
void read_column_(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v, meta::type_list<std::string>) {

    fits::read_table(filename, colname, v);
}

enum class do_next {
    run,
    skip,
    abort
};

do_next check_rows(const std::string& colname, uint_t tmp_row, uint_t& nrow,
    const vec1u& ids) {

    bool first_try = false;
    if (nrow == 0) {
        nrow = tmp_row;
        first_try = true;
    } else if (tmp_row != nrow) {
        warning("incompatible number of rows for column '", colname, "' (got ",
            tmp_row, ", expected ", nrow, ")");
        note("this column will not be serialized");
        return do_next::skip;
    }

    if (first_try && !ids.empty()) {
        uint_t mid = max(ids);
        if (mid >= nrow) {
            error("maximum index in IDS is too large (got ", mid, ", against ",
                nrow, ")");
            return do_next::abort;
        }
    }

    return do_next::run;
}

template<std::size_t D, typename T>
do_next read_column_(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v, uint_t& nrow, const vec1u& ids, meta::type_list<T>) {

    vec<D,T> tmp;
    fits::read_table(filename, colname, tmp);

    auto next = check_rows(colname, tmp.dims[0], nrow, ids);
    if (next != do_next::run) {
        return next;
    }

    if (!ids.empty()) {
        tmp = tmp(ids, repeat<D-1>(_));
    }

    v = to_string_vector(tmp);

    return do_next::run;
}

template<std::size_t D>
do_next read_column_(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v, uint_t& nrow, const vec1u& ids, meta::type_list<std::string>) {

    fits::read_table(filename, colname, v);

    auto next = check_rows(colname, v.dims[0], nrow, ids);
    if (next != do_next::run) {
        return next;
    }

    if (!ids.empty()) {
        v = v(ids, repeat<D-1>(_));
    }

    return do_next::run;
}

template<std::size_t D, typename T>
void read_column(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v) {

    read_column_<D>(filename, colname, v, meta::type_list<T>{});
}

template<std::size_t D, typename T>
do_next read_column(const std::string& filename, const std::string& colname,
    vec<D,std::string>& v, uint_t& nrow, const vec1u& ids) {

    return read_column_<D>(filename, colname, v, nrow, ids, meta::type_list<T>{});
}

template<std::size_t Dim>
do_next read_column(const std::string& in_file, const fits::column_info& cinfo,
    vec<Dim,std::string>& v, std::string& vtype, uint_t& nrow, const vec1u& ids) {

    switch (cinfo.type) {
    case fits::column_info::string : {
        auto next = read_column<Dim,std::string>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "string";
        break;
    }
    case fits::column_info::boolean : {
        auto next = read_column<Dim,bool>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "bool";
        break;
    }
    case fits::column_info::byte : {
        auto next = read_column<Dim,char>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "byte";
        break;
    }
    case fits::column_info::integer : {
        auto next = read_column<Dim,int_t>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "int";
        break;
    }
    case fits::column_info::float_simple : {
        auto next = read_column<Dim,float>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "float";
        break;
    }
    case fits::column_info::float_double : {
        auto next = read_column<Dim,double>(in_file, cinfo.name, v, nrow, ids);
        if (next != do_next::run) return next;
        vtype = "double";
        break;
    }
    }

    return do_next::run;
}

bool regex_match_any_of(const std::string& str, const vec1s& rs) {
    for (auto& r : rs) {
        if (regex_match(str, r)) return true;
    }

    return false;
}

void print_help();

int vif_main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    vec1s tsplit;
    vec1s exclude, include, rename;
    std::string ids_file;
    bool verbose = false;
    uint_t hdu = 1;

    read_args(argc-2, argv+2, arg_list(name(tsplit, "split"), exclude, include,
        rename, name(ids_file, "ids"), verbose, hdu));

    std::string in_file = argv[1];
    std::string out_file = argv[2];
    auto cols = fits::read_table_columns(in_file, hdu);

    vec1s rename_from, rename_to;
    for (auto& n : rename) {
        vec1s spl = split(n, ":");
        if (spl.size() != 2) {
            error("wrong format for renaming '", n, "'");
            note("expected '<old name>:<new name>'");
            return 1;
        }

        rename_from.push_back(trim(spl[0]));
        rename_to.push_back(trim(spl[1]));
    }

    uint_t nsplit = tsplit.size();
    vec1s split_cols(nsplit);
    std::vector<vec1s> split_names(nsplit);
    for (uint_t i : range(nsplit)) {
        vec1s ss = split(tsplit[i], ":");
        split_cols[i] = to_upper(ss[0]);

        auto iterc = std::find_if(cols.begin(), cols.end(),
            [&](const fits::column_info& ci) {
            return ci.name == split_cols[i];
        });

        if (iterc == cols.end()) {
            warning("no column named '", split_cols[i], "'");
            continue;
        } else if (iterc->dims.size() != 2) {
            warning("column '", split_cols[i], "' is not 2-dimensional");
            continue;
        }

        uint_t ncols2 = iterc->dims[1];
        split_names[i].resize(ncols2);

        if (ss.size() == 1) continue;

        for (uint_t j : range(1, ss.size())) {
            if (begins_with(ss[j], "[")) {
                vec1s st = split(erase_end(erase_begin(ss[j], "["), "]"), ",");
                if (st.size() != ncols2) {
                    error("incompatible array for '", split_cols[i], "' split name (got ",
                        st.size(), " elements for ", ncols2, " columns)");
                    return 1;
                }

                split_names[i] += st;
            } else if (begins_with(ss[j], "#")) {
                std::string col = to_upper(erase_begin(ss[j], "#"));
                auto iter = std::find_if(cols.begin(), cols.end(),
                    [&](const fits::column_info& ci) {
                    return ci.name == col;
                });

                if (iter == cols.end()) {
                    error("unknown variable '", col, "'");
                    return 1;
                } else if (iter->dims.size() != 1) {
                    error("naming column '", col, "' must be 1-dimensional");
                    return 1;
                } else if (iter->dims[0] != ncols2) {
                    error("incompatible naming column for '", split_cols[i],
                        "' split name (got ", iter->dims[0], " elements for ", ncols2,
                        " columns)");
                    return 1;
                }

                vec1s st;
                read_column<1,std::string>(in_file, col, st);
                split_names[i] += to_lower(st);
            } else {
                split_names[i] += ss[j];
            }
        }

        vec1u ide = where(empty(split_names[i]));
        split_names[i][ide] = split_cols[i] + '_' + to_string(ide);
    }

    vec1u id_sel;
    if (!ids_file.empty()) {
        fits::read_table(ids_file, "ids", id_sel);
    }

    vec1s names;
    vec1s types;
    std::vector<vec1s> out_cols;
    uint_t nrow_base = 0;

    for (auto& c : cols) {
        if (include.empty()) {
            if (!exclude.empty() && regex_match_any_of(to_lower(c.name), exclude)) continue;
        } else {
            if (!regex_match_any_of(to_lower(c.name), include)) continue;
        }

        if (verbose) {
            print("found column: ", c.name, ", dimension: ", c.dims);
        }

        if (c.dims.size() == 1) {
            std::string type;
            vec1s v;

            auto next = read_column<1>(in_file, c, v, type, nrow_base, id_sel);
            if (next == do_next::skip) continue;
            if (next == do_next::abort) return 1;

            std::string name;
            vec1u idr = where(rename_from == c.name);
            if (!idr.empty()) {
                name = rename_to[idr.back()];
            } else {
                name = c.name;
            }

            out_cols.push_back(v);
            names.push_back(to_lower(name));
            types.push_back("["+type+"]");
        } else if (c.dims.size() == 2) {
            vec1u i2c = where(split_cols == c.name);
            if (i2c.empty()) {
                warning("2-dimensional column '", c.name, "' not treated");
                note("use the 'split' program argument to split it into single columns");
                continue;
            }

            std::string type;
            vec2s v;

            auto next = read_column<2>(in_file, c, v, type, nrow_base, id_sel);
            if (next == do_next::skip) continue;
            if (next == do_next::abort) return 1;

            if (split_names[i2c[0]].size() != v.dims[1]) {
                error("incompatible name list and column dimension (",
                    split_names[i2c[0]].size(), " vs ", v.dims[1], ")");
                return 1;
            }

            append(names, split_names[i2c[0]]);
            append(types, replicate("["+type+"]", v.dims[1]));
            for (uint_t i : range(v.dims[1])) {
                out_cols.push_back(v(_,i));
            }
        } else {
            warning(c.dims.size(), "-dimensional column '", c.name, "' is not supported");
            note("skipping");
        }
    }

    uint_t nrow = nrow_base;
    if (!id_sel.empty()) nrow = id_sel.size();

    if (verbose) note("serializing ", out_cols.size(), " columns of ", nrow, " rows");

    vec1u colsize(out_cols.size());
    for (uint_t ic : range(out_cols.size())) {
        colsize[ic] = max(vec1u{
            names[ic].length(), types[ic].length(), max(length(out_cols[ic]))
        });
    }

    const uint_t padding = 2;
    for (uint_t ic : range(out_cols.size())) {
        names[ic] = align_left(names[ic], colsize[ic]+padding);
        types[ic] = align_left(types[ic], colsize[ic]+padding);
        out_cols[ic] = align_left(out_cols[ic], colsize[ic]+padding);
    }

    std::ofstream out(out_file);
    out << "# " << collapse(names) << "\n";
    out << "# " << collapse(types) << "\n";
    out << "# \n";

    for (uint_t j : range(nrow)) {
        out << "  ";
        for (uint_t ic : range(out_cols.size())) {
            out << out_cols[ic][j];
        }
        out << "\n";
    }

    return 0;
}

void print_help() {
    using namespace terminal_format;

    print("fits2ascii v1.0");
    header("Usage: fits2ascii file.fits out.dat [options]");
    header("Available options:");
    bullet("include=[...]", "By default the program will export all columns to the ASCII file. "
        "You can specify in this variable an array of regular expressions that the program will "
        "use to identify which columns should be exported. The columns not matching the provided "
        "regular expressions will be ignored. Be careful that these are regular expressions, not "
        "simple search strings. In particular, the search pattern 'ra' will include columns named "
        "'ra', but also 'radius' or 'grade'. To only match the 'ra' column, use '^ra$'.");
    bullet("exclude=[...]", "Similar to 'include', but defines which columns to exclude rather "
        "than which columns to include. If both 'include' and 'exclude' are used, 'exclude' is "
        "ignored.");
    bullet("rename=[...]", "Each element of this array should be of the form 'xxx:yyy', which "
        "will change the name of column 'xxx' into 'yyy' in the ASCII file.");
    bullet("split=[...]", "Each element of this array should be of the form 'xxx:y:z:...'. This "
        "will split a 2D column named 'xxx' into multiple 1D columns. The other arguments ('y', "
        "'z', ...) determine the name of these multiple columns by contatenation ('y+z+...'). "
        "These can be: a) a simple string 'foo', b) an array '[foo,bar]' which must contain as "
        "many elements as 1D columns, and c) the name of a variable in the FITS file '#names' "
        "that contains as many strings as 1D columns. For example, if the FITS file contains a "
        "variable 'bands' that holds the name of all the bands in a 2D column 'fluxes', then this "
        "2D column can be exported using 'split=fluxes:flux_:#bands'. If 'bands' contains "
        "['f435w','f606w',...], the resulting 1D columns will be named 'flux_f435w', 'flux_f606w', "
        "etc. If the FITS file does not contain any variable you can use for naming columns, you "
        "can use the array synthax to specify these names manually. Using the same example as "
        "before, one would do instead 'split=fluxes:flux_:[f435w,f606w,...]'.");
    bullet("ids=...", "Must be a path to a valid FITS file containing a column'IDS', which "
        "will be used to define which rows of the catalog to save in the ASCII file. For example, "
        "if IDS=[0,1,2,3] then this program will only output the data for rows 0 to 3.");
    bullet("verbose", "Print the names of the exported columns in the terminal.");
}
