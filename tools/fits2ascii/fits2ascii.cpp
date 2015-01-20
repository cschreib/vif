#include <phypp.hpp>

template<std::size_t D, typename T>
void read_column_(const std::string& filename, const std::string& colname,
    vec_t<D,std::string>& v, type_list<T>) {

    vec_t<D,T> tmp;
    fits::read_table(filename, colname, tmp);
    v = strna(tmp);
}

template<std::size_t D>
void read_column_(const std::string& filename, const std::string& colname,
    vec_t<D,std::string>& v, type_list<std::string>) {

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
    vec_t<D,std::string>& v, uint_t& nrow, const vec1u& ids, type_list<T>) {

    vec_t<D,T> tmp;
    fits::read_table(filename, colname, tmp);

    auto next = check_rows(colname, tmp.dims[0], nrow, ids);
    if (next != do_next::run) {
        return next;
    }

    if (!ids.empty()) {
        tmp = tmp(ids, repeat<D-1>(_));
    }

    v = strna(tmp);

    return do_next::run;
}

template<std::size_t D>
do_next read_column_(const std::string& filename, const std::string& colname,
    vec_t<D,std::string>& v, uint_t& nrow, const vec1u& ids, type_list<std::string>) {

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
    vec_t<D,std::string>& v) {

    read_column_<D>(filename, colname, v, type_list<T>{});
}

template<std::size_t D, typename T>
do_next read_column(const std::string& filename, const std::string& colname,
    vec_t<D,std::string>& v, uint_t& nrow, const vec1u& ids) {

    return read_column_<D>(filename, colname, v, nrow, ids, type_list<T>{});
}

template<std::size_t Dim>
do_next read_column(const std::string& in_file, const fits::column_info& cinfo,
    vec_t<Dim,std::string>& v, std::string& vtype, uint_t& nrow, const vec1u& ids) {

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

void print_help();

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    vec1s tsplit;
    vec1s exclude, include;
    std::string ids_file;
    bool verbose = false;

    read_args(argc-2, argv+2, arg_list(name(tsplit, "split"), exclude, include,
        name(ids_file, "ids"), verbose));

    std::string in_file = argv[1];
    std::string out_file = argv[2];
    auto cols = fits::read_table_columns(in_file);

    uint_t nsplit = tsplit.size();
    vec1s split_cols(nsplit);
    std::vector<vec1s> split_names(nsplit);
    for (uint_t i : range(nsplit)) {
        vec1s ss = split(tsplit[i], ":");
        split_cols[i] = toupper(ss[0]);

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
            if (start_with(ss[j], "[")) {
                vec1s st = split(erase_end(erase_begin(ss[j], "["), "]"), ",");
                if (st.size() != ncols2) {
                    error("incompatible array for '", split_cols[i], "' split name (got ",
                        st.size(), " elements for ", ncols2, " columns)");
                    return 1;
                }

                split_names[i] += st;
            } else if (start_with(ss[j], "#")) {
                std::string col = toupper(erase_begin(ss[j], "#"));
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
                split_names[i] += tolower(st);
            } else {
                split_names[i] += ss[j];
            }
        }

        vec1u ide = where(empty(split_names[i]));
        split_names[i][ide] = split_cols[i] + '_' + strn(ide);
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
            if (!exclude.empty() && match_any_of(tolower(c.name), exclude)) continue;
        } else {
            if (!match_any_of(tolower(c.name), include)) continue;
        }

        if (verbose) {
            print("found column: ", c.name);
        }

        if (c.dims.size() == 1) {
            std::string type;
            vec1s v;

            auto next = read_column<1>(in_file, c, v, type, nrow_base, id_sel);
            if (next == do_next::skip) continue;
            if (next == do_next::abort) return 1;

            out_cols.push_back(v);
            names.push_back(tolower(c.name));
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
    print("fits2ascii v1.0");
    print("[WIP]");
}
