#include <phypp.hpp>

template<std::size_t D, typename T>
vec_t<D,std::string> read_column_(const std::string& filename, const std::string& colname, type_list<T>) {
    vec_t<D,T> v;
    fits::read_table(filename, colname, v);
    return strna(v);
}

template<std::size_t D>
vec_t<D,std::string> read_column_(const std::string& filename, const std::string& colname, type_list<std::string>) {
    vec_t<D,std::string> v;
    fits::read_table(filename, colname, v);
    return v;
}

template<std::size_t D, typename T>
vec_t<D,std::string> read_column(const std::string& filename, const std::string& colname) {
    return read_column_<D>(filename, colname, type_list<T>{});
}

void print_help();

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    vec1s tsplit;
    vec1s exclude, include;

    read_args(argc-2, argv+2, arg_list(name(tsplit, "split"), exclude, include));

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

                vec1s st = read_column<1,std::string>(in_file, col);
                split_names[i] += tolower(st);
            } else {
                split_names[i] += ss[j];
            }
        }

        vec1u ide = where(empty(split_names[i]));
        split_names[i][ide] = split_cols[i] + '_' + strn(ide);
    }

    vec1s names;
    vec1s types;
    std::vector<vec1s> out_cols;
    for (auto& c : cols) {
        if (include.empty()) {
            if (!exclude.empty() && match_any_of(tolower(c.name), exclude)) continue;
        } else {
            if (!match_any_of(tolower(c.name), include)) continue;
        }

        if (c.dims.size() == 1) {
            names.push_back(tolower(c.name));
            std::string type;
            switch (c.type) {
            case fits::column_info::string :
                out_cols.push_back(read_column<1,std::string>(in_file, c.name));
                type = "string";
                break;
            case fits::column_info::boolean :
                out_cols.push_back(read_column<1,bool>(in_file, c.name));
                type = "bool";
                break;
            case fits::column_info::byte :
                out_cols.push_back(read_column<1,char>(in_file, c.name));
                type = "byte";
                break;
            case fits::column_info::integer :
                out_cols.push_back(read_column<1,int_t>(in_file, c.name));
                type = "int";
                break;
            case fits::column_info::float_simple :
                out_cols.push_back(read_column<1,float>(in_file, c.name));
                type = "float";
                break;
            case fits::column_info::float_double :
                out_cols.push_back(read_column<1,double>(in_file, c.name));
                type = "double";
                break;
            }

            types.push_back("["+type+"]");
        } else if (c.dims.size() == 2) {
            vec1u i2c = where(split_cols == c.name);
            if (i2c.empty()) {
                warning("2-dimensional column '", c.name, "' not treated");
                note("use the 'split' program argument to split it into single columns");
                continue;
            }

            vec2s cols2;
            std::string type;
            switch (c.type) {
            case fits::column_info::string :
                cols2 = read_column<2,std::string>(in_file, c.name);
                type = "string";
                break;
            case fits::column_info::boolean :
                cols2 = read_column<2,bool>(in_file, c.name);
                type = "bool";
                break;
            case fits::column_info::byte :
                cols2 = read_column<2,char>(in_file, c.name);
                type = "byte";
                break;
            case fits::column_info::integer :
                cols2 = read_column<2,int_t>(in_file, c.name);
                type = "int";
                break;
            case fits::column_info::float_simple :
                cols2 = read_column<2,float>(in_file, c.name);
                type = "float";
                break;
            case fits::column_info::float_double :
                cols2 = read_column<2,double>(in_file, c.name);
                type = "double";
                break;
            }

            if (split_names[i2c[0]].size() != cols2.dims[1]) {
                error("incompatible name list and column dimension (",
                    split_names[i2c[0]].size(), " vs ", cols2.dims[1], ")");
                return 1;
            }

            append(names, split_names[i2c[0]]);
            append(types, replicate("["+type+"]", cols2.dims[1]));
            for (uint_t i : range(cols2.dims[1])) {
                out_cols.push_back(cols2(_,i));
            }
        } else {
            warning(c.dims.size(), "-dimensional column '", c.name, "' is not supported");
            note("skipping");
        }
    }

    vec1b skip(out_cols.size());
    uint_t nrow = 0;
    vec1u colsize(out_cols.size());
    for (uint_t ic : range(out_cols.size())) {
        if (nrow == 0) {
            nrow = out_cols[ic].size();
        } else if (out_cols[ic].size() != nrow) {
            warning("incompatible number of rows for column '", names[ic], "' (got ",
                out_cols[ic].size(), ", expected ", nrow, ")");
            note("skipping");
            skip[ic] = true;
            continue;
        }

        colsize[ic] = max(vec1u{
            names[ic].length(), types[ic].length(), max(length(out_cols[ic]))
        });
    }

    const uint_t padding = 2;
    for (uint_t ic : range(out_cols.size())) {
        if (skip[ic]) continue;
        names[ic] = align_left(names[ic], colsize[ic]+padding);
        types[ic] = align_left(types[ic], colsize[ic]+padding);
        out_cols[ic] = align_left(out_cols[ic], colsize[ic]+padding);
    }

    vec1u idi = where(!skip);

    std::ofstream out(out_file);
    out << "# " << collapse(names[idi]) << "\n";
    out << "# " << collapse(types[idi]) << "\n";
    out << "# \n";

    for (uint_t j : range(nrow)) {
        out << "  ";
        for (uint_t ic : range(out_cols.size())) {
            if (skip[ic]) continue;
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
