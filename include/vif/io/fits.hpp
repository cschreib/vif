#ifndef PHYPP_IO_FITS_HPP
#define PHYPP_IO_FITS_HPP

#include "phypp/io/fits/base.hpp"
#include "phypp/io/fits/table.hpp"
#include "phypp/io/fits/image.hpp"
#include "phypp/utility/os.hpp"

#ifndef NO_CFITSIO

namespace phypp {
namespace fits {
    // Return the number of dimensions of a FITS file
    // Note: will return 0 for FITS tables
    inline uint_t file_axes(const std::string& filename) {
        return fits::input_image(filename).axis_count();
    }

    // Return the dimensions of a FITS image
    // Note: will return an empty array for FITS tables
    inline vec1u file_dimensions(const std::string& filename) {
        return fits::input_image(filename).image_dims();
    }

    inline bool is_cube(const std::string& filename) {
        return fits::input_image(filename).is_cube();
    }

    inline bool is_image(const std::string& filename) {
        return fits::input_image(filename).is_image();
    }

    // Load the content of a FITS file into an array.
    template<std::size_t Dim, typename Type>
    void read_hdu(const std::string& filename, vec<Dim,Type>& v, uint_t hdu, fits::header& hdr) {
        fits::input_image img(filename, hdu);
        hdr = img.read_header();
        img.read(v);
    }

    template<std::size_t Dim, typename Type>
    void read_hdu(const std::string& filename, vec<Dim, Type>& v, uint_t hdu) {
        fits::input_image(filename, hdu).read(v);
    }

    template<std::size_t Dim, typename Type>
    void read(const std::string& filename, vec<Dim,Type>& v, fits::header& hdr) {
        fits::input_image img(filename);
        hdr = img.read_header();
        img.read(v);
    }

    template<std::size_t Dim, typename Type>
    void read(const std::string& filename, vec<Dim, Type>& v) {
        fits::input_image(filename).read(v);
    }

    template<std::size_t Dim = 2, typename Type = double>
    vec<Dim, Type> read(const std::string& filename) {
        vec<Dim, Type> v;
        read(filename, v);
        return v;
    }

    inline vec1s read_sectfits(const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            throw fits::exception("could not find file '"+filename+"'");
        }

        std::string dir = file::get_directory(filename);

        vec1s files;
        std::string line;
        while (std::getline(file, line)) {
            if (!line.empty() && line[0] != '#') {
                files.push_back(dir+line);
            }
        }

        if (files.empty()) {
            throw fits::exception("empty sectfits '"+filename+"'");
        }

        return files;
    }

    inline fits::header read_header(const std::string& filename) {
        return fits::generic_file(filename).read_header();
    }

    inline fits::header read_header_hdu(const std::string& filename, uint_t hdu) {
        return fits::generic_file(filename, hdu).read_header();
    }

    inline fits::header read_header_sectfits(const std::string& filename, uint_t sect) {
        if (ends_with(filename, ".sectfits")) {
            vec1s sects = read_sectfits(filename);
            phypp_check(sect < sects.size(), "no section ", sect, " in '", filename,
                "' (only ", sects.size()," available)");

            return read_header(sects[sect]);
        } else {
            return read_header(filename);
        }
    }

    inline fits::header read_header_sectfits_hdu(const std::string& filename,
        uint_t sect, uint_t hdu) {

        if (ends_with(filename, ".sectfits")) {
            vec1s sects = read_sectfits(filename);
            phypp_check(sect < sects.size(), "no section ", sect, " in '", filename,
                "' (only ", sects.size()," available)");

            return read_header_hdu(sects[sect], hdu);
        } else {
            return read_header_hdu(filename, hdu);
        }
    }

    // Write an image in a FITS file
    template<std::size_t Dim, typename Type>
    void write(const std::string& filename, const vec<Dim,Type>& v, const fits::header& hdr) {
        fits::output_image img(filename);
        img.write(v);
        img.write_header(hdr);
    }

    template<std::size_t Dim, typename Type>
    void write(const std::string& filename, const vec<Dim,Type>& v) {
        fits::output_image(filename).write(v);
    }

    // Write an image in a FITS file
    template<std::size_t Dim, typename Type>
    void update_hdu(const std::string& filename, const vec<Dim,Type>& v, uint_t hdu) {
        fits::image(filename, hdu).update(v);
    }

    // Read information about the columns of a FITS table
    inline vec<1,column_info> read_table_columns(const std::string& filename) {
        return fits::input_table(filename).read_columns_info();
    }
    inline vec<1,column_info> read_table_columns(const std::string& filename, uint_t hdu) {
        fits::input_table itbl(filename);
        itbl.reach_hdu(hdu);
        return itbl.read_columns_info();
    }

    // Read several columns in a FITS file.
    template<typename ... Args>
    void read_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fits::input_table(filename).read_columns(name, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void read_table_loose(const std::string& filename, const std::string& name, Args&& ... args) {
        fits::input_table(filename).read_columns(fits::missing, name, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void read_table(const std::string& filename, impl::ascii_impl::macroed_t,
        const std::string& names, Args&& ... args) {
        fits::input_table(filename).read_columns(impl::ascii_impl::macroed_t{},
            names, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void read_table_loose(const std::string& filename, impl::ascii_impl::macroed_t,
        const std::string& names, Args&& ... args) {
        fits::input_table(filename).read_columns(fits::missing, impl::ascii_impl::macroed_t{},
            names, std::forward<Args>(args)...);
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void read_table(const std::string& filename, T& t) {
        fits::input_table(filename).read_columns(t);
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void read_table_loose(const std::string& filename, T& t) {
        fits::input_table(filename).read_columns(fits::missing, t);
    }

    // Write several columns in a FITS file.
    template<typename ... Args>
    void write_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fits::output_table(filename).write_columns(name, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void write_table(const std::string& filename, impl::ascii_impl::macroed_t,
        const std::string& names, Args&& ... args) {
        fits::output_table(filename).write_columns(impl::ascii_impl::macroed_t{},
            names, std::forward<Args>(args)...);
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void write_table(const std::string& filename, const T& t) {
        fits::output_table(filename).write_columns(t);
    }

    // Append an new column to an existing FITS table
    // Note: if the column already exsits in the file, it will be overwritten

    template<typename ... Args>
    void update_table(const std::string& filename, const std::string& name, Args&& ... args) {
        fits::table(filename).update_columns(name, std::forward<Args>(args)...);
    }

    template<typename ... Args>
    void update_table(const std::string& filename, impl::ascii_impl::macroed_t,
        const std::string& names, Args&& ... args) {
        fits::table(filename).update_columns(impl::ascii_impl::macroed_t{},
            names, std::forward<Args>(args)...);
    }

    template<typename T, typename enable = typename std::enable_if<reflex::enabled<T>::value>::type>
    void update_table(const std::string& filename, const T& t) {
        fits::table(filename).update_columns(t);
    }
}
}

#endif
#endif

