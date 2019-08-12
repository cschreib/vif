#ifndef VIF_IO_FITS_FILE_HPP
#define VIF_IO_FITS_FILE_HPP

#include "vif/io/fits/image.hpp"
#include "vif/io/fits/table.hpp"

#ifndef NO_CFITSIO

namespace vif {
namespace fits {
    class input_file : public fits::input_image, public fits::input_table {
    public :
        explicit input_file(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::generic_file, filename, impl::fits_impl::read_only),
            input_image(filename), input_table(filename) {}

        input_file(input_file&&) noexcept = default;
    };

    class output_file : public fits::output_image, public fits::output_table {
    public :
        explicit output_file(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::generic_file, filename, impl::fits_impl::write_only),
            output_image(filename), output_table(filename) {}

        output_file(output_file&&) noexcept = default;
    };

    class file : public fits::image, public fits::table {
    public :
        explicit file(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::generic_file, filename, impl::fits_impl::read_write),
            image(filename), table(filename) {}

        file(file&&) noexcept = default;
    };
}
}

#endif
#endif
