#ifndef PHYPP_FITS_IMAGE_HPP
#define PHYPP_FITS_IMAGE_HPP

#include "phypp/fits_base.hpp"

namespace fits {
    // FITS input table (read only)
    class input_image : public virtual impl::file_base {
    public :
        explicit input_image(const std::string& filename) :
            impl::file_base(impl::image_file, filename, impl::read_only) {}
        explicit input_image(const std::string& filename, uint_t hdu) :
            impl::file_base(impl::image_file, filename, impl::read_only) {
            reach_hdu(hdu);
        }

        input_image(input_image&&) = default;
        input_image(const input_image&) = delete;
        input_image& operator = (input_image&&) = delete;
        input_image& operator = (const input_image&&) = delete;

        template<std::size_t Dim, typename Type>
        void read(vec<Dim,Type>& v) const {
            status_ = 0;

            int naxis;
            fits_get_img_dim(fptr_, &naxis, &status_);
            phypp_check_fits(naxis == Dim, "FITS file has wrong number of dimensions "
                "(expected "+strn(Dim)+", got "+strn(naxis)+")");

            int bitpix;
            std::vector<long> naxes(naxis);
            fits_get_img_param(fptr_, naxis, &bitpix, &naxis, naxes.data(), &status_);

            int type = bitpix_to_type(bitpix);
            phypp_check_fits(traits<Type>::is_convertible(type), "wrong image type "
                "(expected "+pretty_type_t(Type)+", got "+type_to_string_(type)+")");

            type = traits<Type>::ttype;

            for (uint_t i : range(naxis)) {
                v.dims[i] = naxes[naxis-1-i];
            }

            v.resize();

            Type def = traits<Type>::def();
            int anynul;
            fits_read_img(fptr_, type, 1, v.size(), &def, v.data.data(), &anynul, &status_);
        }
    };

    // Output FITS table (write only, overwrites existing files)
    class output_image : public virtual impl::file_base {
    protected :

        explicit output_image(const std::string& filename, impl::readwrite_tag_t) :
            impl::file_base(impl::image_file, filename, impl::write_only) {}

    public :

        explicit output_image(const std::string& filename) :
            impl::file_base(impl::image_file, filename, impl::write_only) {}

        output_image(output_image&&) = default;
        output_image(const output_image&) = delete;
        output_image& operator = (output_image&&) = delete;
        output_image& operator = (const output_image&&) = delete;

    protected :

        template<std::size_t Dim, typename Type>
        void write_impl_(const vec<Dim,Type>& v) {
            fits_write_img(fptr_, traits<Type>::ttype, 1, v.size(),
                const_cast<typename vec<Dim,Type>::dtype*>(v.data.data()), &status_);
        }

    public :

        template<std::size_t Dim, typename Type>
        void write(const vec<Dim,Type>& v) {
            status_ = 0;

            std::array<long,Dim> naxes;
            for (uint_t i : range(Dim)) {
                naxes[i] = v.dims[Dim-1-i];
            }

            fits_create_img(fptr_, traits<rtype_t<Type>>::image_type, Dim,
                naxes.data(), &status_);

            write_impl_(v.concretise());
        }

        void write_empty() {
            status_ = 0;
            long naxes = 0;
            fits_create_img(fptr_, traits<float>::image_type, 0, &naxes, &status_);
        }
    };

    // Input/output FITS table (read & write, modifies existing files)
    class image : public output_image, public input_image {
    public :
        explicit image(const std::string& filename) :
            impl::file_base(impl::image_file, filename, impl::read_write),
            output_image(filename, impl::readwrite_tag), input_image(filename) {}

        explicit image(const std::string& filename, uint_t hdu) :
            impl::file_base(impl::image_file, filename, impl::read_write),
            output_image(filename, impl::readwrite_tag), input_image(filename) {
            reach_hdu(hdu);
        }

        image(image&&) = default;
        image(const image&) = delete;
        image& operator = (image&&) = delete;
        image& operator = (const image&&) = delete;

        template<std::size_t Dim, typename Type>
        void update(const vec<Dim,Type>& v) {
            status_ = 0;

            int naxis;
            fits_get_img_dim(fptr_, &naxis, &status_);
            phypp_check_fits(naxis == Dim, "FITS file has wrong number of dimensions "
                "(expected "+strn(Dim)+", got "+strn(naxis)+")");

            int bitpix;
            std::vector<long> naxes(naxis);
            fits_get_img_param(fptr_, naxis, &bitpix, &naxis, naxes.data(), &status_);

            int type = bitpix_to_type(bitpix);
            phypp_check_fits(traits<rtype_t<Type>>::is_convertible(type), "wrong image type "
                "(expected "+pretty_type_t(rtype_t<Type>)+", got "+type_to_string_(type)+")");

            std::array<uint_t,Dim> d;
            for (uint_t i : range(Dim)) {
                d[i] = naxes[Dim-1-i];
            }

            phypp_check_fits(v.dims == d, "incompatible array dimensions ("+strn(v.dims)+
                " vs. "+strn(d)+")");

            write_impl_(v.concretise());
        }
    };
}

#endif
