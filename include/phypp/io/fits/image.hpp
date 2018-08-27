#ifndef PHYPP_IO_FITS_IMAGE_HPP
#define PHYPP_IO_FITS_IMAGE_HPP

#include "phypp/io/fits/base.hpp"

#ifndef NO_CFITSIO

namespace phypp {
namespace fits {
    // FITS input table (read only)
    class input_image : public virtual impl::fits_impl::file_base {
    public :
        input_image() :
            impl::fits_impl::file_base(impl::fits_impl::image_file, impl::fits_impl::read_only) {}

        explicit input_image(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::image_file, filename, impl::fits_impl::read_only) {}
        explicit input_image(const std::string& filename, uint_t hdu) :
            impl::fits_impl::file_base(impl::fits_impl::image_file, filename, impl::fits_impl::read_only) {
            reach_hdu(hdu);
        }

        input_image(input_image&&) noexcept = default;
        input_image(const input_image&) noexcept = delete;
        input_image& operator = (input_image&&) noexcept = delete;
        input_image& operator = (const input_image&&) noexcept = delete;

    protected:
        template<typename Type>
        void read_prep_(uint_t rdims, int& naxis, std::vector<long>& naxes, int& type) const {
            fits_get_img_dim(fptr_, &naxis, &status_);
            fits::phypp_check_cfitsio(status_, "could not read dimensions of HDU");
            phypp_check_fits(naxis == int(rdims), "FITS file has wrong number of dimensions "
                "(expected "+to_string(rdims)+", got "+to_string(naxis)+")");

            int bitpix;
            naxes.resize(naxis);
            fits_get_img_param(fptr_, naxis, &bitpix, &naxis, naxes.data(), &status_);
            fits::phypp_check_cfitsio(status_, "could not read image parameters of HDU");

            type = impl::fits_impl::bitpix_to_type(bitpix);
            phypp_check_fits(impl::fits_impl::traits<Type>::is_convertible(type), "wrong image type "
                "(expected "+pretty_type_t(Type)+", got "+impl::fits_impl::type_to_string_(type)+")");

            type = impl::fits_impl::traits<Type>::ttype;
        }

        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel) const {}

        template <typename ... Args>
        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel,
            impl::range_impl::full_range_t, const Args& ... args) const {

            fpixel[naxes.size()-1-idim] = 1;
            lpixel[naxes.size()-1-idim] = naxes[naxes.size()-1-idim];

            make_indices_(idim+1, naxes, fpixel, lpixel, args...);
        }

        template <typename ... Args>
        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel,
            impl::range_impl::left_range_t r, const Args& ... args) const {

            phypp_check_fits(r.last < uint_t(naxes[naxes.size()-1-idim]), "image subset goes outside of "
                "the image boundaries (axis "+to_string(idim)+": "+to_string(r.last)+" vs. "+
                to_string(naxes[naxes.size()-1-idim]));

            fpixel[naxes.size()-1-idim] = 1;
            lpixel[naxes.size()-1-idim] = r.last+1;

            make_indices_(idim+1, naxes, fpixel, lpixel, args...);
        }

        template <typename ... Args>
        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel,
            impl::range_impl::right_range_t r, const Args& ... args) const {

            fpixel[naxes.size()-1-idim] = r.first+1;
            lpixel[naxes.size()-1-idim] = naxes[naxes.size()-1-idim];

            make_indices_(idim+1, naxes, fpixel, lpixel, args...);
        }

        template <typename ... Args>
        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel,
            impl::range_impl::left_right_range_t r, const Args& ... args) const {

            phypp_check_fits(r.last < uint_t(naxes[naxes.size()-1-idim]), "image subset goes outside of "
                "the image boundaries (axis "+to_string(idim)+": "+to_string(r.last)+" vs. "+
                to_string(naxes[naxes.size()-1-idim]));

            fpixel[naxes.size()-1-idim] = r.first+1;
            lpixel[naxes.size()-1-idim] = r.last+1;

            make_indices_(idim+1, naxes, fpixel, lpixel, args...);
        }

        template <typename ... Args>
        void make_indices_(uint_t idim, const std::vector<long>& naxes,
            std::vector<long>& fpixel, std::vector<long>& lpixel,
            uint_t i, const Args& ... args) const {

            phypp_check_fits(i < uint_t(naxes[naxes.size()-1-idim]), "image subset goes outside of "
                "the image boundaries (axis "+to_string(idim)+": "+to_string(i)+" vs. "+
                to_string(naxes[naxes.size()-1-idim]));

            fpixel[naxes.size()-1-idim] = i+1;
            lpixel[naxes.size()-1-idim] = i+1;

            make_indices_(idim+1, naxes, fpixel, lpixel, args...);
        }

    public:
        template<std::size_t Dim, typename Type>
        void read(vec<Dim,Type>& v) const {
            check_is_open_();

            int naxis, type;
            std::vector<long> naxes;
            read_prep_<Type>(Dim, naxis, naxes, type);

            for (uint_t i : range(naxis)) {
                v.dims[i] = naxes[naxis-1-i];
            }

            v.resize();

            Type def = impl::fits_impl::traits<Type>::def();
            int anynul;
            fits_read_img(fptr_, type, 1, v.size(), &def, v.raw_data(), &anynul, &status_);
            fits::phypp_check_cfitsio(status_, "could not read image from HDU");
        }

        template<std::size_t Dim, typename Type, typename ... Args>
        void read_subset(vec<Dim,Type>& v, const Args& ... args) const {
            static_assert(Dim == sizeof...(Args), "incompatible subset and vector dimensions");

            check_is_open_();

            int naxis, type;
            std::vector<long> naxes;
            read_prep_<Type>(Dim, naxis, naxes, type);

            std::vector<long> fpixel(Dim), lpixel(Dim);
            make_indices_(0, naxes, fpixel, lpixel, args...);

            for (uint_t i : range(naxis)) {
                v.dims[i] = (lpixel[naxis-1-i]-fpixel[naxis-1-i])+1;
            }

            v.resize();

            Type def = impl::fits_impl::traits<Type>::def();
            int anynul;
            std::vector<long> inc(Dim, 1);
            fits_read_subset(fptr_, type, fpixel.data(), lpixel.data(), inc.data(),
                &def, v.raw_data(), &anynul, &status_);
            fits::phypp_check_cfitsio(status_, "could not read subset image from HDU");
        }

        template<typename Type = double>
        Type read_pixel(vec1u p) const {
            check_is_open_();

            int naxis, type;
            std::vector<long> naxes;
            read_prep_<Type>(p.size(), naxis, naxes, type);

            uint_t ppos = 1;
            uint_t pitch = 1;
            bool no_error = true;
            for (uint_t i : range(p)) {
                no_error = no_error && p.safe[naxis-1-i] < uint_t(naxes[i]);
                ppos += p.safe[naxis-1-i]*pitch;
                pitch *= naxes[i];
            }

            if (!no_error) {
                vec<1,long> nn; nn.data = naxes; nn.dims = p.dims;
                nn = reverse(nn);
                phypp_check_fits(no_error,
                    "FITS file has too small dimensions (reading pixel "+to_string(p)+
                    " in dims "+to_string(nn)+")");
            }

            type = impl::fits_impl::traits<Type>::ttype;

            Type val;
            Type def = impl::fits_impl::traits<Type>::def();
            int anynul;
            fits_read_img(fptr_, type, ppos, 1, &def, &val, &anynul, &status_);
            fits::phypp_check_cfitsio(status_, "could not pixel from HDU");

            return val;
        }
    };

    // Output FITS table (write only, overwrites existing files)
    class output_image : public virtual impl::fits_impl::file_base {
    public :
        output_image() :
            impl::fits_impl::file_base(impl::fits_impl::image_file, impl::fits_impl::write_only) {}

        explicit output_image(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::image_file, filename, impl::fits_impl::write_only) {}

        output_image(output_image&&) noexcept = default;
        output_image(const output_image&) = delete;
        output_image& operator = (output_image&&) noexcept = delete;
        output_image& operator = (const output_image&&) = delete;

    protected :

        template<std::size_t Dim, typename Type>
        void write_impl_(const vec<Dim,Type>& v) {
            fits_write_img(fptr_, impl::fits_impl::traits<Type>::ttype, 1, v.size(),
                const_cast<typename vec<Dim,Type>::dtype*>(v.raw_data()), &status_);
            fits::phypp_check_cfitsio(status_, "could not write image to HDU");
        }

    public :

        template<std::size_t Dim, typename Type>
        void write(const vec<Dim,Type>& v) {
            check_is_open_();

            std::array<long,Dim> naxes;
            for (uint_t i : range(Dim)) {
                naxes[i] = v.dims[Dim-1-i];
            }

            if (hdu_count() > 0) {
                // Check current HDU is empty or null
                auto type = hdu_type();
                phypp_check(type == fits::null_hdu || type == fits::empty_hdu,
                    "cannot write image, there is already data in this HDU");

                // Then resize it
                fits_resize_img(fptr_, impl::fits_impl::traits<meta::rtype_t<Type>>::image_type, Dim,
                    naxes.data(), &status_);
                fits::phypp_check_cfitsio(status_, "could not create image HDU");
            } else {
                // No HDU yet, just create the image in the primary array
                fits_insert_img(fptr_, impl::fits_impl::traits<meta::rtype_t<Type>>::image_type, Dim,
                    naxes.data(), &status_);
                fits::phypp_check_cfitsio(status_, "could not create image HDU");
            }

            // Finally write the data
            write_impl_(v.concretise());
        }

        void write_empty() {
            check_is_open_();

            long naxes = 0;

            if (hdu_count() > 0) {
                // Check current HDU is empty or null
                auto type = hdu_type();
                phypp_check(type == fits::null_hdu || type == fits::empty_hdu,
                    "cannot write image, there is already data in this HDU");

                // Then resize it
                fits_resize_img(fptr_, impl::fits_impl::traits<float>::image_type, 0, &naxes, &status_);
                fits::phypp_check_cfitsio(status_, "could not create image HDU");
            } else {
                // And create the new HDU
                fits_insert_img(fptr_, impl::fits_impl::traits<float>::image_type, 0, &naxes, &status_);
                fits::phypp_check_cfitsio(status_, "could not create empty image HDU");
            }

            // Update internal structures, because CFITSIO won't do that by itself
            fits_set_hdustruc(fptr_, &status_);
            fits::phypp_check_cfitsio(status_, "could not create empty image HDU");
        }
    };

    // Input/output FITS table (read & write, modifies existing files)
    class image : public output_image, public input_image {
    public :
        image() :
            impl::fits_impl::file_base(impl::fits_impl::image_file, impl::fits_impl::read_write),
            output_image(), input_image() {}

        explicit image(const std::string& filename) :
            impl::fits_impl::file_base(impl::fits_impl::image_file, filename, impl::fits_impl::read_write),
            output_image(filename), input_image(filename) {}

        explicit image(const std::string& filename, uint_t hdu) :
            impl::fits_impl::file_base(impl::fits_impl::image_file, filename, impl::fits_impl::read_write),
            output_image(filename), input_image(filename) {
            reach_hdu(hdu);
        }

        image(image&&) noexcept = default;
        image(const image&) = delete;
        image& operator = (image&&) noexcept = delete;
        image& operator = (const image&&) = delete;

        template<std::size_t Dim, typename Type>
        void update(const vec<Dim,Type>& v) {
            check_is_open_();

            int naxis;
            fits_get_img_dim(fptr_, &naxis, &status_);
            fits::phypp_check_cfitsio(status_, "could not read dimensions of HDU");
            phypp_check_fits(naxis == Dim, "FITS file has wrong number of dimensions "
                "(expected "+to_string(Dim)+", got "+to_string(naxis)+")");

            int bitpix;
            std::vector<long> naxes(naxis);
            fits_get_img_param(fptr_, naxis, &bitpix, &naxis, naxes.data(), &status_);
            fits::phypp_check_cfitsio(status_, "could not read image parameters of HDU");

            int type = impl::fits_impl::bitpix_to_type(bitpix);
            phypp_check_fits(impl::fits_impl::traits<meta::rtype_t<Type>>::is_convertible(type), "wrong image type "
                "(expected "+pretty_type_t(meta::rtype_t<Type>)+", got "+impl::fits_impl::type_to_string_(type)+")");

            std::array<uint_t,Dim> d;
            for (uint_t i : range(Dim)) {
                d[i] = naxes[Dim-1-i];
            }

            phypp_check_fits(v.dims == d, "incompatible array dimensions ("+to_string(v.dims)+
                " vs. "+to_string(d)+")");

            write_impl_(v.concretise());
        }
    };
}
}

#endif
#endif
