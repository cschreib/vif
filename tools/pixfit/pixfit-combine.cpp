#include "pixfit-common.hpp"

int vif_main(int argc, char* argv[]) {
    std::string map_file;
    std::string incat;
    std::string outcat;

    read_args(argc, argv, arg_list(name(map_file, "maps"), incat, outcat));

    bool bad = false;
    if (map_file.empty()) {
        error("please provide the map list file in maps=...");
        bad = true;
    }
    if (incat.empty()) {
        error("please provide the input flux catalog in incat=...");
        bad = true;
    }
    if (outcat.empty()) {
        error("please provide the output flux catalog in outcat=...");
        bad = true;
    }

    if (bad) return 1;

    std::vector<map_info> maps;
    if (!read_maps(map_file, maps)) return 1;

    file::copy(incat, outcat);

    fits::table tbl(outcat);

    uint_t ngal; {
        vec1u id;
        tbl.read_column("id", id);
        ngal = id.size();
    }

    vec2f flux, flux_err;
    vec2u flux_group;
    vec2f flux_group_cov;
    vec1f group_flux, group_flux_err;
    vec1s bands;
    tbl.read_columns(fits::missing, "bands", bands, "flux", flux, "flux_err", flux_err);

    flux_group = replicate(npos, flux.dims);
    flux_group_cov.resize(flux.dims);

    for (auto& map : maps) {
        std::string band = map.band;
        std::string cat = map.band+"-fit.fits";

        bands.push_back(band);

        fits::input_table itbl(cat);

        // Merge individual flux measurements
        vec1f tflx, terr;
        itbl.read_columns("flux", tflx, "flux_err", terr);

        append<1>(flux,     reform(tflx, ngal, 1));
        append<1>(flux_err, reform(terr, ngal, 1));

        // Merge group membership
        vec1u group_aper_id;
        vec1f group_cov;

        itbl.read_columns(fits::missing, "group_cov", group_cov, "group_aper_id", group_aper_id);

        if (group_cov.empty()) {
            group_cov = replicate(0.0, ngal);
            group_aper_id = replicate(0, ngal);
        }

        vec1u idg = where(group_cov > 0);
        if (!idg.empty()) {
            // Adjust IDs to refer to the merged groups
            vec1u idb = complement(group_cov, idg);
            uint_t i0 = min(group_aper_id[idg]);
            group_aper_id[idg] += group_flux.size() - i0;
            group_aper_id[idb] = npos;
        }

        append<1>(flux_group,     reform(group_aper_id, ngal, 1));
        append<1>(flux_group_cov, reform(group_cov,     ngal, 1));

        // Merge flux groups
        vec1b gfit;
        vec1f tgflx, tgerr;

        itbl.read_columns(fits::missing, "group_fit", gfit,
            "group_flux", tgflx, "group_flux_err", tgerr);

        vec1u ida = where(!gfit);

        if (!ida.empty()) {
            append(group_flux,     tgflx[ida]);
            append(group_flux_err, tgerr[ida]);
        }
    }

    tbl.update_column("bands", bands);
    tbl.update_column("flux", flux);
    tbl.update_column("flux_err", flux_err);
    tbl.update_column("flux_group", flux_group);
    tbl.update_column("flux_group_cov", flux_group_cov);
    tbl.update_column("group_flux", group_flux);
    tbl.update_column("group_flux_err", group_flux_err);

    return 0;
}
