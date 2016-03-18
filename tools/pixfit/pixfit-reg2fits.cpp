#include <phypp.hpp>

bool read_ds9_region_circles(std::string file_name, vec2d& regs, vec1s& text, bool& physical, std::string color = "") {
    std::ifstream file(file_name);

    std::string global_color = "green";
    physical = false;
    std::string line;
    uint_t l = 0;
    while (std::getline(file, line)) {
        ++l;
        if (line.empty() || trim(line).empty() || trim(line)[0] == '#') continue;

        if (start_with(line, "global")) {
            std::string key = "color=";
            auto pos = line.find(key);
            if (pos != line.npos) {
                pos += key.size();
                global_color = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
            }
            continue;
        }

        auto spos = line.find_first_of('(');
        if (spos == line.npos) {
            if (trim(line) == "fk5") physical = false;
            if (trim(line) == "physical") physical = true;
            continue;
        }

        std::string type = trim(line.substr(0, spos));
        if (type != "circle") continue;

        auto epos = line.find_first_of(')', spos+1);
        std::string targs = line.substr(spos+1, epos-(spos+1));
        vec1s args = split(targs, ",");
        if (args.size() != 3) {
            error(file_name, ":", l, ": ",
                "ill formed 'circle' line, expecting 3 arguments, got ", args.size());
            return false;
        }

        double ra, dec, rad;
        args = trim(args);
        if (args[0].find_first_of(':') != args[0].npos) {
            if (!sex2deg(args[0], args[1], ra, dec)) {
                error(file_name, ":", l, ": ",
                    "could not convert sexagesimal coordinates to degrees");
                return false;
            }

            if (!end_with(args[2], "\"")) {
                error(file_name, ":", l, ": expected radius in arcsec");
                return false;
            }
        } else {
            if (!from_string(args[0], ra) || !from_string(args[1], dec)) {
                error(file_name, ":", l, ": ",
                    "could not read coordinates to ", (physical ? "(x,y)" : "degrees"));
                return false;
            }
        }

        if (physical) {
            if (!from_string(args[2], rad)) {
                error(file_name, ":", l, ": could not read radius in pixels");
                return false;
            }
        } else {
            args[2] = erase_end(args[2], "\"");
            if (!from_string(args[2], rad)) {
                error(file_name, ":", l, ": could not read radius in arcsec");
                return false;
            }
        }

        if (!color.empty()) {
            std::string rcol = global_color;
            spos = line.find_first_of('#', epos);
            if (spos != line.npos) {
                std::string key = "color=";
                auto pos = line.find(key, spos+1);
                if (pos != line.npos) {
                    pos += key.size();
                    rcol = trim(line.substr(pos, line.find_first_of(" \t", pos)-pos));
                }
            }

            if (rcol != color) continue;
        }

        std::string txt;
        spos = line.find_first_of('#', epos);
        if (spos != line.npos) {
            std::string key = "text=";
            auto pos = line.find(key, spos+1);
            if (pos != line.npos) {
                pos += key.size()+1;
                txt = line.substr(pos, line.find_first_of("}", pos)-pos);
            }
        }

        text.push_back(txt);
        append<0>(regs, vec2d{{ra, dec, rad}});
    }

    return true;
}

int main(int argc, char* argv[]) {
    vec2d regs;
    vec1s text;

    bool physical = false;
    if (!read_ds9_region_circles(argv[1], regs, text, physical)) {
        return 1;
    }

    if (physical) {
        error("must be WCS coordinates");
        return 1;
    }

    vec1u id;
    vec1f z, m;
    for (std::string txt : text) {
        vec1s spl = split(txt,",");
        if (spl.size() != 3) {
            error("need three elements: ID, z and M*");
            return 1;
        }

        uint_t tid;
        float tz, tm;
        if (!from_string(spl[0], tid) || !from_string(spl[1], tz) || !from_string(spl[2], tm)) {
            error("could not convert ID, z and M* into numbers");
            return 1;
        }

        id.push_back(tid);
        z.push_back(tz);
        m.push_back(tm);
    }

    fits::write_table(file::remove_extension(argv[1])+".fits",
        "id", id, "ra", vec1d{regs(_,0)}, "dec", vec1d{regs(_,1)}, "z", z, "m", m
    );

    return 0;
}
