#include <phypp.hpp>

void print_help();

bool convolve(int argc, char* argv[]);
bool multiply(int argc, char* argv[]);

void print_convolve_help();
void print_multiply_help();

int main(int argc, char* argv[]) {
    if (argc < 3) {
        print_help();
        return 0;
    }

    std::string op = tolower(argv[1]);

    if (op == "help") {
        op = tolower(argv[2]);
        if (op == "convolve") {
            print_convolve_help();
        } else if (op == "multiply") {
            print_multiply_help();
        } else {
            error("unknown operation '", op, "'");
        }
    } else {
        if (op == "convolve") {
            convolve(argc-1, argv+1);
        } else if (op == "multiply") {
            multiply(argc-1, argv+1);
        } else {
            error("unknown operation '", op, "'");
        }
    }

    return 0;
}

void print_convolve_help() {
    using namespace format;

    paragraph("The program will convolve the first image (map) with the second image (kernel).");

    header("List of available command line options:");
    bullet("normalize", "[flag] normalize kernel to unit integral before convolution");
    bullet("help", "[flag] print this text");
    print("");
}

bool convolve(int argc, char* argv[]) {
    if (argc < 4) {
        error("usage: imgtool convolve [map.fits] [kernel.fits] [output.fits]");
        return false;
    }

    bool help = false;
    bool normalize = false;
    read_args(argc-3, argv+3, arg_list(help, normalize));

    if (help) {
        print_convolve_help();
        return true;
    }

    vec2d map = fits::read(argv[1]);
    vec2d kernel = fits::read(argv[2]);
    if (normalize) {
        kernel /= total(kernel);
    }

    map = convolve2d(map, kernel);

    file::mkdir(file::get_directory(argv[3]));
    fits::write(argv[3], map);

    return true;
}

void print_multiply_help() {
    using namespace format;

    paragraph("The program will multiply the first image by the provided value.");

    header("List of available command line options:");
    bullet("help", "[flag] print this text");
    print("");
}

bool multiply(int argc, char* argv[]) {
    if (argc < 4) {
        error("usage: imgtool multiply [img.fits] [value] [output.fits]");
        return false;
    }

    bool help = false;
    read_args(argc-3, argv+3, arg_list(help));

    if (help) {
        print_multiply_help();
        return true;
    }

    vec2d map = fits::read(argv[1]);
    double value;
    if (!from_string(argv[2], value)) {
        error("could not parse multiply value '", argv[2], "'");
        return false;
    }

    map *= value;

    file::mkdir(file::get_directory(argv[3]));
    fits::write(argv[3], map);

    return true;
}

void print_help() {
    using namespace format;

    print("imgtool v1.0");
    header("Usage: imgtool operation [options]");
    header("Available operations:");
    bullet("convolve", "convolve a 2D kernel to a map");
    print("");
    header("To learn more about each operations and see the list of avilable options, run "
        "'fitstool operation help'.");

    print("");
    paragraph("Copyright (c) 2014 C. Schreiber (corentin.schreiber@cea.fr)");

    paragraph("This software is provided 'as-is', without any express or implied warranty. In no "
        "event will the authors be held liable for any damages arising from the use of this "
        "software.");

    paragraph("Permission is granted to anyone to use this software for any purpose, including "
        "commercial applications, and to alter it and redistribute it freely, subject to the "
        "following restrictions:");

    bullet("1", "The origin of this software must not be misrepresented; you must not claim that "
        "you wrote the original software. If you use this software in a product, an acknowledgment "
        "in the product documentation would be appreciated but is not required.");
    bullet("2", "Altered source versions must be plainly marked as such, and must not be "
        "misrepresented as being the original software.");
    bullet("3", "This notice may not be removed or altered from any source distribution.");

    print("");
}
