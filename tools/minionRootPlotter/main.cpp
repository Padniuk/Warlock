/**
 * @file main.cpp
 * @brief CLI entry point for minionRootPlotter - see H5ToRootConverter.hpp for the conversion logic.
 */

#include "H5ToRootConverter.hpp"
#include <iostream>

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cout << "Usage: minionRootPlotter <input.h5> [output.root]" << std::endl;
        return 1;
    }

    std::string in = argv[1];
    std::string out = (argc > 2) ? argv[2] : in.substr(0, in.find_last_of('.')) + ".root";

    framework::H5ToRootConverter converter;
    if (converter.convert(in, out)) {
        return 0;
    }

    return 1;
}