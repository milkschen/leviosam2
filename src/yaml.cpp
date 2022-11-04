/*
 * yaml.hpp
 *
 * Simple functions to handle YAML config files.
 *
 * The core library is rapidyaml: https://github.com/biojppm/rapidyaml
 *
 * Author: Nae-Chyun Chen
 * Dept. of Computer Science, Johns Hopkins University
 *
 * Distributed under the MIT license
 * https://github.com/milkschen/leviosam2
 */
#include "yaml.hpp"

#include <iostream>
#define RYML_SINGLE_HDR_DEFINE_NOW
#include "rapidyaml.hpp"

namespace Yaml {

std::string read_yaml_raw(std::ifstream& f) {
    return std::string((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
}

std::string read_yaml_raw(const std::string& fn) {
    std::ifstream f(fn);
    return std::string((std::istreambuf_iterator<char>(f)),
                       std::istreambuf_iterator<char>());
}

ryml::Tree parse_yaml(const std::string& filename) {
    const char* c = filename.c_str();
    std::string contents = file_get_contents<std::string>(c);
    ryml::Tree tree = ryml::parse_in_arena(ryml::to_csubstr(contents));
    return tree;
}

};  // namespace Yaml

// int main(int argc, char** argv) {
//     std::string filename = "pacbio.yaml";
//     ryml::Tree tree = parse_yaml(filename);
//     std::cout << tree["q"].val() << "\n";
//     std::cout << tree["q2"].val() << "\n";
//     std::cout << tree["flag"].val() << "\n";
//     std::cout << tree["alignment"]["t"].val() << "\n";
// }
