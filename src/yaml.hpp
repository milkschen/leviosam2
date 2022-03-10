/*
 * yaml.hpp
 *
 * Simple functions to handle YAML config files. 
 *
 * The core library is rapidyaml: https://github.com/biojppm/rapidyaml
 *
 * Authors: Nae-Chyun Chen
 *
 * Distributed under the MIT license
 * https://github.com/alshai/levioSAM
 */
#ifndef YAML_HPP
#define YAML_HPP

#include <string>
#include <fstream>
#include "rapidyaml.hpp"

namespace Yaml {

template<class CharContainer>
CharContainer file_get_contents(const std::string& filename);

// std::string read_yaml_raw(std::ifstream& f);
std::string read_yaml_raw(std::ifstream &f);
// std::string read_yaml_raw(const std::string& fn);
ryml::Tree parse_yaml(const std::string& filename);

};

#endif

