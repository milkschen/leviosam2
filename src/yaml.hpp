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
#include "rapidyaml.h"

namespace Yaml {

ryml::Tree parse_yaml(const std::string& filename);

};

#endif

