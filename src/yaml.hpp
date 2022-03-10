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

template<class CharContainer> CharContainer file_get_contents(const char *filename);
template<class CharContainer> size_t        file_get_contents(const char *filename, CharContainer *v);

template<class CharContainer>
CharContainer file_get_contents(const char* filename);

// helper functions for sample_parse_file()

// C4_SUPPRESS_WARNING_MSVC_WITH_PUSH(4996) // fopen: this function or variable may be unsafe
/** load a file from disk and return a newly created CharContainer */
template<class CharContainer>
size_t file_get_contents(const char* filename, CharContainer *v)
{
    ::FILE *fp = ::fopen(filename, "rb");
    C4_CHECK_MSG(fp != nullptr, "[E::yaml::file_get_contents] Could not open file %s", filename);
    ::fseek(fp, 0, SEEK_END);
    long sz = ::ftell(fp);
    v->resize(static_cast<typename CharContainer::size_type>(sz));
    if(sz)
    {
        ::rewind(fp);
        size_t ret = ::fread(&(*v)[0], 1, v->size(), fp);
        C4_CHECK(ret == (size_t)sz);
    }
    ::fclose(fp);
    return v->size();
}

/** load a file from disk into an existing CharContainer */
template<class CharContainer>
CharContainer file_get_contents(const char* filename)
{
    CharContainer cc;
    file_get_contents(filename, &cc);
    return cc;
}

// std::string read_yaml_raw(std::ifstream& f);
std::string read_yaml_raw(std::ifstream &f);
// std::string read_yaml_raw(const std::string& fn);
ryml::Tree parse_yaml(const char* filename);

};

#endif

