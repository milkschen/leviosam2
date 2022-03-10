#include <iostream>
#include "yaml.hpp"
#define RYML_SINGLE_HDR_DEFINE_NOW
#include "rapidyaml.hpp"

namespace Yaml {

template<class CharContainer> CharContainer file_get_contents(const char *filename);
template<class CharContainer> size_t        file_get_contents(const char *filename, CharContainer *v);

// helper functions for sample_parse_file()

C4_SUPPRESS_WARNING_MSVC_WITH_PUSH(4996) // fopen: this function or variable may be unsafe
/** load a file from disk and return a newly created CharContainer */
template<class CharContainer>
size_t file_get_contents(const std::string& filename, CharContainer *v)
{
    ::FILE *fp = ::fopen(filename.data(), "rb");
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
CharContainer file_get_contents(const std::string& filename)
{
    CharContainer cc;
    file_get_contents(filename, &cc);
    return cc;
}

std::string read_yaml_raw(std::ifstream &f) {
    return std::string(
        (std::istreambuf_iterator<char>(f)),
        std::istreambuf_iterator<char>());
}

std::string read_yaml_raw(const std::string& fn) {
    std::ifstream f(fn);
    return std::string(
        (std::istreambuf_iterator<char>(f)),
        std::istreambuf_iterator<char>());
}

ryml::Tree parse_yaml(const std::string& filename) {
    std::string contents = file_get_contents<std::string>(filename);
    ryml::Tree tree = ryml::parse_in_arena(ryml::to_csubstr(contents));
    return tree;
}

};

// int main(int argc, char** argv) {
//     std::string filename = "pacbio.yaml";
//     ryml::Tree tree = parse_yaml(filename);
//     std::cout << tree["q"].val() << "\n";
//     std::cout << tree["q2"].val() << "\n";
//     std::cout << tree["flag"].val() << "\n";
//     std::cout << tree["alignment"]["t"].val() << "\n";
// }
