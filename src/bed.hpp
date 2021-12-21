#ifndef BED_HPP
#define BED_HPP

#include <iostream>
#include <unordered_map>
#include "IITree.h"
#include "robin_hood.h"

namespace BedUtils {

using BedMap = robin_hood::unordered_map<std::string, IITree<std::size_t, bool>>;

class Bed {
    public:
        Bed();
        Bed(const std::string& fn);
        void init(const std::string& fn);

        int index();
        bool add_interval(const std::string &line);
        bool intersect(
            const std::string &contig,
            const size_t &pos1, const size_t &pos2);
        bool intersect(const std::string &contig, const size_t &pos);
        BedMap get_intervals();
        std::string get_fn();

    private:
        BedMap intervals;
        std::string bed_fn = "";
        bool is_valid = false;

}; // Bed class

}; // namespace

#endif
