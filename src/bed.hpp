#ifndef BED_HPP
#define BED_HPP

#include <iostream>
#include <unordered_map>
//#include "IntervalTree.h"
#include "IITree.h"
#include "robin_hood.h"

namespace BedUtils {

using BedMap = robin_hood::unordered_map<std::string, IITree<std::size_t, bool>>;

class Bed {
    public:
        Bed() {};
        Bed(const std::string &fn);
        // Bed(std::ifstream& in);

        int index();
        bool add_interval(const std::string &line);
        bool intersect(const std::string &contig, const size_t &pos);

        BedMap get_intervals();

    private:
        BedMap intervals;
        // robin_hood::unordered_map<std::string, IITree<std::size_t, bool>> intervals;

}; // Bed class

}; // namespace

#endif
