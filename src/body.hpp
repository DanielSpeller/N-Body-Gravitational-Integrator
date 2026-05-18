#pragma once

#include <array>
#include <string>

struct Body {
    std::string           name;
    double                mass;
    std::array<double, 3> pos;   // position (AU)
    std::array<double, 3> vel;   // velocity (AU/yr)
};
