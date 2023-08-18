#pragma once

#include"define.h"
#include<vector>


struct Atom {
  static constexpr double m=1.0; // mass
  static constexpr double a=1.0; // diameter
  vec q; // position
  vec p; // momentum
  vec F; // force
};
using Atoms = std::vector<Atom>;
