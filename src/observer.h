#pragma once

#include"atoms.h"
#include"simbox.h"
#include<cstdio>

struct Observer {

private:
  FILE* _fp;

public:
  Observer();
  ~Observer();

  void watch(const int step, Atoms& atoms, SimBox& simbox);

  vec average_momentum(Atoms& atoms);

  vec potential_energy(Atoms& atoms, SimBox& simbox);

  double kinetic_energy(Atoms& atoms);

  inline double temperature(Atoms& atoms){
    return kinetic_energy(atoms)/1.5;
  }

};
