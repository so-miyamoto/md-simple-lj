#pragma once

#include"atoms.h"
#include"simbox.h"

#include<cstdio>
#include<string>

struct SystemIO {

  void dump(FILE* fp, int step, Atoms& atoms, SimBox& simbox);

  void load_state(std::string fname, Atoms& atoms, SimBox& simbox);

};
