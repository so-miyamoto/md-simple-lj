#pragma once

#include"define.h"
#include"atoms.h"
#include<cstdio>

#include<vector>

struct SimBox {
  vec L; // [0.0, L[0]] x [0.0, L[1]] x [0.0, L[2]]

  double cutoff;
  double cutoff2;

  // cell address method
  double margine;
  double lifetime;
  int num_cells;
  double len_cell;
  std::vector<std::vector<int>> atoms_lists;

  inline double size(){
    return L[0]*L[1]*L[2];
  }
  void initialize(const double icutoff,const double imargine);
  int get_neighbors_idx(vec q);
  void initialkick_neighbor_list(Atoms& atoms);
  void update_addresses(const double dq_max,Atoms& atoms);
private:
  void _register_atoms(Atoms& atoms);
  void _register_atom(const int i, const int ix, const int iy, const int iz);

public:
  void perform_boundary_condition(Atoms& atoms);
  vec relative_position(vec dq);
};
// ------------------------------------------------------------
