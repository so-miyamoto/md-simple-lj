
#include"simbox.h"

void SimBox::initialize(const double icutoff,const double imargine){
  cutoff = icutoff;
  cutoff2 = cutoff*cutoff;
  margine = imargine;
  return;
};

int SimBox::get_neighbors_idx(vec q){
  int ix = q[0]/len_cell;      
  int iy = q[1]/len_cell;      
  int iz = q[2]/len_cell; 
  return ix*num_cells*num_cells+iy*num_cells+iz;
}

void SimBox::initialkick_neighbor_list(Atoms& atoms){
  num_cells = (L[0]/(cutoff+margine)) + 1;
  if( num_cells <= 2 ){
    std::fprintf(stderr,"too small simulation box to construct pair lists in cells");
    exit(1);
  }
  len_cell = L[0]/num_cells;
  atoms_lists.resize(num_cells*num_cells*num_cells);
  lifetime = -1.0;
  update_addresses(0.0,atoms);
  return;
}

void SimBox::update_addresses(const double dq_max,Atoms& atoms){
  lifetime -= dq_max;
  if( lifetime < 0.0 ){
    // std::fprintf(stdout,"construct neighbor list\n");
    _register_atoms(atoms);
    lifetime = margine;
  } 
  return;
}

void SimBox::_register_atoms(Atoms& atoms){
  for(int i = 0; i < atoms_lists.size(); i++){
    atoms_lists[i].clear();
  }
  for(int i = 0; i < atoms.size(); i++){
    vec& q = atoms[i].q;
    int ix = q[0]/len_cell;      
    int iy = q[1]/len_cell;      
    int iz = q[2]/len_cell; 
    for(int ni=-1; ni<=1; ni++){
      for(int nj=-1; nj<=1; nj++){
        for(int nk=-1; nk<=1; nk++){
          _register_atom(i,ix+ni,iy+nj,iz+nk);      
        }
      }
    }
  }
  return;
}

void SimBox::_register_atom(const int i, const int ix, const int iy, const int iz){
  int idx = ((ix+num_cells)%num_cells)*num_cells*num_cells 
          + ((iy+num_cells)%num_cells)*num_cells 
          + ((iz+num_cells)%num_cells);
  atoms_lists[idx].emplace_back(i);     
  return;
}


void SimBox::perform_boundary_condition(Atoms& atoms){
  constexpr int x=0, y=1, z=2;
  for(int i = 0; i < atoms.size(); i++){
    vec& q = atoms[i].q;
    if( q[x] <  0.0  ) q[x] += L[x]; // left x
    if( q[y] <  0.0  ) q[y] += L[y]; // left y
    if( q[z] <  0.0  ) q[z] += L[z]; // left z

    if( q[x] >= L[x]  ) q[x] -= L[x]; // right x
    if( q[y] >= L[y]  ) q[y] -= L[y]; // right y
    if( q[z] >= L[z]  ) q[z] -= L[z]; // right z
  }
  return;
}

vec SimBox::relative_position(vec dq){
  constexpr int x=0, y=1, z=2;
  if( dq[x] <  -0.5*L[x]  ) dq[x] += L[x]; // x
  if( dq[y] <  -0.5*L[y]  ) dq[y] += L[y]; // y
  if( dq[z] <  -0.5*L[z]  ) dq[z] += L[z]; // z

  if( dq[x] >= +0.5*L[x]  ) dq[x] -= L[x]; // x
  if( dq[y] >= +0.5*L[y]  ) dq[y] -= L[y]; // y
  if( dq[z] >= +0.5*L[z]  ) dq[z] -= L[z]; // z
  return dq;
}
// ------------------------------------------------------------
