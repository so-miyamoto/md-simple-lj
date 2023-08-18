#include"observer.h"


Observer::Observer(){
  std::fprintf(stdout,"# observer initialized\n");
  _fp = std::fopen("stat.log","w");
  std::fprintf(_fp,"# V  K  Eall  T  P\n");
  std::fprintf(stdout,"# observer open a new file: stat.log\n");
}

Observer::~Observer(){
  std::fclose(_fp);
  std::fprintf(stdout,"# observer close a file: stat.log\n");
}

vec Observer::average_momentum(Atoms& atoms){
  vec pm = {0.0,0.0,0.0};
  for(auto &a: atoms){
    pm[0] += a.p[0];
    pm[1] += a.p[1];
    pm[2] += a.p[2];
  }
  pm[0] /= atoms.size();
  pm[1] /= atoms.size();
  pm[2] /= atoms.size();
  return pm;
}

double Observer::kinetic_energy(Atoms& atoms){
  double K = 0.0;
  for(auto &a: atoms){
    K += (a.p[0]*a.p[0] + a.p[1]*a.p[1] + a.p[2]*a.p[2])/(a.m);
  }
  return 0.5*K/atoms.size();
}

void Observer::watch(const int step, Atoms& atoms, SimBox& simbox){
  std::vector<double> obs;
  vec V_P_T = potential_energy(atoms,simbox);
  const double V = V_P_T[0];
  const double K = V_P_T[2]*1.5;
  const double P = V_P_T[1];
  const double T = V_P_T[2];
  std::fprintf(_fp,"%d %f %f %f %f %f\n",step,V,K,V+K,T,P);
  return;
}


vec Observer::potential_energy(Atoms& atoms, SimBox& simbox){
  constexpr int x=0, y=1, z=2;
  const double rc2 = simbox.cutoff2;
  const double rc6 = rc2*rc2*rc2;
  const double C0 = - 4.0*( 1.0 - rc6 ) / (rc6*rc6);

  const int N = atoms.size();
  const double Vol = simbox.size();
  const double density = N/Vol;
  const double T = temperature(atoms);

  double V = 0.0;
  double virial = 0.0;

  for(int i = 0; i < N; i++){
    Atom& ai = atoms[i];
    auto& jlists = simbox.atoms_lists[simbox.get_neighbors_idx(ai.q)];
    // for(int j = 0; j < N; j++){ 
    for(auto& j: jlists){ 
      if( !( i < j ) ) continue;
      Atom& aj = atoms[j];
      vec qij;
      qij[x] = ai.q[x] - aj.q[x];
      qij[y] = ai.q[y] - aj.q[y];
      qij[z] = ai.q[z] - aj.q[z];
      qij = simbox.relative_position(qij);
      double r2 = qij[x]*qij[x] + qij[y]*qij[y] + qij[z]*qij[z];
      if( r2 > rc2 ) continue;
      double r6 = r2*r2*r2;
      V += 4.0*( 1.0 - r6 ) / (r6*r6) + C0; 
      virial += (48.0 - 24.0*r6)/(r6*r6); 
    }
  }
  virial /= (3.0*Vol);
  return {V/N  , density*T + virial, T};
}

