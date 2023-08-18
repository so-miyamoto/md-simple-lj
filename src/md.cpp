#include"define.h"
#include"atoms.h"
#include"simbox.h"
#include"observer.h"
#include"systemIO.h"

#include<iostream>
#include<vector>
#include<cmath>
#include<chrono>

class MDsystem {

public:

  MDsystem(const double tau_nose){
    _zeta_nose = 0.0;
    _tau_nose  = tau_nose;
  }

  void run(const int max_step, const double dt, const int obs_step,const double aimed_temperature, Atoms& atoms, SimBox& simbox){
    Observer observer;
    SystemIO sysio;

    FILE* fp = std::fopen("LJrun.dump","w");
    simbox.perform_boundary_condition(atoms);
    simbox.initialkick_neighbor_list(atoms);

    for(int step=0; step<max_step; step++){

      double dq_max = _update_positions(atoms,dt); // q(t) -> q(t + dt/2)
      simbox.perform_boundary_condition(atoms);
      simbox.update_addresses(dq_max,atoms);

      _compute_force(atoms,simbox); // F(t+dt/2) from q(t+dt/2)

      _update_momentums(atoms,dt); // p(t) -> p(t + dt)

#ifdef NVT_ENSEMBLE
      _nose_hoover(atoms,observer.kinetic_energy(atoms)/1.5,aimed_temperature,dt);
#endif

      double dq_max_ = _update_positions(atoms,dt); // q(t + dt/2) -> q(t + dt)
      simbox.perform_boundary_condition(atoms);
      simbox.update_addresses(dq_max_,atoms);

      if( step % obs_step == 0 ){
        std::fprintf(stdout,"step=%d\n",step);
        observer.watch(step,atoms,simbox);
        sysio.dump(fp,step,atoms,simbox);
      }
    }
    std::fclose(fp);
    return;
  }

private:

  double _zeta_nose;
  double _tau_nose;

  void _nose_hoover(Atoms& atoms,const double temperature ,const double aimed_temperature, const double dt){
    _zeta_nose += (temperature - aimed_temperature)/(_tau_nose*_tau_nose)*dt;
    for(auto& a:atoms){
      a.p[0] -= a.p[0] * _zeta_nose*dt;
      a.p[1] -= a.p[1] * _zeta_nose*dt;
      a.p[2] -= a.p[2] * _zeta_nose*dt;
    }
    return;
  }

  void _compute_force(Atoms& atoms, SimBox& simbox){
    constexpr int x=0, y=1, z=2;    
    const int N = atoms.size();
    for(auto& a: atoms){ a.F = {0.0,0.0,0.0};}
    for(int i = 0; i < N; i++){
      Atom& ai = atoms[i];

      auto& jlists = simbox.atoms_lists[simbox.get_neighbors_idx(ai.q)];
      // for(int j = 0; j < N; j++){ 
      for(auto& j: jlists){ 
        if( !(i < j) ) continue;
        Atom& aj = atoms[j];
        vec qij;
        qij[x] = ai.q[x] - aj.q[x];
        qij[y] = ai.q[y] - aj.q[y];
        qij[z] = ai.q[z] - aj.q[z];
        qij = simbox.relative_position(qij);
        double r2 = qij[x]*qij[x] + qij[y]*qij[y] + qij[z]*qij[z];
        if( r2 > simbox.cutoff2 ) continue;
        double r6 = r2*r2*r2;
        double df = (48.0 - 24.0*r6) / (r6*r6*r2); 
        ai.F[x] += df*qij[x];
        ai.F[y] += df*qij[y];
        ai.F[z] += df*qij[z];
        aj.F[x] -= df*qij[x];
        aj.F[y] -= df*qij[y];
        aj.F[z] -= df*qij[z];
      }
    }
    return;
  }

  double _update_positions(Atoms& atoms, const double dt){
    constexpr int x=0, y=1, z=2;    
    double dq_max2 = 0.0;
    for(int i = 0; i < atoms.size(); i++){
      vec dq;
      dq[x] = (0.5*dt/atoms[i].m)*atoms[i].p[x];
      dq[y] = (0.5*dt/atoms[i].m)*atoms[i].p[y];
      dq[z] = (0.5*dt/atoms[i].m)*atoms[i].p[z];
      atoms[i].q[x] += dq[x];
      atoms[i].q[y] += dq[y];
      atoms[i].q[z] += dq[z];
      dq_max2 = std::max(dq_max2, dq[x]*dq[x]+dq[y]*dq[y]+dq[z]*dq[z]);
    }
    return std::sqrt(dq_max2);
  }
  void _update_momentums(Atoms& atoms, const double dt){
    constexpr int x=0, y=1, z=2;    
    for(int i = 0; i < atoms.size(); i++){
      atoms[i].p[x] += dt*atoms[i].F[x];
      atoms[i].p[y] += dt*atoms[i].F[y];
      atoms[i].p[z] += dt*atoms[i].F[z];
    }
    return;
  }

};






int main(){
  double cutoff = 2.5;
  double margine = 1.0;
  int max_step = 10000;
  double dt = 0.01;
  int obs_step = 50;
  double aimed_temperature=1.0;
  double tau_nose = 0.10;

  FILE* fp = std::fopen("config.md","r");
  std::fscanf(fp,"%lf\n",&dt);
  std::fscanf(fp,"%d\n",&max_step);
  std::fscanf(fp,"%d\n",&obs_step);
  std::fscanf(fp,"%lf\n",&aimed_temperature);
  std::fscanf(fp,"%lf\n",&tau_nose);
  std::fscanf(fp,"%lf\n",&cutoff);
  std::fscanf(fp,"%lf\n",&margine);
  std::fclose(fp);

  std::fprintf(stdout,"dt      =%lf\n",dt);
  std::fprintf(stdout,"max_step=%d\n",max_step);
  std::fprintf(stdout,"obs_step=%d\n",obs_step);
  std::fprintf(stdout,"target T=%lf\n",aimed_temperature);
  std::fprintf(stdout,"tau_nose=%lf\n",tau_nose);
  std::fprintf(stdout,"r_cutoff=%lf\n",cutoff);
  std::fprintf(stdout,"margine =%lf\n",margine);

  // initialize
  Atoms atoms;
  SimBox simbox;
  simbox.initialize(cutoff,margine);

  SystemIO sysio;
  sysio.load_state(std::string("init.dump"),atoms,simbox);

  // runnung
  std::chrono::system_clock::time_point  start, end;
  std::time_t time_stamp;
  start = std::chrono::system_clock::now();

  MDsystem md(tau_nose);
  md.run(max_step,dt,obs_step,aimed_temperature,atoms,simbox);

  end = std::chrono::system_clock::now(); 
  auto time = end - start;
  time_stamp = std::chrono::system_clock::to_time_t(start);
  std::cout << std::ctime(&time_stamp);
  auto msec = std::chrono::duration_cast<std::chrono::milliseconds>(time).count();
  std::cout << msec << " msec" << std::endl;
  std::fprintf(stdout,"# simulation completed :)\n");
  
  return 0;
}
