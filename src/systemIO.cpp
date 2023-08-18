#include"systemIO.h"


void SystemIO::dump(FILE* fp, int step, Atoms& atoms, SimBox& simbox){
  std::fprintf(fp,"ITEM: TIMESTEP\n");
  std::fprintf(fp,"%d\n",step);
  std::fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  std::fprintf(fp,"%lu\n",atoms.size());
  std::fprintf(fp,"ITEM: BOX BOUNDS xx yy zz\n");
  std::fprintf(fp,"0.0 %f\n",simbox.L[0]);
  std::fprintf(fp,"0.0 %f\n",simbox.L[1]);
  std::fprintf(fp,"0.0 %f\n",simbox.L[2]);
  std::fprintf(fp,"ITEM: ATOMS id type x y z vx vy vz diameter\n");
  for(int i = 0; i < atoms.size(); i++){ 
    Atom& a = atoms[i];
    std::fprintf(fp,"%d %d %f %f %f %f %f %f %f\n"
      ,i,i
      ,a.q[0],a.q[1],a.q[2]
      ,a.p[0],a.p[1],a.p[2],a.a);          
  }
  return;
}

void SystemIO::load_state(std::string fname, Atoms& atoms, SimBox& simbox){
  FILE* fp = std::fopen(fname.c_str(),"r");
  char dummy[64];
  std::fprintf(stdout,"# reading %s\n",fname.c_str());

  int n,m,k; double f;
  std::fgets(dummy,64,fp); // std::fprintf(fp,"ITEM: TIMESTEP\n");
  std::fscanf(fp,"%d\n",&n);// std::fprintf(fp,"%d\n",step);
  std::fgets(dummy,64,fp); // std::fprintf(fp,"ITEM: NUMBER OF ATOMS\n");
  std::fscanf(fp,"%d\n",&n); // std::fprintf(fp,"%lu\n",atoms.size());
  std::fgets(dummy,64,fp); // std::fprintf(fp,"ITEM: BOX BOUNDS xx yy zz\n");
  std::fscanf(fp,"%lf %lf\n",&f,&simbox.L[0]); // std::fprintf(fp,"0.0 %f\n",simbox.L[0]);
  std::fscanf(fp,"%lf %lf\n",&f,&simbox.L[1]); // std::fprintf(fp,"0.0 %f\n",simbox.L[1]);
  std::fscanf(fp,"%lf %lf\n",&f,&simbox.L[2]); // std::fprintf(fp,"0.0 %f\n",simbox.L[2]);
  std::fgets(dummy,64,fp); // std::fprintf(fp,"ITEM: ATOMS id type x y z vx vy vz diameter\n");
  atoms.resize(n);
  for(int i = 0; i < n; i++){ 
    double qx,qy,qz;
    double px,py,pz;
    std::fscanf(fp,"%d %d %lf %lf %lf %lf %lf %lf %lf\n"
      ,&m,&k
      ,&qx,&qy,&qz
      ,&px,&py,&pz,&f);          
    Atom& a = atoms[i];
    a.q = {qx,qy,qz};
    a.p = {px,py,pz};
  }
  std::fclose(fp);

  std::fprintf(stdout,"# atom size = %d\n",n);
  std::fprintf(stdout,"# box size = {%f, %f, %f}\n",simbox.L[0],simbox.L[1],simbox.L[2]);
  std::fprintf(stdout,"# box volume = %f\n",simbox.size());
  std::fprintf(stdout,"# number density = %f\n",n/simbox.size());
  std::fprintf(stdout,"# finish loading\n");
  return;
}

