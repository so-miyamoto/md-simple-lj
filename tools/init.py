

# https://zenn.dev/kaityo256/articles/md_initial_condition 


import numpy as np
import yaml
import random


def get_lattice_number(L, rho):
    m = np.ceil((L**3 * rho / 4.0)**(1.0 / 3.0))
    return int(m)


def make_fcc_pure(L, rho):
    m = get_lattice_number(L, rho)
    a = L / m
    ha = a * 0.5
    atoms = []
    for i in range(m**3):
        ix = i % m
        iy = (i // m) % m
        iz = i // (m * m)
        x = ix * a
        y = iy * a
        z = iz * a
        atoms.append((x, y, z))
        atoms.append((x + ha, y + ha, z))
        atoms.append((x + ha, y, z + ha))
        atoms.append((x, y + ha, z + ha))
    return atoms


def make_fcc_defect(L, rho):
    atoms = make_fcc_pure(L, rho)
    n = int(rho * L**3)
    return random.sample(atoms, n)

def scaling_momentum(T,vel):
    Kmean = (0.5/N)*np.sum(vel*vel)
    Tobs = Kmean/1.5
    alpha = np.sqrt(T/Tobs)
    vel *= alpha
    return

with open("params.yml","r") as file:
    print("opened params.yml")
    params = yaml.safe_load(file)
    L   = float(params["system"]["box_size"])
    rho = float(params["system"]["density"])
    T   = float(params["system"]["temperature"])

pos = make_fcc_defect(L, rho)
N = len(pos)
print(f"atoms={N}")

vel = np.random.rand(N*3).reshape(N,3)-0.5
vel -= np.mean(vel,axis=0)
scaling_momentum(T,vel)

with open("init.dump","w") as file:
    print("ITEM: TIMESTEP",file=file)
    print(f"{0}",file=file)
    print("ITEM: NUMBER OF ATOMS",file=file)
    print(f"{N}",file=file)
    print("ITEM: BOX BOUNDS xx yy zz",file=file)
    print(f"0.0 {L}",file=file)
    print(f"0.0 {L}",file=file)
    print(f"0.0 {L}",file=file)
    print("ITEM: ATOMS id type x y z vx vy vz radius",file=file)
    for i in range(N):
        a = pos[i]
        v = vel[i]
        print(f"{i} {i} {a[0]} {a[1]} {a[2]} {v[0]} {v[1]} {v[2]} {0.5}",file=file)
print("write init.dump")

with open("config.md","w") as file:
    print(params["integrator"]["dt"],file=file)
    print(params["integrator"]["max_timestep"],file=file)
    print(params["integrator"]["obs_timestep"],file=file)
    print(params["system"]["temperature"],file=file)
    print(params["control"]["nose_hoover_ctime"],file=file)
    print(params["hyparameter"]["cutoff"],file=file)
    print(params["hyparameter"]["margine"],file=file)
print("write config.md")

with open("run.sh","w") as file:
    if params["ensemble"]=="nve":
        print("../bin/nve.out",file=file)
    elif params["ensemble"]=="nvt":
        print("../bin/nvt.out",file=file)
    else:
        print("Unknown ensemble!")
print("create run.sh")

