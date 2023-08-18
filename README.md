

# MD Simulator for 3D - Lennard Jones fluid

- 3D Bulk LJ system with Periodic BC
- Nose-hoover thermostat for NVT
- Formatting Lammps .dump file

## Usage

Build the source
~~~bash
make 
~~~

Move in working directory
~~~bash
cd hoge
~~~
	
Create initial configuration
~~~bash
python3 ../app/init.py 
~~~
	
run.sh simulate NVE or NVT system
~~~bash
./run.sh
~~~
	 
See the results
~~~bash
gnuplot ../app/stat.plt
open results.png
~~~
        
Ovito can display the trajectry
~~~bash
open -a Ovito LJrun.dump
~~~

 
