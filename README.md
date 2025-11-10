Very basic instructions for running something (can't really remember what) on Hamilton 8 using this version of LARE, which is vanilla except fo replacing sdf filesystem with netcdf, and creating the 'run.py' python wrapper which sets up the variables to be read in at runtime. 

# To install

- Navigate to where you want it to live

```
git clone https://github.com/oekrice/lare3d_netcdf.git
module load intel
module load intempi
module load netcdf
```
Test compiling manually with 
`make`

To set up parameters (resolutions etc.) to be read-in at run time, use
```
python run.py 0
```
where here 0 is the 'run number' (will be used in the output files if you want to do multiple runs at once). Obviously more can be hard-coded if necessary.

I believe this is set up to create an initial condition for the Pariat jet model, which is done in Python in 'init_bfield.py', and then this raw field is read in to Lare proper.

```
To actually run:
mpiexec -n 1 ./bin/lare3d 0
```

where here 1 is the number of processors and 0 is the run number. Outputs are saved to './Data/' by default. The file 'read_data.py' contains a script which is set up to read these outputs and plot them in some manner which I forget. It looks like the dimensions will need to be changed manually with this particular version.

Hope that works!
