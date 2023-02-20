# Generator

A simple set of function to generate particles with given distributions for
particle tracking.
 ...

## How to run

Add the command ```import generatorTool.py``` in your script. you can then use various function to generate particle distribution. If you need to write openPMD format you also need to include the command ```import makeopenpmd```


## Directory structure

- [ ] sobol_lib.py: peudo-random generator based on the Sobol sequencing (this is a legacy script should be removed once the main scripts uses the sobol sequence available in ```scipy``` 
- [ ] gen_cath_r_multi_t_gaus_p_iso.py: generate a distribution at the cathode surface (z=0 and varying time)
- [ ] gen_baneDistrib.py: generate Bane's temporal distribution
- [ ] gen_gaussian-WarpX.ipynb: generate a Gaussian beam for input into WarpX

