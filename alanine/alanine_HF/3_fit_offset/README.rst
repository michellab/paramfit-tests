Usage
=====

Here we have to compute the offset between quantum energies and  amber energies
This is done with paramfit through a job.in file

For each folder we have we have to execute:
```
paramfit -i fit_K.in -p ../2_quantum_extract/topology/fit.prmtop -c ../2_quantum_extract/topology/all.mdcrd  -q ../2_quantum_extract/quantum_energies.dat > fit_K.out
```

with
```
ipyhton notebook fit_K_amber.ipynb
```
the general agreement between quantum energy and amber energy can be seen

N.B. The topology we are using here is with Phi and Psi amplitued set to 0.0 
