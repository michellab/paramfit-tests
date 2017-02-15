About
=====

Here I am trying to develop :

1) A robust protocol for dihedral fitting - may be with Paramfit or mgdx

2) A robust script which use gradients techniques to dihedral fitting problem
which can be implemented (an well documented) in Sire

Folders
=======

alanine/ :
test with paramfit and my own script on alanine ff99SB
Here I want to parametrize Phi and Psi dihedral and see the adiabatic map
the adiabatic map will be compared with  ff99SB and quantum mechanics results
In particular, for the moment as QM calculation we are using
- HF/6-31G : how much could this simple parametrization be effective?
- b3lyp/6-31 G: something more advance
- b3lyp/6-31 G optimization + MP2/cc-pVDZ  single point calculation
For each of them (at the moment) I have compute the Phi and Psi amplitudes with
Paramfit, by keeping the same multiplicity as ff99SB
I am testing now the adiabatic surface - FE surface

paramfit/ :
Here I have done some  tests on capped-trialanine ( or tetraalanine) in order to
reprodue Paramfit paper results on fitting the Phi and Psi amplitudes

stochastic_gradient_descent:
I would  say "some test" what is important here is to take the code, make it general
and try to apply to alanine case (remember to clean the directory)

scripts:
test on offsets

mgdx_tutorial:
to fill very soon (hopefully) by testing mgdx
