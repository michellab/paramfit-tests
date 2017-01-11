About
=====

In this directory you can find all the script and files to test Paramfit with
capped trialanine, starting with a quantum scan at HF level of theory.

To replicate ff99SB and study the stability/consistency of Paramfit I made 4 tests:

- original:
            all the folders under this name has a starting topology with the
            correct ff99SB amplitude for phi and psi dihedral angles.
            In this way, when Paramfit tries to fit quantum and molecular
            mechanical energies the starting point for phi and psi is exactly
            the same as in the forcefield

- multiplicity_3 :
            rather than having 4 multiplicities as in ff99SB we try to have
            only 3 multiplicities. This may sound wrong, but it seems to be what
            Paramfit authors have done in the paper.
            The initial amplitude for phi and psi is set to 0.0 kcal/mol

- multiplicity_3_amplitude_1 :
            as above but the initial amplitude is set to 1.0
            kcal/mol . This is what authors of Paramfit papars have done, but the
            initial value of the amplitude SHOULD NOT influence the final fit

- multiplicity_3_amplitude_2 :
            as above but with initial amplitude of 2.0 kcal/mol
            In this way we can check  the consistency of Paramfit

If the initial choice does not influence the final results ( and it should since
it is a LLS problem) we should end with the same value for all the multiplicity_*
topologies.

Folders
=====

1_inputfiles:

with dihedral_scan.py script input and submission files for Gaussian optimization
are created

2_quantum_extract:

here we extract the quantum energies form quantum output files, we calculate the
respective molecular mechanical energy with Paramfit and we create one mdcrd file,
which contains all the optimized quantum structures (necessary for Paramfit)

3__fit_offset:

Once the molecular mechanics and quantum energies are retrieved we need to calculate
the offset to procede with the fitting. The offset is calculated with Paramfit

4_fit_amplitude:

This is the real fitting. Here we try to fit each of the starting topology we have
with the quantum output. In this case we proceed normally with amber vs quantum
Finally, a comparison between quantum and mm adiabatic maps is done

5_fit_weight:

Here I put weights on the quantum energies. As written in Paramfit paper, you can
give a low weight (0.0) to the highest quantum energies, in order to fit better
the minima (weight = 10.0). Furthermore, I gave weight also to the intermediate
quantum energy from 3.0 to 5.0.

6_fith_weight_max:

Here I weighted 0.0 the max energy values and 10.0 the minima and no weights to
the intermediate.
