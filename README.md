Hopefully some of this info answers some of the questions you might have about the code.  I've tried my best but it's been a while since I did this so don't assume that everything stated below is explained 100% accurrately.

### Process for running a coarse grain simulation for three different types of particle in a mixture:
---
1. Manually create a .pdb file for each of your particles, such as the `CG_Water.pdb` file
2. Once you have these made you then create a .inp file, such as `Example_Structure.inp` which contains the information about the number of each particles and the volume they take up (units are in nm).
3. Install Packmol for Julia and run `packmol_script.jl` with this file to generate an evenly packed mixture. `Example_Structure.inp` would generate `Example_Structure.
4. .pdb` provided all the individual .pdb files for each particle exists.
5. GROMACS requires an additional line in the .pdb structure file that specifies the dimensions of the volume (cubic) that packmol doesn't provide, which is line 5 in `Example_Structure.pdb`:
```
CRYST1   x(nm)  y(nm)  z(nm)  90.00  90.00  90.00 P 1           1
```
5. You'll also need a topology file (.top) which stores the information on how particles interact. See `Example_Structure.top`. **Note**: the `[ nonbond_params ]` data allows you to tweak interaction parameters between particular particles independantly and overwrites anything under `[ atomtypes ]`. Some other lines might be redundant but this seemed to work for me.
6. Then you'll need the .mdp files which will dictate how the mixture reaches equilbrium. Examples of these with comments can be found in the `.mdp files` directory. A lot of the lines in the example files are redundant.
7. Then organise these however you want so that the file paths coincide with those in `HPC Scripts/eAP-eAB.pbs`. Change the `# Define file and folder names
` section to match your setup. The only reason for the particular script name was because I was varying epsilons between particles at the time. 
8. Running the .pbs file on the HPC should output energies (corresponding to the integers on line 61 and 84), pressure, volume and density during both equilbriation and the production run, as well as the radial distribution functions between a selected particle and a reference particle during production run.

9. ### Process for running finding free energy of solvation of a mixture of two types of coarse grain particles:
---
The idea is taken from: https://tutorials.gromacs.org/docs/free-energy-of-solvation.html
1. Repeat steps 1-5 as above, but produce a range of .pdb files that have different ratios of particle A:B for instance (1:999, 100:900, 500:500, 900:100, 999:1 etc.).
2. You'll need a new .top file for each of the different .pdb files specifying the number of each molecule but keeping everything else identical
3. Create the same .mdp files as before but replace with `prod_run.mdp` with `lambda.mdp`
4. The `lambda.mdp` file contains extra info regarding the how strongly the Van der Waals forces are between the two particles (coupled to $\lambda$). The example file contains a list of 14 different strengths ranging from 0 (no force) to 1 (usual force).
5. Repeat step 7 but using the `HPC Scripts/lambda.pbs` script for each ratio of A:B you've set up. This script should actually perform all the equilbriation *and* production runs at each lambda but I only have it set up to perform production runs at each lambda. This won't have much effect for small:large ratios but will probably cause a jump in conditions otherwise, so this needs to be fixed.
6.  After looping over all the specified lambdas it uses the BAR method to find the free energy of solvation of the mixture at that ratio from all the different lambda_i data. This outputs three additional files:
     * `bar.xvg`, containing the values of $\Delta G (KT)$ between each value of lambda
     * `barint.xvg` which is simply the cumulative version of the previous data
     * `histo.xvg` contains histograms of the  free energy distributions at adjacent lambda states.
7. Summing $\Delta G (KT)$ across all the $\lambda$ points will give you the free energy of solvation for a particular ratio, which you can then use to plot the free energy of solvation as a function of particle ratios.

