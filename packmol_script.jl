import Pkg; Pkg.activate("Packmol", shared=true) # Activate Packmol environment
using Packmol
input_file = raw"\\rds.imperial.ac.uk\rds\user\mdo21\home\Lambda\0201\600-401\CG_WAT_TEA.inp"
run_packmol(input_file)