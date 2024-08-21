import Pkg; Pkg.activate("Packmol", shared=true) # Activate Packmol environment
using Packmol
input_file = raw"file_path"
run_packmol(input_file)
