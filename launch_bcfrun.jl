using Pkg

#Checks all the required packages are installed
function checkInstall()
    mypkg = Pkg.installed()
    pkgs = ["GeneticVariation" "TableReader" "CSV" "GLM" "Distributions" "DataFrames" "DataStructures"]
    [ Pkg.add(p) for p in filter(x->!haskey(mypkg, x), pkgs) ]
end
checkInstall()

using   GeneticVariation, TableReader,
        CSV, GLM,
        Distributions, DataFrames,
        DataStructures, Distributions

include(Base.source_dir()*"/cmd_parse.jl")
include(Base.source_dir()*"/bcf_utils.jl")
include(Base.source_dir()*"/bcf_gwas.jl")
include(Base.source_dir()*"/bcf_hwe.jl")
include(Base.source_dir()*"/launch_bcfrun.jl")

println("\nCOMMAND:")
println(join(ARGS, " "))
@time process_bcf(ARGS)
