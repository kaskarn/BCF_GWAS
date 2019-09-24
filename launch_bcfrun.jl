using Pkg

#Checks all the required packages are installed
function checkInstall()
    mypkg = Pkg.installed()
    pkgs = ["GeneticVariation" "TableReader" "CSV" "GLM" "Distributions" "DataFrames" "DataStructures" "MixedModels" "StatsModels" "BlockArrays"]
    [ Pkg.add(p) for p in filter(x->!haskey(mypkg, x), pkgs) ]
end
checkInstall()

using   GeneticVariation, TableReader,
        CSV, GLM,
        Distributions, DataFrames,
        DataStructures, Distributions,
        MixedModels, StatsModels, BlockArrays,
        LinearAlgebra

include(Base.source_dir()*"/src/bcf_utils.jl")
include(Base.source_dir()*"/src/cmd_parse.jl")
include(Base.source_dir()*"/src/bcf_gwas.jl")
include(Base.source_dir()*"/src/bcf_hwe.jl")
include(Base.source_dir()*"/src/bcf_ldsc.jl")
include(Base.source_dir()*"/src/bcf_run.jl")

println("\nCOMMAND:")
println(join(ARGS, " "))
@time process_bcf(ARGS)
