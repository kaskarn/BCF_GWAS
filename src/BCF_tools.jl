module BCF_tools

using   GeneticVariation, TableReader,
        CSV, GLM,
        Distributions, DataFrames,
        DataStructures, Distributions,
        MixedModels, StatsModels, BlockArrays

include("bcf_utils.jl")
include("cmd_parse.jl")
include("bcf_gwas.jl")
include("bcf_hwe.jl")
include("bcf_ldsc.jl")
include("bcf_run.jl")


end # module
