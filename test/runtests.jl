using   GeneticVariation, TableReader,
        CSV, GLM,
        Distributions, DataFrames,
        DataStructures, Distributions,
        MixedModels, StatsModels, BlockArrays,
        LinearAlgebra

include("../src/bcf_utils.jl")
include("../src/cmd_parse.jl")
include("../src/bcf_gwas.jl")
include("../src/bcf_hwe.jl")
include("../src/bcf_ldsc.jl")
include("../src/bcf_run.jl")

incmd = "--phenid analysis_id
--gwas pwav=age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10+rr_d
--bcf test/test22_aricaa.bcf.gz
--phepath test/ecg_whi_aric_sol_all.txt
--outfile test/test_aa22.txt
--famid analysis_fid
--ldsc"

pinputs = parse_cmd(incmd)
run_ldsc = haskey(pinputs, "ldsc")
run_gwas = haskey(pinputs, "gwas")
islmm = haskey(pinputs, "famid") & run_gwas
run_hwe = !haskey(pinputs, "nohwe")
maf = get(pinputs, "maf", 0.)
ldsc_maf = get(pinputs, "ldsc_maf", 0.01)
gwas_maf = get(pinputs, "gwas_maf", 0.)
sep = get(pinputs, "sep", '\t')[1]

#get list of BCF files
bcflist = get_bcflist(pinputs["bcf"])

#open output file
out = open(pinputs["outfile"], "w")
print_fun(v::GWAS_variant) = print_bcf(
  out, v; sep = sep, gwas = run_gwas, hwe = run_hwe, ldsc = run_ldsc, lmm = islmm
)
print_header(out; sep = '\t', gwas = run_gwas, hwe = run_hwe, ldsc = run_ldsc)


bcf = bcflist[1]
println("Now processing $bcf")

#Set up
vcf_phen, vcfind, Xs, y, form = prep_gwas(pinputs, bcf)
reader, key, vec, varnow = prep_run(bcf, vcfind)
var_deque = Deque{GWAS_variant}()
var_deque_lowmaf = Deque{GWAS_variant}()

#crappy temporary solution to find chromosome # from contig field
#needs updating for BCF files with multiple CHR
h = header(reader)
chloc = filter(x->metainfotag(h.metainfo[x])=="contig", eachindex(h.metainfo))[1]
charr = Char.(h.metainfo[chloc].data[h.metainfo[chloc].dictval[1]])
chrnow = String(charr)

#for testing
i = haskey(pinputs, "test") ? 20000 : -1

vcfnow = BCF.Record()
#Loop over BCF contents

read!(reader, vcfnow)
load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
# varnow.caf, varnow.qual
#skip if below global MAF threshold
# maf < varnow.caf < (1-maf) || continue

#overwrites incorrect chromosome number
varnow.chrom = chrnow

# Profile.clear
# Juno.@profiler process_var_lmm!(varnow, Xs, y, form)

#Runs GWAS
if run_gwas && ( gwas_maf < varnow.caf < (1-gwas_maf) )
    if islmm
        @time process_var_lmm!(varnow, Xs, y, form)
    else
        process_var_glm!(varnow, Xs, y)
    end
end

#Calculates HWE
run_hwe && process_var_hwe!(varnow, vcfind)

#Calculates LDSC
if run_ldsc
  process_var_ldsc!(varnow, var_deque, vardeque_lowmaf, print_fun)
else
    #Prints to file if no LDSC
    print_fun(varnow)
end
