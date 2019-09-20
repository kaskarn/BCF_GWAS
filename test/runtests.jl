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
--famid analysis_fid"
# --ldsc"

pinputs = parse_cmd(incmd)
run_ldsc = haskey(pinputs, "ldsc")
bcflist = get_bcflist(pinputs["bcf"])

out = open(pinputs["outfile"], "w")
print_header(out; sep = '\t', ldsc = run_ldsc)

bcf = bcflist[1]
println("Now processing $bcf")

d = pinputs
# phenpath, bcf, gwas, phen_id = d["phepath"], d["bcf"], d["gwas"], d["phenid"]

vcfind, Xs, y, form = prep_gwas(pinputs)
reader, key, vec, varnow = prep_run(bcf, vcfind)
var_deque = Deque{GWAS_variant}()
# var_deque_lowmaf = Deque{GWAS_variant}()

i = haskey(pinputs, "test") ? 20000 : -1

#crappy temporary solution to find chromosome # from contig field
h = header(reader)
chloc = filter(x->metainfotag(h.metainfo[x])=="contig", eachindex(h.metainfo))[1]
charr = Int.(h.metainfo[chloc].data[h.metainfo[chloc].dictval[1]] .- 0x30)
chrnow = sum([10^(i-1)*charr[i] for i in reverse(eachindex(charr))])
vcfnow = read(reader)

load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
varnow.chrom = chrnow #overwrites incorrect chromosome number

process_var_glm!(varnow, Xs, y)
process_var_hwe!(varnow, vcfind)
if run_ldsc && (0.05 < varnow.caf < 0.95)
  process_var_ldsc!(varnow, var_deque, out) ## WIP
else
  print_bcf(out, varnow, sep = '\t', ldsc = run_ldsc)
end



@time yresp = GLM.LmResp(y[varnow.ind]);
@time cp = GLM.cholpred(Xs[1], false);
@time mnow = LinearModel(yresp, cp)
@time fit!(mnow);

@time lm(Xs[1][varnow.ind,:], y[varnow.ind]);

@benchmark fit!(LinearMixedModel(formula, pheno))
@benchmark refit!(mnow)
@benchmark process_var_lmm!(varnow, Xs, y, form)

process_var_glm!(varnow, Xs, y)



process_var_hwe!(varnow, vcfind)
if run_ldsc && (0.05 < varnow.caf < 0.95)
  process_var_ldsc!(varnow, var_deque, out) ## WIP
else
  print_bcf(out, varnow, sep = '\t', ldsc = run_ldsc)
end
i % 1000 == 0 && println("$(-i); deque: $(length(var_deque))") #for testing
i -= 1
i == 0 && break
