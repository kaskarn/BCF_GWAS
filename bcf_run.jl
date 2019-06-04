using Pkg
function checkInstall()
    mypkg = Pkg.installed()
    pkgs = ["GeneticVariation" "TableReader" "CSV" "GLM" "Distributions" "DataFrames" "StaticArrays" "DataStructures"]
    [ Pkg.add(p) for p in filter(x->!haskey(mypkg, x), pkgs) ]
end

checkInstall()

using GeneticVariation, TableReader, CSV, GLM, Distributions, DataFrames, StaticArrays, DataStructures

include(Base.source_dir()*"/cmd_parse.jl")
include(Base.source_dir()*"/bcf_utils.jl")
include(Base.source_dir()*"/bcf_gwas.jl")

function process_bcf(incmd)
  key = 3
  println("\nCOMMAND:")
  println(incmd)

  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")
  display(pinputs)

  gwas, hwe, ldsc = (haskey(pinputs, k) for k in ["gwas" "ldsc" "hwe"])

  Xmat, y, vcfind = prep_gwas(pinputs)
  reader, vec, varnow = prep_run(pinputs["bcf"], vcfind, key)

  var_deque = Deque{GWAS_variant}();

  out = open(pinputs["outfile"], "w+")
  print_header(out, sep = '\t', gwas = gwas, ldsc = ldsc, hwe = hwe)

  vcfnow = read(reader)

  i = pinputs["test"] ? 1000 : -1
  for vcfnow in reader
    load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
    gwas && (0.01 < varnow.caf < 0.99) && process_var_glm!(varnow, Xmat, y)
    hwe && process_var_hwe!(varnow, vcfind)
    # ldsc && process_var_ldsc!(varnow, var_deque)
    print_bcf(out, varnow, sep = '\t', gwas = gwas, hwe = hwe, ldcs = ldsc)
    i -= 1
    i == 0 && return("Test run (1000 variants) completed!")
  end

  return("Run completed!")
end


incmd = "
--gwas twav=age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10
--bcf /proj/epi/CVDGeneNas/antoine/PAGE_BCF/ARIC_AA_gwas_frz3/ARIC_AA_gwas_frz3.chr1.bcf.gz
--phepath /proj/epi/CVDGeneNas/antoine/ecg_sixtraits_gwas_project/phenotypes/all/ecg_whi_aric_sol_all.txt
--phenid analysis_id
--outfile test.txt
--hwe
--test"



#GO!
@time process_bcf(incmd)
