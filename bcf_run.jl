using GeneticVariation, BenchmarkTools, TableReader, CSV, GLM, Distributions, DataFrames, StaticArrays
using DataStructures


source("cmd_parse.jl")
source("bcf_utils.jl")
source("bcf_gwas.jl")

function process_bcf(incmd)
  key = 3
  println("\nCOMMAND:")
  println(incmd)

  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")
  @show pinputs

  gwas, hwe, ldsc = (haskey(pinputs, k) for k in ["gwas" "ldsc" "hwe"])

  gwas && (Xmat, y, vcfind = prep_gwas(pinputs))
  reader, vec, varnow = prep_run(pinputs["bcfpath"], vcfind, key)

  var_deque = Deque{GWAS_variant}();

  out = open(pinputs["outfile"], "w+")
  print_header(out, sep = '\t', gwas = gwas, ldsc = ldsc, hwe = hwe)

  vcfnow = read(reader)

  i = pinputs["test"] ? 1000 : -1
  for vcfnow in reader
    load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
    gwas && (0.01 < varnow.caf < 0.99) && process_var_glm!(varnow, Xmat, y)
    # hwe && process_var_hwe!(varnow)
    # ldsc && process_var_ldsc!(varnow, var_deque)
    print(out, varnow, sep = '\t', gwas = gwas, hwe = hwe, ldcs = ldsc)
    i -= 1
    i == 0 && return("Test run (1000 variants) completed!")
  end

  return("Run completed!")
end


incmd = "
--gwas twav=age+EV1+EV2+EV3+EV4+EV5+EV6+EV7+EV8+EV9+EV10
--study ARIC_AA
--phepath /proj/epi/CVDGeneNas/antoine/ecg_sixtraits_gwas_project/phenotypes/all/ecg_whi_aric_sol_all.txt
--phenid analysis_id
--outfile test.txt
--hwe"



#GO!
@time process_bcf(incmd)
