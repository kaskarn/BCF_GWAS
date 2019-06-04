using Pkg
function checkInstall()
    mypkg = Pkg.installed()
    pkgs = ["GeneticVariation" "TableReader" "CSV" "GLM" "Distributions" "DataFrames" "StaticArrays" "DataStructures"]
    [ Pkg.add(p) for p in filter(x->!haskey(mypkg, x), pkgs) ]
end

checkInstall()

using GeneticVariation, TableReader, CSV, GLM, Distributions, DataFrames, StaticArrays, DataStructures, Distributions

include(Base.source_dir()*"/cmd_parse.jl")
include(Base.source_dir()*"/bcf_utils.jl")
include(Base.source_dir()*"/bcf_gwas.jl")
include(Base.source_dir()*"/bcf_hwe.jl")

#var_deque = Deque{GWAS_variant}();

function process_bcf(incmd)
  println("\nCOMMAND:")
  println(incmd)

  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")
  display(pinputs)

  gwas, hwe, ldsc = (haskey(pinputs, k) for k in ["gwas" "hwe" "ldsc"])
  
  bcflist = isdir(pinputs["bcf"]) ? filter(x->match(r"bcf.gz$", x)!= nothing, readdir(pinputs["bcf"])) : [ pinputs["bcf"] ]
  if length(bcflist) > 1
    fint_ind = [ parse(Int, match(r"chr([0-9]+)", s).captures[1]) for s in bcflist ]
    f_int = bcflist[sortperm(fint_ind)]
    bcflist = [ pinputs["bcf"]*"/"*i for i in [f_int..., setdiff(bcflist, f_int)...] ]
  elseif length(bcflist) == 0
    throw("$(pinputs["bcf"]) is empty!\n")
  else
    isfile(bcflist[1]) || throw("$(bcflist[1]) does not exist\n")
  end
  
  println(bcflist)
  for bcf in bcflist
    pinputs["bcf"] = bcf
    println(pinputs["bcf"])
    Xmat, y, vcfind = prep_gwas(pinputs)
    reader, key, vec, varnow = prep_run(pinputs["bcf"], vcfind)

    out = open(pinputs["outfile"], "w+")
    print_header(out, sep = '\t', gwas = gwas, ldsc = ldsc, hwe = hwe)

    vcfnow = read(reader)

    i = pinputs["test"] ? 10000 : -1
    for vcfnow in reader
      load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
      gwas && (0.01 < varnow.caf < 0.99) && process_var_glm!(varnow, Xmat, y)
      hwe && process_var_hwe!(varnow, vcfind)
      # ldsc && process_var_ldsc!(varnow, var_deque)
      print_bcf(out, varnow, sep = '\t', gwas = gwas, hwe = hwe, ldsc = ldsc)
      i -= 1
      i == 0 && break
    end
  end

  return("Run completed!")
end

#GO!
@time process_bcf(ARGS)
