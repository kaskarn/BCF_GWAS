function prep_gwas(phenpath, bcf, gwas, phen_id)
  reader = BCF.Reader(open(bcf, "r"))
  pheno = CSV.File(String(phenpath); delim = '\t') |> DataFrame

  vcfids = DataFrame(id_rdy = reader.header.sampleID, ind = 1:length(reader.header.sampleID))
  close(reader)

  rename!(pheno, Symbol(phen_id) => :id_rdy)
  vcf_phen = join(pheno, vcfids, on = :id_rdy, kind = :inner)
  sort!(vcf_phen , order(:ind))
  vcfind = vcf_phen[:ind]

  lhs_s, rhs_s = split(gwas, '=')

  rhs = Symbol.(split("G+"*String(rhs_s), "+"))
  lhs = Symbol(lhs_s)
  formula = @eval(@formula($lhs ~ (+)(1, $(rhs...))))

  gind = trues(length(vcfind))
  vcf_phen[:G] = randn(size(vcf_phen,1))
  tmdf = fit(LinearModel, formula, vcf_phen).mf.df
  Xmat = convert(Matrix{Float64}, tmdf)
  Xmat[:,1] .= 1.
  y = convert(Vector{Float64}, vcf_phen[lhs] )

  return Xmat, y, vcfind
end

function prep_gwas(d::Dict)
  prep_gwas(d["phepath"], d["bcf"], d["gwas"], d["phenid"])
end

function process_var_glm!(Xmat, y, dsvec, gind)
  Xmat[:,2] = dsvec

  mnow = lm(Xmat[gind,:],y[gind])
  b = coef(mnow)[2]
  se = stderror(mnow)[2]
  p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(b/se))
  return b, se, p
end

function process_var_glm!(v::GWAS_variant, Xmat, y)
  Xmat[:,2] = v.ds

  mnow = lm(Xmat[v.ind,:],y[v.ind])
  v.b = coef(mnow)[2]
  v.se = stderror(mnow)[2]
  v.p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(v.b/v.se))
end
