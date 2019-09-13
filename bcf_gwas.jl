#Sets up GWAS objects
prep_gwas(d::Dict) = prep_gwas(d["phepath"], d["bcf"], d["gwas"], d["phenid"])

# mutable struct FeMat{T,S<:AbstractMatrix}
#     x::S
#     wtx::S
#     piv::Vector{Int}
#     rank::Int
#     cnames::Vector{String}
# end

# mutable struct ReMat{T,S} <: AbstractMatrix{T}
#     trm::CategoricalTerm
#     refs::Vector{Int32}
#     cnames::Vector{String}
#     z::Matrix{T}
#     wtz::Matrix{T}
#     λ::LowerTriangular{T,Matrix{T}}
#     inds::Vector{Int}
#     adjA::SparseMatrixCSC{T,Int32}
#     scratch::Matrix{T}
# end

function prep_gwas(
    phenpath, bcf, gwas, phen_id, fam_id = phen_id;
    model = fam_id == phen_id ? LinearModel : LinearMixedModel
)

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

  # gind = trues(length(vcfind))
  vcf_phen[:G] = randn(size(vcf_phen,1))
  #
  # tmdf = fit(LinearModel, formula, vcf_phen).mf.df
  # Xmat = convert(Matrix{Float64}, tmdf)
  # Xmat[:,1] .= 1.
  # y = convert(Vector{Float64}, vcf_phen[lhs] )

  #Below adapted from MixedModels source code
  #originally by Douglas Bates

  form = apply_schema(formula, schema(formula, vcf_phen), model)
  y, Xs = StatsModels.modelcols(form, vcf_phen)
  y = reshape(float(y), (:, 1)) # y as a floating-point matrix
  T = eltype(y)

  reterms = ReMat{T}[]
  feterms = MixedModels.FeMat{T}[]

  for (i,x) in enumerate(Xs)
    if isa(x, ReMat{T})
        push!(reterms, x)
    else
        cnames = coefnames(form.rhs[i])
        push!(feterms, MixedModels.FeMat(x, isa(cnames, String) ? [cnames] : collect(cnames)))
    end
  end

  push!(feterms, MixedModels.FeMat(y, [""]))
  sort!(reterms, by=MixedModels.nranef, rev=true) #shoudn't change a thing

  terms = vcat(reterms, feterms)
  k = length(terms)
  sz = append!(size.(reterms, 2), rank.(feterms))
  A = BlockArrays._BlockArray(AbstractMatrix{T}, sz, sz)
  L = BlockArrays._BlockArray(AbstractMatrix{T}, sz, sz)

  optsum = OptSummary(1.0, 0.0, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))

  #no ReMat to set if Linear Model
  (model == LinearModel) && return feterms, reterms, A, L, optsum

  #Else, finish setup
  for j in 1:k
      for i in j:k
          Lij = L[Block(i,j)] = MixedModels.densify(terms[i]'terms[j])
          A[Block(i,j)] = deepcopy(isa(Lij, BlockedSparse) ? Lij.cscmat : Lij)
      end
  end
  for i in 2:length(reterms) # check for fill-in due to non-nested grouping factors
      ci = reterms[i]
      for j in 1:(i - 1)
          cj = reterms[j]
          if !isnested(cj, ci)
              for l in i:k
                  L[Block(l, i)] = Matrix(L[Block(l, i)])
              end
              break
          end
      end
  end

  return feterms, reterms, A, L, optsum

end

#Run GLM using GWAS_variant
function process_var_glm!(v::GWAS_variant, Xmat, y)
  Xmat[:,2] = v.ds

  mnow = lm(Xmat[v.ind,:],y[v.ind])
  v.b = coef(mnow)[2]
  v.se = stderror(mnow)[2]
  v.p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(v.b/v.se))
end

# function process_var_glm!(Xmat, y, dsvec, gind)
#   Xmat[:,2] = dsvec
#
#   mnow = lm(Xmat[gind,:],y[gind])
#   b = coef(mnow)[2]
#   se = stderror(mnow)[2]
#   p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(b/se))
#   return b, se, p
# end

process_var_glm!(v::GWAS_variant, X::MixedModels.FeMat) = process_var_glm!(v,X[1].x,X[2].x)
