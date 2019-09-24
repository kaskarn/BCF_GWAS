#Sets up GWAS objects
prep_gwas(d::Dict) = prep_gwas(
    d["phepath"], d["bcf"], d["gwas"], d["phenid"],
    get(d, "famid", d["phenid"])
)

prep_gwas(d::Dict, bcf::String) = prep_gwas(
    d["phepath"], bcf, d["gwas"], d["phenid"],
    get(d, "famid", d["phenid"])
)

function prep_gwas(
    phenpath, bcf, gwas, phen_id,
    fam_id = phen_id;
    model = fam_id == phen_id ? LinearModel : LinearMixedModel
)
  reader = BCF.Reader(open(bcf, "r"))
  pheno = CSV.File(String(phenpath); delim = '\t', missingstring = "NA") |> DataFrame

  vcfids = DataFrame(id_rdy = reader.header.sampleID, ind = 1:length(reader.header.sampleID))
  close(reader)

  rename!(pheno, Symbol(phen_id) => :id_rdy)
  vcf_phen = join(pheno, vcfids, on = :id_rdy, kind = :inner)
  sort!(vcf_phen , order(:ind))

  lhs_s, rhs_s = split(gwas, '=')
  [ dropmissing!(vcf_phen, Symbol(x), disallowmissing=true) for x in split(rhs_s, "+") ]
  dropmissing!(vcf_phen, Symbol(lhs_s), disallowmissing=true)
  dropmissing!(vcf_phen, :ind, disallowmissing=true)
  if model == LinearMixedModel
      rhs_s = rhs_s*"+(1|"*fam_id*")"
      categorical!(vcf_phen, Symbol(fam_id))
  end
  formula = eval(Meta.parse("@formula $lhs_s ~ G + $rhs_s"))
  vcf_phen.G = randn(size(vcf_phen,1))

  #Below adapted from MixedModels source code
  #originally by Douglas Bates

  form = apply_schema(formula, schema(formula, vcf_phen), model)
  y, Xs = StatsModels.modelcols(form, vcf_phen)
  y = reshape(float(y), (:, 1)) # y as a floating-point matrix

  vcfind = vcf_phen.ind
  return vcfind, Xs, y, form
end

function process_var_lmm!(v::GWAS_variant, Xs, y, form, wts = [])
    Xs[1][:,2] = v.ds
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
    lbd = foldl(vcat, lowerbd(c) for c in reterms)
    θ = foldl(vcat, MixedModels.getθ(c) for c in reterms)
    optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
    fill!(optsum.xtol_abs, 1.0e-10)
    mnow = fit!(LinearMixedModel(form, reterms, feterms, sqrt.(convert(Vector{T}, wts)), A, L, optsum))
    v.b = fixef(mnow)[2]
    v.se = sqrt(vcov(mnow)[2,2])
    v.p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(v.b/v.se))
    v.b, v.se, v.p
end

#Run GLM using GWAS_variant
function process_var_glm!(v::GWAS_variant,
    Xmat::Matrix, y)

    Xmat[:,2] = v.ds

    mnow = lm(Xmat[v.ind,:],y[v.ind])
    v.b = coef(mnow)[2]
    v.se = stderror(mnow)[2]
    v.p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(v.b/v.se))
    v.b, v.se, v.p
end

function process_var_glm!(v::GWAS_variant,
    Xmat::Tuple{Array{Float64,2},ReMat{Float64,1}}, y)

    process_var_glm!(v, Xmat[1], y)
end

process_var_glm!(v::GWAS_variant, X::MixedModels.FeMat) = process_var_glm!(v,X[1].x,X[2].x)

# function process_var_glm!(Xmat, y, dsvec, gind)
#   Xmat[:,2] = dsvec
#
#   mnow = lm(Xmat[gind,:],y[gind])
#   b = coef(mnow)[2]
#   se = stderror(mnow)[2]
#   p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(b/se))
#   return b, se, p
# end





# function prep_lmm(Xs, y, form, wts = [])
#     T = eltype(y)
#
#     reterms = ReMat{T}[]
#     feterms = MixedModels.FeMat{T}[]
#
#     for (i,x) in enumerate(Xs)
#       if isa(x, ReMat{T})
#           push!(reterms, x)
#       else
#           cnames = coefnames(form.rhs[i])
#           push!(feterms, MixedModels.FeMat(x, isa(cnames, String) ? [cnames] : collect(cnames)))
#       end
#     end
#
#     push!(feterms, MixedModels.FeMat(y, [""]))
#     sort!(reterms, by=MixedModels.nranef, rev=true) #shoudn't change a thing
#
#     terms = vcat(reterms, feterms)
#     k = length(terms)
#     sz = append!(size.(reterms, 2), rank.(feterms))
#     A = BlockArrays._BlockArray(AbstractMatrix{T}, sz, sz)
#     L = BlockArrays._BlockArray(AbstractMatrix{T}, sz, sz)
#
#     for j in 1:k
#         for i in j:k
#             Lij = L[Block(i,j)] = MixedModels.densify(terms[i]'terms[j])
#             A[Block(i,j)] = deepcopy(isa(Lij, BlockedSparse) ? Lij.cscmat : Lij)
#         end
#     end
#     for i in 2:length(reterms) # check for fill-in due to non-nested grouping factors
#         ci = reterms[i]
#         for j in 1:(i - 1)
#             cj = reterms[j]
#             if !isnested(cj, ci)
#                 for l in i:k
#                     L[Block(l, i)] = Matrix(L[Block(l, i)])
#                 end
#                 break
#             end
#         end
#     end
#     lbd = foldl(vcat, lowerbd(c) for c in reterms)
#     θ = foldl(vcat, MixedModels.getθ(c) for c in reterms)
#     optsum = OptSummary(θ, lbd, :LN_BOBYQA, ftol_rel = T(1.0e-12), ftol_abs = T(1.0e-8))
#     fill!(optsum.xtol_abs, 1.0e-10)
#
#     return LinearMixedModel(form, reterms, feterms, sqrt.(convert(Vector{T}, wts)), A, L, optsum)
# end
#
# function process_var_lmm!(reterms, feterms, A, L, optsum)
#
#
#     mnow = fit!(LinearMixedModel(form, reterms, feterms, sqrt.(convert(Vector{T}, wts)), A, L, optsum))
#     v.b = fixef(mnow)[2]
#     v.se = sqrt(vcov(mnow)[2,2])
#     v.p = GLM.ccdf(GLM.FDist(1,GLM.dof_residual(mnow)),abs2(v.b/v.se))
#     v.b, v.se, v.p
# end
