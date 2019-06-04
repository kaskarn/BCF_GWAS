mutable struct GWAS_variant
    id::String
    id2::String
    pg::Matrix{Float64}
    ds::Vector{Float64}
    ind::BitArray{1}
    chrom::Int64
    pos::Int64
    ref::String
    alt::String
    qual::Float64
    n::Int64
    ac::Float64
    caf::Float64
    b::Float64
    se::Float64
    p::Float64
    n_aa::Float64
    n_Aa::Float64
    n_AA::Float64
    hwe::Float64
    sumr2::Float64
end

function GWAS_variant(n)
  varnow = GWAS_variant(
    "", "",
    zeros(Float64, 3, n), zeros(Float64, n), falses(n),
    -1, -1, "", "", NaN,
    -1, NaN, NaN, #alleles N
    NaN, NaN, NaN, #gwas
    NaN, NaN, NaN, NaN, #hwe
    NaN
  )
end

function bcftype2(b::UInt8)
  b &= 0x0f #take right 4 bits
  t::DataType = b == 0x01 ? Int8 :
      b == 0x02 ? Int16 :
      b == 0x03 ? Int32 :
      b == 0x05 ? Float32 :
      b == 0x07 ? Char : Nothing
  t
end

function genotype2(record::GeneticVariation.BCF.Record, key::Integer)
    GeneticVariation.BCF.checkfilled(record)
    N = GeneticVariation.BCF.n_sample(record)
    offset::Int = record.sharedlen
    for j in 1:GeneticVariation.BCF.n_format(record)
        k, offset = GeneticVariation.BCF.loadvec(record.data, offset)
        @assert length(k) == 1
        head, offset = GeneticVariation.BCF.loadvechead(record.data, offset)
        typ = bcftype2(record.data[offset])
        s = head[2]*sizeof(typ)
        if k[1] == key
            vals = Matrix{typ}(undef, head[2], N)
            ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), pointer(vals), pointer(record.data, offset+1), N*s)
            return vals
        end
        offset += N*s
    end
    throw(KeyError(key))
end


function genotype2!(record::GeneticVariation.BCF.Record, key::Integer, vals::Array)
    GeneticVariation.BCF.checkfilled(record)
    N = GeneticVariation.BCF.n_sample(record)
    offset::Int = record.sharedlen
    for j in 1:GeneticVariation.BCF.n_format(record)
        k, offset = GeneticVariation.BCF.loadvec(record.data, offset)
        @assert length(k) == 1
        head, offset = GeneticVariation.BCF.loadvechead(record.data, offset)
        typ = bcftype2(record.data[offset])
        s = head[2]*sizeof(typ)
        if k[1] == key
            sizeof(vals) != N*s && throw(ErrorException("Target array incorrect size. Is $(sizeof(vals)) bytes, should be $(N*s)"))
            ccall(:memcpy, Ptr{Cvoid}, (Ptr{Cvoid}, Ptr{Cvoid}, Csize_t), pointer(vals), pointer(record.data, offset+1), N*s)
            return offset
        end
        offset += N*s
    end
    throw(KeyError(key))
end

function bcf_findkey(reader, s::String)
  h = header(reader)
  format_ind = filter(x->metainfotag(h.metainfo[x])=="FORMAT", eachindex(h.metainfo))
  gp_ind = filter(y->h.metainfo[y]["ID"] == s, format_ind)[1]
  gp_key_s = h.metainfo[gp_ind]["IDX"]
  parse(Int, gp_key_s)
end

function prep_run(bcf::AbstractString, vcfind::Vector{Int})
  reader = BCF.Reader(open(bcf, "r"))
  v = read(reader)
  key = bcf_findkey(reader, "GP")
  vec = genotype2(v, key)
  close(reader)
  return BCF.Reader(open(bcf, "r")), key, vec, GWAS_variant(length(vcfind))
end

function load_bcf_variant!(v::GWAS_variant, vec, vcfnow, key, vcfind)
  genotype2!(vcfnow, key, vec)
  v.n, v.ac = 0, 0.0
  @inbounds for (i, k) in enumerate(vcfind)
    v.pg[:,i] = vec[:,k]
    ds = v.pg[2,i] + 2*v.pg[3,i]
    v.ds[i] = ds
    v.ind[i] = !isnan(ds) #missing BCF2 standard
    if v.ind[i]
      v.n += 1
      v.ac += ds
    end
  end
  v.caf = v.ac/v.n

  #init variant info
  v.id = BCF.id(vcfnow)
  v.id2 = BCF.info(vcfnow)[1][2]
  v.pos = BCF.pos(vcfnow)
  v.chrom = BCF.chrom(vcfnow)
  v.ref = BCF.ref(vcfnow)
  v.alt = BCF.alt(vcfnow)[1]
  v.qual = BCF.qual(vcfnow)
  
  #init GLM results
  v.b = NaN
  v.se = NaN
  v.p = NaN

  #init LDS results
  v.sumr2 = 0.

  #init HWE results
  v.hwe = NaN

  #return 0
  0.
end

function print_header(io = stdout; sep = '\t', gwas = false, hwe = false, ldsc = false)
  join(io, ["VCF_ID", "ID2", "CHROM", "POS", "REF", "ALT", "ALT_AF", "ALT_AC", "N_INFORMATIVE", "QUAL"], sep)
  gwas && join(io, ["", "BETA", "SE", "PVALUE"], sep)
  hwe && join(io, ["", "N_aa", "N_Aa", "N_AA", "HWE_p"], sep)
  ldsc && join(io, ["", "LDSC"], sep)
  write(io, "\n")
end


function print_bcf(io, v::GWAS_variant; sep = '\t', gwas = false, hwe = false, ldsc = false)
  join(io, [v.id, v.id2, v.chrom, v.pos, v.ref, v.alt, v.caf, v.ac, v.n, v.qual], sep)
  gwas && join(io, ["", v.b, v.se, v.p], sep)
  hwe && join(io, ["", v.n_aa, v.n_Aa, v.n_AA, v.hwe], sep)
  ldsc && join(io, ["", v.sumr2], sep)
  write(io, '\n')
end
