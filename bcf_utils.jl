mutable struct GWAS_variant
    id::String
    pg::Matrix{Float64}
    ds::Vector{Float64}
    ind::BitArray{1}
    chrom::Int64
    pos::Int64
    ref::String
    alt::String
    n::Int64
    ac::Float64
    caf::Float64
    b::Float64
    se::Float64
    p::Float64
    hwe::Float64
    sumr2::Float64
end

function GWAS_variant(n)
  varnow = GWAS_variant(
    "",
    zeros(Float64, 3, n), zeros(Float64, n), falses(n),
    0, 0, "", "",
    0, 0, 0.0,
    0.0, 0.0, 0.0,
    0.0, 0.0
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

function prep_run(bcf, vcfind, key = 3)
  reader = BCF.Reader(open(bcf, "r"))
  vcfnow = read(reader)
  vec = genotype2(vcfnow, key)

  close(reader)
  return BCF.Reader(open(bcf, "r")), vec, GWAS_variant(length(vcfind))
end

function load_bcf_variant!(v, vec, vcfnow, key, vcfind)
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
  v.pos = BCF.pos(vcfnow)
  v.chrom = BCF.chrom(vcfnow)
  v.ref = BCF.ref(vcfnow)
  v.alt = BCF.alt(vcfnow)[1]

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
  join(io, ["VCF_ID", "CHROM", "POS", "REF", "ALT", "ALT_AF", "ALT_AC", "N_INFORMATIVE"], sep)
  gwas && join(io, ["", "BETA", "SE", "PVALUE"], sep)
  hwe && write(io, "HWE_p")
  ldsc && write(io, "LDSC")
  write(io, "\n")
end


function print(io, v::GWAS_variant; sep = '\t', gwas = false, hwe = false, ldsc = false)
  join(io, [v.id, v.chrom, v.pos, v.ref, v.alt, v.caf, v.ac, v.n], sep)
  gwas && join(io, ["", v.b, v.se, v.p], sep)
  hwe && join(io, ["", v.hwe], sep)
  ldsc && join(io, ["", v.sumr2], sep)
  write(io, '\n')
end
