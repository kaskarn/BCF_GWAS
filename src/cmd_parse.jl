#Takes in arg string, parses into dictionary objects
function parse_cmd(incmd::AbstractString)
  pinputs = Dict()
  [ #split into dictionary elements
    get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true)
    for i in split(incmd, "--")[2:end]
  ]

  pinputs
end
parse_cmd(incmd::Array) = parse_cmd(join(incmd, " "))


#Returns list of bcf files, ordered by chromosome
#finds .bcf.gz files in specified directory
#pattern chr# used to order files -> (chr1, chr2, ..., chr22, ..., chrX)
function get_bcflist(path::AbstractString)
    bcflist = isdir(path) ? filter(x->!isnothing(match(r"bcf.gz$", x)), readdir(path)) : [path]
    if length(bcflist) > 1
      fint_ind = [ parse(Int, match(r"chr([0-9]+)", s).captures[1]) for s in bcflist ]
      f_int = bcflist[sortperm(fint_ind)]
      bcflist = [ path*"/"*i for i in [f_int..., setdiff(bcflist, f_int)...] ]
    elseif length(bcflist) == 0
      throw("$(path) has no BCF files (.bcf.gz) !\n")
    else
      isfile(bcflist[1]) || throw("$(bcflist[1]) does not exist\n")
    end
    bcflist
end
