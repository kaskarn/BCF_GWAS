function parse_cmd(incmd::String)
  pinputs = Dict()
  flags = ["hwe" "ldsc" "test"]
  for f in flags
    if occursin("--$(f)", incmd)
      replace(incmd, " --($f)" => "")
      get!(pinputs, f, true)
    else
      get!(pinputs, f, false)
    end
  end
  [get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true) for i in split(incmd, "--")[2:end]]

  pinputs
end
parse_cmd(incmd::Array) = parse_cmd(join(incmd, " "))
