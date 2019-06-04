function parse_cmd(incmd::String)
  pinputs = Dict()
  [ get!(pinputs, split(i)[1], size(split(i),1) > 1 ? split(i)[2] : true) for i in split(incmd, "--")[2:end] ]
  pinputs
end
parse_cmd(incmd::Array) = parse_cmd(join(incmd, " "))
