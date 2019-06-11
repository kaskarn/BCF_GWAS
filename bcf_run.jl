#Main function, calls upon all pieces
function process_bcf(incmd)
  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")
  display(pinputs)

  run_ldsc = haskey(pinputs, "ldsc")
  bcflist = get_bcflist(pinputs["bcf"])

  println(bcflist)
  for bcf in bcflist
    println("Now processing $bcf")

    Xmat, y, vcfind = prep_gwas(pinputs["phepath"], bcf, d["gwas"], d["phenid"])
    reader, key, vec, varnow = prep_run(bcf, vcfind)
    var_deque = Deque{GWAS_variant}()

    out = open(pinputs["outfile"], "w+")
    print_header(out; sep = '\t', ldsc = run_ldsc)

    vcfnow = read(reader)

    i = pinputs["test"] ? 10000 : -1
    for vcfnow in reader
      load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
      (0.01 < varnow.caf < 0.99) && process_var_glm!(varnow, Xmat, y)
      process_var_hwe!(varnow, vcfind)
      if run_ldsc && (0.05 < varnow.caf < 0.95)
          process_var_ldsc!(varnow, var_deque, out) ## WIP
      else
          print_bcf(out, varnow, sep = '\t')
      end

      i -= 1
      i == 0 && break
    end
  end
  return("Run completed!")
end
