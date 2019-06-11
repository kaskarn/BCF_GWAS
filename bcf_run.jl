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

    Xmat, y, vcfind = prep_gwas(pinputs)
    reader, key, vec, varnow = prep_run(bcf, vcfind)
    var_deque = Deque{GWAS_variant}()
    var_deque_lowmaf = Deque{GWAS_variant}()

    out = open(pinputs["outfile"], "w+")
    print_header(out; sep = '\t', ldsc = run_ldsc)

    i = pinputs["test"] ? 10000 : -1
    @time for vcfnow in reader
      load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
      (0.01 < varnow.caf < 0.99) && process_var_glm!(varnow, Xmat, y)
      process_var_hwe!(varnow, vcfind)
      if run_ldsc
          process_var_ldsc!(varnow, var_deque, var_deque_lowmaf, out) ## WIP
      else
          print_bcf(out, varnow, sep = '\t', ldsc = run_ldsc)
      end
      i -= 1
      i == 0 && break
    end
    run_ldsc && ldsc_clear_deque!(out, var_deque, var_deque_lowmaf)
  end
  close(out)
  return("Run completed!")
end
