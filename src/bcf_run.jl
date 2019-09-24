#Main function, calls upon all pieces
function process_bcf(incmd)
  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")
  display(pinputs)

  islmm = haskey(pinputs, "famid")
  run_ldsc = haskey(pinputs, "ldsc")
  bcflist = get_bcflist(pinputs["bcf"])

  out = open(pinputs["outfile"], "w")
  print_header(out; sep = '\t', ldsc = run_ldsc)

  println(bcflist)
  for bcf in bcflist
    println("Now processing $bcf")

    vcfind, Xs, y, form = prep_gwas(pinputs)
    reader, key, vec, varnow = prep_run(bcf, vcfind)
    var_deque = Deque{GWAS_variant}()
    # var_deque_lowmaf = Deque{GWAS_variant}()

    i = haskey(pinputs, "test") ? 20000 : -1

    #crappy temporary solution to find chromosome # from contig field
    h = header(reader)
    chloc = filter(x->metainfotag(h.metainfo[x])=="contig", eachindex(h.metainfo))[1]
    charr = Char.(h.metainfo[chloc].data[h.metainfo[chloc].dictval[1]])
    # chrnow = sum([10^(i-1)*charr[i] for i in reverse(eachindex(charr))])
    chrnow = String(charr)

    for vcfnow in reader
        load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)
        varnow.chrom = chrnow #overwrites incorrect chromosome number
        if 0.001 < varnow.caf < 0.999
            if islmm
                process_var_lmm!(Xs, y, form)
            else
                process_var_glm!(varnow, Xmat, y)
            end
        end
        process_var_hwe!(varnow, vcfind)
        if run_ldsc && (0.05 < varnow.caf < 0.95)
          process_var_ldsc!(varnow, var_deque, out) ## WIP
        else
          print_bcf(out, varnow, sep = '\t', ldsc = run_ldsc)
        end
        i % 1000 == 0 && println("$(-i); deque: $(length(var_deque))") #for testing
        i -= 1
        i == 0 && break
    end
    run_ldsc && ldsc_clear_deque!(out, var_deque)
  end
  close(out)
  return("Run completed!")
end
