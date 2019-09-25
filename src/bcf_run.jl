#Main function, calls upon all pieces
function process_bcf(incmd)
  pinputs = parse_cmd(incmd)
  println("\nParsed Arguments:")

  display(pinputs)
  println("\n")

  #Grab arguments from dictionary
  islmm = haskey(pinputs, "famid")
  run_ldsc = haskey(pinputs, "ldsc")
  run_gwas = haskey(pinputs, "gwas")
  run_hwe = !haskey(pinputs, "nohwe")
  maf = get(pinputs, "maf", 0.)
  ldsc_maf = get(pinputs, "ldsc_maf", 0.01)
  gwas_maf = get(pinputs, "gwas_maf", 0.)
  sep = get(pinputs, "sep", '\t')[1]

  #get list of BCF files
  bcflist = get_bcflist(pinputs["bcf"])

  #open output file
  out = open(pinputs["outfile"], "w")
  print_fun(v::GWAS_variant) = print_bcf(
    out, v; sep = sep, gwas = run_gwas, hwe = run_hwe, ldsc = run_ldsc
  )
  print_header(out; sep = '\t', gwas = run_gwas, hwe = run_hwe, ldsc = run_ldsc)

  #Loop over list of BCF files
  for bcf in bcflist
    println("Now processing $bcf")

    #Set up
    vcfind, Xs, y, form = prep_gwas(pinputs, bcf)
    reader, key, vec, varnow = prep_run(bcf, vcfind)
    var_deque = Deque{GWAS_variant}()
    var_deque_lowmaf = Deque{GWAS_variant}()

    #crappy temporary solution to find chromosome # from contig field
    #needs updating for BCF files with multiple CHR
    h = header(reader)
    chloc = filter(x->metainfotag(h.metainfo[x])=="contig", eachindex(h.metainfo))[1]
    charr = Char.(h.metainfo[chloc].data[h.metainfo[chloc].dictval[1]])
    chrnow = String(charr)

    #for testing
    i = haskey(pinputs, "test") ? 20000 : -1

    vcfnow = BCF.Record()
    #Loop over BCF contents
    while !eof(reader)
        read!(reader, vcfnow)
        load_bcf_variant!(varnow, vec, vcfnow, key, vcfind)

        #skip if below global MAF threshold
        maf < varnow.caf < (1-maf) && continue

        #overwrites incorrect chromosome number
        varnow.chrom = chrnow

        #Runs GWAS
        if run_gwas && ( gwas_maf < varnow.caf < (1-gwas_maf) )
            if islmm
                process_var_lmm!(varnow, Xs, y, form)
            else
                process_var_glm!(varnow, Xs, y)
            end
        end

        #Calculates HWE
        run_hwe && process_var_hwe!(varnow, vcfind)

        #Calculates LDSC
        if run_ldsc
          process_var_ldsc!(varnow, var_deque, vardeque_lowmaf, print_fun)
        else
            #Prints to file if no LDSC
            print_fun(varnow)
        end

        #for testing
        # i % 1000 == 0 && println("$(-i); deque: $(length(var_deque))") #for testing
        i -= 1
        i == 0 && break
    end
    run_ldsc && ldsc_clear_deque!(var_deque, print_fun)
  end
  close(out)
  return("Run completed!")
end
