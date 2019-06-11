#Update deque ldsc values with new variant
function ldsc_update_r2!(deque, varnow, pos_now)
    tot = 0.0
    @inbounds for var in deque
        newind = varnow.ind .& var.ind
        n = sum(newind)
        r2 = cor(varnow.ds[newind], var.ds[newind])
        r2_adj = (r2 - (1-r2)/(n-2)) #N-adjusted formula from LDSC paper
    	var.sumr2 += r2_adj
    	tot += r2_adj
    end
    varnow.sumr2 += tot

    return tot
end

#Clear deque variants outside window
function ldsc_update_deque_window!(out, vardeque, pos_now; win=1000000)
    while front(vardeque).pos + win < pos_now
        printbcf(out, popfirst!(vardeque); ldsc = true)
    end
    return 0
end
#Fully clear deque of variants
function ldsc_clear_deque!(out, vardeque)
    @inbounds for i in 1:length(vardeque)
    	printbcf(out, popfirst!(vardeque); ldsc = true)
    end
    return 0
end

#LDSC workhorse function
function process_var_ldsc!(varnow, vardeque, out, win=1000000)
    r2tot = 0.0
    if !isempty(vardeque)
        #If new chromosome, clear deque
        varnow.chrom != front(vardeque).chr && ldsc_clear_deque!(out, vardeque)

        #Remove variants outside window
        ldsc_update_deque_window!(out, vardeque, varnow.pos)

        #Update variant R2
        ldsc_update_r2!(vardeque, varnow.ds, varnow.pos)
    end
    push!(vardeque, varnow)

    return 0
end
