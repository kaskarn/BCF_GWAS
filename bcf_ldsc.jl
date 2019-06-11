###############
#functions with 1 Deque
##############

function ldsc_update_r2!(deque, varnow, pos_now)
    tot = 0.0
    @inbounds for var in deque
        0.05 < var.caf < 0.95 && continue
        newind = varnow.ind .& var.ind
        n = sum(newind)
        r2 = cor(varnow.ds[newind], var.ds[newind])
        r2_adj = (r2 - (1-r2)/(n-2)) #N-adjusted formula from LDSC paper
    	var.sumr2 += r2_adj
    	tot += r2_adj
    end
    varnow.sumr2 = tot

    return tot
end

#Clear deque variants outside window
function ldsc_update_deque_window!(out, vardeque, pos_now; win=1000000)
    while front(vardeque).pos + win < pos_now
        print_bcf(out, popfirst!(vardeque); ldsc = true)
    end
    return 0
end
#Fully clear deque of variants
function ldsc_clear_deque!(out, vardeque)
    @inbounds for i in 1:length(vardeque)
    	print_bcf(out, popfirst!(vardeque); ldsc = true)
    end
    return 0
end

#LDSC workhorse function
function process_var_ldsc!(varnow, vardeque, out; win=1000000)
    r2tot = 0.0
    if !isempty(vardeque)
        #If new chromosome, clear deque
        varnow.chrom != front(vardeque).chrom && ldsc_clear_deque!(out, vardeque)

        #Remove variants outside window
        ldsc_update_deque_window!(out, vardeque, varnow.pos)

        #Update variant R2
        if 0.05 < varnow.caf < 0.95
            ldsc_update_r2!(vardeque, varnow, varnow.pos)
        end
    end
    push!(vardeque, copy(varnow))

    return 0
end

###############
#functions with 2 deques
##############

#Clear deque variants outside window
function ldsc_update_deque_window!(out, vardeque, vardeque_lowmaf, pos_now; win=1000000)
    while front(vardeque).pos + win < pos_now
        varnow = popfirst!(vardeque)
        while(first(vardeque_lowmaf).pos < varnow.pos)
            print_bcf(out, popfirst!(vardeque_lowmaf), ldsc = true)
        end
        print_bcf(out, popfirst!(vardeque); ldsc = true)
    end
    return 0
end
#Fully clear deque of variants
function ldsc_clear_deque!(out, vardeque, vardeque_lowmaf)
    for i in 1:length(vardeque)
        varnow = popfirst!(vardeque)
        while(first(vardeque_lowmaf).pos < varnow.pos)
            print_bcf(out, popfirst!(vardeque_lowmaf), ldsc = true)
        end
    	print_bcf(out, varnow; ldsc = true)
    end
    ldsc_clear_deque!(out, vardeque_lowmaf)
    return 0
end

function process_var_ldsc!(varnow, vardeque, vardeque_lowmaf, out; win=1000000)
    r2tot = 0.0
    if !isempty(vardeque)
        #If new chromosome, clear deque
        varnow.chrom != front(vardeque).chrom && ldsc_clear_deque!(out, vardeque)

        #Remove variants outside window
        ldsc_update_deque_window!(out, vardeque, vardeque_lowmaf, varnow.pos)

        #Update variant R2
        if 0.05 < varnow.caf < 0.95
            ldsc_update_r2!(vardeque, varnow, varnow.pos)
        end
    end
    if 0.05 < varnow.caf < 0.05
        push!(vardeque, copy(varnow))
    else
        push!(vardeque_lowmaf, copy(varnow))
    end

    return 0
end
