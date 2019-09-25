###############
#functions with 1 Deque
##############

# function Base.popfirst!(s::Stack{T}) where (T <: Any)
#     pop!(s)
# end
#
# function DataStructures.front(s::Stack{T}) where (T <: Any)
#     top(s)
# end

#faster computation of correlation
#basically only computes dot product
function var_cor(v1::GWAS_variant, v2::GWAS_variant)
    ss = 0.0
    n = 0
    @inbounds for i in eachindex(v1.ds)
        if v1.ind[i] && v2.ind[i]
            ss += v1.ds[i]*v2.ds[i]
            n += 1
        end
    end
    cvar = ss/n - v1.caf*v2.caf
    @fastmath v_corr = cvar / sqrt(v1.v*v2.v)
    return v_corr, n
end

# function ldsc_update_r2!(deque, varnow)
#     tot = 0.0
#     @inbounds for var in deque
#         0.05 < var.caf < 0.95 || continue
#         newind = varnow.ind .& var.ind
#         n = sum(newind)
#         r2 = cor(varnow.ds[newind], var.ds[newind])
#         r2_adj = (r2 - (1-r2)/(n-2)) #N-adjusted formula from LDSC paper
#     	var.sumr2 += r2_adj
#     	tot += r2_adj
#     end
#     varnow.sumr2 = tot
#
#     return tot
# end

#update the whole stack
function ldsc_update_r2!(deque, varnow)
    tot = 0.0
    @inbounds for var in deque
        0.05 < var.caf < 0.95 || continue
        r2, n = var_cor(var, varnow)
        r2_adj = (r2 - (1-r2)/(n-2)) #N-adjusted formula from LDSC paper
    	var.sumr2 += r2_adj
    	tot += r2_adj
    end
    varnow.sumr2 = tot
    return tot
end

#Clear deque variants outside window
function ldsc_update_deque_window!(vardeque, pos_now, print_fun; win=1000000)
    while !isempty(vardeque) && (front(vardeque).pos + win < pos_now)
        print_fun(popfirst!(vardeque))
    end
    return 0
end
#Fully clear deque of variants
function ldsc_clear_deque!(vardeque, print_fun)
    while !isempty(vardeque)
    	print_fun(popfirst!(vardeque))
    end
    return 0
end

#LDSC workhorse function
function process_var_ldsc!(varnow, vardeque, print_fun; win=1000000)
    if !isempty(vardeque)
        #If new chromosome, clear deque
        varnow.chrom != front(vardeque).chrom && ldsc_clear_deque!(vardeque, print_fun)

        #Remove variants outside window
        ldsc_update_deque_window!(vardeque, varnow.pos, print_fun)

        #Update variant R2
        ldsc_update_r2!(vardeque, varnow)
    end
    push!(vardeque, copy(varnow))

    return 0
end

###############
#functions with 2 deques
##############

#Clear deque variants outside window
function ldsc_update_deque_window!(vardeque, vardeque_lowmaf, pos_now, print_fun; win=1000000)
    while !isempty(vardeque) && front(vardeque).pos + win < pos_now
        varnow = popfirst!(vardeque)
        while(first(vardeque_lowmaf).pos < varnow.pos)
            print_fun(popfirst!(vardeque_lowmaf))
        end
        print_fun(varnow)
    end
    return 0
end

#Fully clear deque of variants
function ldsc_clear_deque!(vardeque, vardeque_lowmaf, print_fun)
    while !isempty(vardeque)
        varnow = popfirst!(vardeque)
        while !isempty(vardeque_lowmaf) && first(vardeque_lowmaf).pos < varnow.pos
            print_fun(popfirst!(vardeque_lowmaf))
        end
    	print_fun(varnow)
    end
    ldsc_clear_deque!(vardeque_lowmaf, print_fun)
    return 0
end

function process_var_ldsc!(varnow, vardeque, vardeque_lowmaf, print_fun; ldsc_maf = 0.01, win=1000000)
    if !isempty(vardeque)
        #If new chromosome, clear deque
        varnow.chrom != front(vardeque).chrom && ldsc_clear_deque!(vardeque, print_fun)

        #Remove variants outside window
        ldsc_update_deque_window!(vardeque, vardeque_lowmaf, varnow.pos, print_fun)

        #Update variant R2
        if ldsc_maf < varnow.caf < (1-ldsc_maf)
            ldsc_update_r2!(vardeque, varnow)
        end
    end
    if ldsc_maf < varnow.caf < (1-ldsc_maf)
        varnow.sumr2 = 0.0
        push!(vardeque, copy(varnow))
    else
        varnow.sumr2 = NaN
        push!(vardeque_lowmaf, copy(varnow))
    end

    return 0
end
