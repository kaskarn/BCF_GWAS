function process_var_hwe!(varnow, vcfind)
    Naa, NAa, NAA = mapslices(sum, varnow.pg, dims = 2)
    chi2_score = length(vcfind) * ((4*NAA*Naa - (NAa)^2) / ((2NAA + NAa)*(2Naa+NAa)))^2
    varnow.hwe = 1 - cdf(Chisq(1), chi2_score)
    varnow.n_AA, varnow.n_Aa, varnow.n_aa = NAA, NAa, Naa
end
