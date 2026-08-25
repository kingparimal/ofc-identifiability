function c = binocdf_l(k, n, p)
%BINOCDF_L  Binomial CDF without the Statistics Toolbox.
    c = 0;
    for i = 0:k, c = c + nchoosek(n,i)*p^i*(1-p)^(n-i); end
end
