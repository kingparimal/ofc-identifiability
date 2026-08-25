function x = chi2q(p, k)
%CHI2Q  Chi-square quantile. Uses chi2inv when available, else bisects gammainc.
    if exist('chi2inv','file'), x = chi2inv(p,k); return, end
    lo=0; hi=1;
    while gammainc(hi/2, k/2) < p, hi = hi*2; end
    for i=1:200
        m=(lo+hi)/2;
        if gammainc(m/2, k/2) < p, lo=m; else, hi=m; end
    end
    x=(lo+hi)/2;
end
