function x = finvbisect(p, d1, d2)
    lo = 0; hi = 1;
    while betainc(d1*hi/(d1*hi+d2), d1/2, d2/2) < p, hi = hi*2; end
    for i = 1:200
        m = (lo+hi)/2;
        if betainc(d1*m/(d1*m+d2), d1/2, d2/2) < p, lo = m; else, hi = m; end
    end
    x = (lo+hi)/2;
end
