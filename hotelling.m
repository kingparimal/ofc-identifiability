function x = hotelling(p, nRep)
%HOTELLING  95% threshold for a Mahalanobis distance with ESTIMATED covariance.
    if exist('finv','file'), F = finv(0.95, p, nRep-p); else, F = finvbisect(0.95, p, nRep-p); end
    x = (p*(nRep-1)/(nRep-p)) * F;
end
