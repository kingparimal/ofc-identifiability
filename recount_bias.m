function recount_bias(bestR, r_range, r_true)
%RECOUNT_BIAS  Count best-fit r above / at / below truth, correctly.
%
%   recount_bias(bestR, r_range, r_true)
%
%   bestR   : per-subject best-fit r (each value is an element of r_range)
%   r_range : the swept grid
%   r_true  : the generating value, e.g. theta0(4)
%
%   WHY THIS EXISTS.  A grid built with logspace does not place the truth on a
%   node exactly: logspace(log10(2e-4), log10(2e-2), 15) puts the nominal
%   r = 0.002 node 4.34e-19 ABOVE 0.002.  Comparing bestR > r_true therefore
%   scores every subject who recovered the truth exactly as having overshot it.
%   That is what produced the apparent upward bias in r.
%
%   The fix is to compare against the NODE, r_range(iTrue), not against r_true.
%   Equality is safe because bestR values are drawn from r_range itself.
%   Better still, store the grid INDEX and compare indices; then the failure
%   mode cannot recur.

    bestR = bestR(:);
    [~, iTrue] = min(abs(r_range - r_true));
    node = r_range(iTrue);

    fprintf('--- grid diagnostic ---\n');
    fprintf('  truth %.6e at node %d of %d\n', r_true, iTrue, numel(r_range));
    fprintf('  node - truth   : %.3e\n', node - r_true);
    fprintf('  node == truth  : %d   (0 means the naive comparison is unsafe)\n', ...
            node == r_true);
    fprintf('  nodes below/above truth: %d / %d\n\n', iTrue-1, numel(r_range)-iTrue);

    nA = sum(bestR >  node);
    nE = sum(bestR == node);
    nB = sum(bestR <  node);
    nAbad = sum(bestR > r_true);          % the naive count, for comparison

    fprintf('--- best-fit r vs truth, %d subjects ---\n', numel(bestR));
    fprintf('  CORRECT (vs node) : %d above / %d at / %d below\n', nA, nE, nB);
    fprintf('  NAIVE   (vs %.3e) : %d above  <-- inflated by the %d at truth\n', ...
            r_true, nAbad, nAbad - nA);

    n = nA + nB;
    if n > 0
        k = min(nA, nB);
        p = 2 * sum(arrayfun(@(i) nchoosek(n,i), 0:k)) / 2^n;
        p = min(p, 1);
        fprintf('  two-sided sign test : p = %.4f  (n = %d, ties excluded)\n', p, n);
    else
        fprintf('  two-sided sign test : undefined, every subject at truth\n');
    end
end
