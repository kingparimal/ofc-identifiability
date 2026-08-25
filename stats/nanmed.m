function y = nanmed(x)
%NANMED  Median ignoring NaN. Octave's median lacks the 'omitnan' option.
    x = x(~isnan(x));  y = median(x);
end
