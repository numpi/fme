function [am, ap] = fractional_symbol(alpha, n)
%FRACTIONAL_SYMBOL Construct the symbol of the Grunwald-Letkinov derivative
%
% [AM, AP] = FRACTIONAL_SYMBOL(ALPHA, N) construct the negative and
%     positive parts of the symbol of the Toeplitz matrix discretizing the
%     fractional derivative by means of the Grunwald-Letkinov shifted
%     formulas. 
%

v = zeros(n+2, 1);
v(1) = 1;

v = -cumprod([1,1-((alpha + 1)./(1:n))]);
am = v(2:end);
ap = [v(2), v(1)];

end
