function [am, ap] = symbol_FE2(alpha, n)
%SYMBOL_FE2 Stable computation of the symbol for the FE discretization. 
%
% [AM, AP] = SYMBOL_FE2(ALPHA, N) constructs the symbol of the Finite
%     Element discretization described in the paper "Fast Solvers for 2D
%     Fractional Differential Equations using rank structured matrices", by
%     Massei, Mazza, and Robol. 
%
%     The formulas are based on the paper "A finite element formulation 
%     preserving symmetric and banded diffusion stiffness matrix 
%     characteristics for fractional differential equations" by Lin and
%     Wang, but has been modified to be numerically stable to evaluate. 

omega0=1;
omega1=-4+2^(3-alpha);
omega2=6-2^(5-alpha)+3^(3-alpha);

% We use gamma to avoid repeating 3 - alpha a lot of times in the code. 
gamma = (3 - alpha);

% For indices smaller than this guy, we use the direct formulation.
% Otherwise, we resort to the series expansion. 
jh = 10;

% The formula given in the paper is only good for small j, in our case we
% choose only j<= 10, which should guarantee at least 12 correct digits, if
% I am not mistaken. 
j1 = 3 : jh;
omegaj1 = (j1+1).^(gamma) - 4*j1.^(gamma)+6*(j1-1).^(gamma) ...
    -4*(j1-2).^(gamma) + (j1-3).^(gamma);

% For the other indices (j > 10), we choose to compute the truncated
% series. Only 16 terms are necessary, because of the polynomial decay
% of order at least 1. But we sum 32 nonetheless just to be sure. 
j2 = (jh+1) : (n-2);

% Fractional binomial coefficient
bc = gamma * (gamma - 1) * (gamma - 2) / 6;

% Coefficient of the expansion, depending on s: we might want to actually
% compute them in a stable way, but it does not matter much because when
% this becomes unstable, the powers of j are small enough anyway. 
c = @(s) (1 + 6 * (-1)^s - 4 * (-2)^s + (-3)^s);

% Power of j2, accumulated
j2powers = j2.^(gamma - 3);

omegaj2 =  c(3) * bc * j2powers;
for s  = 4 : 32
    % Update the binomial coefficient and the power of j2
    bc = bc * (gamma - s + 1) / s;
    j2powers = j2powers ./ j2;
    omegaj2 = omegaj2 + c(s) * bc * j2powers;
end

am=[omega1, omega2, omegaj1, omegaj2];
ap=[omega1, omega0];
