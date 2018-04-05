function [x, it] = toeplitz_system(am, ap, d1, d2, ni, b, P, gmres_tol)
%TOEPLITZ_SYSTEM Solve a Toeplitz linear system by means of PGMRES. 
%
% [X, IT] = TOEPLITZ_SYSTEM(AM, AP, D1, D2, NI, B, P) solves the linear
%     G X = B where the matrix G is given by 
%
%       G = NI * eye(N) + DIAG(D1) * TOEPLITZ(AM, AP) + ...
%               DIAG(D2) * TOEPLITZ(AP, AM); 
%
%     The solution X is returned along with the number of iterations
%     required. 

    if ~exist('gmres_tol', 'var')
        gmres_tol = 1e-8;
    end

	[x,~,~,it] = gmres(@(x) mat_mul1D(am, ap, d1, d2, ni, x), ...
			b, [], gmres_tol, 1000, P); 
        
	it = it(2);
end

