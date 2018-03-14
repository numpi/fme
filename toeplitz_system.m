function [x, it] = toeplitz_system(am, ap, d1, d2, ni, b, P)
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

	[x,~,~,it] = gmres(@(x) mat_mul1D(am, ap, d1, d2, ni, x), ...
			b, [], 1e-7, 1000, P); 
        
	it = it(2);
end

