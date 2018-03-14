function v = mat_mul1D(am, ap, d1, d2, ni, x)
%MAT_MUL1D Fast matvec multiplication for the 1D discretization using FD. 
%
% V = MAT_MUL1D(AM, AP, D1, D2, NI, X) performs the matrix vector
%     multplication by the matrix G, defined as follows: 
%
%       G = NI * eye(N) + DIAG(D1) * TOEPLITZ(AM, AP) + ...
%               DIAG(D2) * TOEPLITZ(AP, AM); 
%
%     The multiplication takes O(N \log N) flops, where N is the size of
%     the matrix. The vector containing the symbol are automatically
%     resized if they are shorted than N. 

	n = size(x,1);
	d1 = reshape(d1, n, 1);
	d2 = reshape(d2, n, 1);

	v = ni * x + d1 .* toepmult_fft(am, ap, n, n, x) + d2 .* toepmult_fft(ap, am, n, n, x);
    
end
