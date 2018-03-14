function [U, V, nrm] = compress_low_rank(U, V, tol)
%COMPRESS_LOW_RANK Recompress a low-rank factorization. 
%
% [U, V, NRM] = COMPRESS_LOW_RANK(U, V, TOL) recovers an optimal (in size)
%     factorization of U * V', up to a certain tolerance TOL. The tolerance
%     is given in the Frobenius norm, relative to the norm of U * V', which
%     is returned in NRM. 
% 


[QU, RU] = qr(U, 0);
[QV, RV] = qr(V, 0);

[U1, S, V1] = svd(RU * RV');

sv = diag(S);
sv = sv(end:-1:1);
rk = sum(cumsum(sv) > sum(sv) * tol);

U = QU * U1(:,1:rk) * sqrt(S(1:rk,1:rk));
V = QV * V1(:,1:rk) * sqrt(S(1:rk,1:rk));

nrm = sum(sv); % S(1,1);

end

