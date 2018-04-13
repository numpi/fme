function Experiment1D

    % Values of N in the paper
    Ns = 2.^(13 : 17);
    
    hmoption('threshold', 1e-7);
    hmoption('block-size', 256);
        
    data = zeros(length(Ns), 7);
    
    for experiment = [1, 2]

        if experiment == 1
            alpha = 1.8;
            
            % Select GMRES accuracy to match HODLR one
            gmres_tol = 1e-6;
        else
            alpha = 1.2;
            
            % Select GMRES accuracy to match HODLR one
            gmres_tol = 1e-7;
        end
        for i = 1 : length(Ns)
            n = Ns(i);
            h = 2 / (n+2);
            ni = h^(alpha - 1); % time step is assumed dt equal to h

            t = linspace(h, 2-h, n);
            d1 = gamma(3 - alpha) .* t.^alpha;
            d2 = gamma(3 - alpha) .* (2 - t).^alpha;
            [am, ap] = fractional_symbol(alpha, n);

            A = hm('toeplitz', am, ap, n);
            A = hm('diagonal', ni*ones(n,1)) + hm('diagonal', d1) * A + hm('diagonal', d2) * A';

            xi = t;
            b=h^(alpha)*(-32*(xi.^2+1/8*(2-xi).^2.*(8+xi.^2)-3/(3-alpha)*(xi.^3+(2-xi).^3)+3/((4-alpha)*(3-alpha))*(xi.^4+(2-xi).^4)))';

            if experiment == 1
                P = ( spdiags(reshape(d1 + d2,n,1), 0, n, n) * spdiags(ones(n,1) * [-1 2 -1], -1:1, n, n) + spdiags(ni*ones(n,1), 0, n, n));
            else
                P = spdiags(ones(n,1) * [1 -1], 0:1, n, n);
                P = spdiags(reshape(d1,n,1), 0, n, n) * P + spdiags(reshape(d2,n,1), 0, n, n) * P' + spdiags(ni*ones(n,1), 0, n, n);
            end

            [L, U] = lu(A);
            data(i,6) = timeit(@() lu(A));
            data(i, 4) = timeit(@() U \ (L \ b));
            x = U \ (L \ b);
            data(i, 5) = norm(A*x-b)/norm(x);
            data(i, 1) = timeit(@() toeplitz_system(am, ap, d1, d2, ni, b, P, gmres_tol));
            [x, data(i, 2)] = toeplitz_system(am, ap, d1, d2, ni, b, P, gmres_tol);
            data(i,3) = norm(mat_mul1D(am, ap, d1, d2, ni, x) - b)/norm(x);
            data(i,7) = qsrank(A);
        end
        
        if experiment == 1
            dlmwrite('e1.dat', [Ns', data], '\t');
        else
            dlmwrite('e2.dat', [Ns', data], '\t');
        end
    end
end
