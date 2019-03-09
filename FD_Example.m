function FD_Example

hmoption('threshold', 1e-8);

Ns = 2.^( 9 : 16 );

for k = 1 : 4

	switch mod(k, 2)
		case 1
			beta1 = 1.3;
			beta2 = 1.7;
		case 0
			beta1 = 1.7;
			beta2 = 1.9;
	end

	times = zeros(1, length(Ns));
	ranks = zeros(1, length(Ns));
    qsranks = zeros(1, length(Ns));

	for i = 1 : length(Ns)
		n = Ns(i);
	
		h = 1 / (n+2);
        
        [am1, ap1] = fractional_symbol(beta1, n);
        [am2, ap2] = fractional_symbol(beta2, n);
        
        L1 = hm('toeplitz', am1, ap1, n);
        L2 = hm('toeplitz', am2, ap2, n);
        
        L1 = L1 + L1';
        L2 = L2 + L2';
        
        % This is the same time step that they have in the code Stoll-Breiten-Simoncini
        dt = 1;
        tau1 = -dt / h^beta1; 
        tau2 = -dt / h^beta2;
         
        L1 = tau1 * L1 + .5 * hm('diagonal', ones(n,1));
        L2 = tau2 * L2 + .5 * hm('diagonal', ones(n,1)); 
		
		nrm = max( norm(L1), norm(L2) );
        % nrm = 10 / h;
        
        % Avoid any normalization, for now
        nrm = 1;
		    
		L1 = L1 / nrm;
		L2 = L2 / nrm;
	
		t = linspace(0, 1, n);
		f1 = dt * 100 * sin(10 * pi * t)' / nrm;
		f2 = cos(pi * t)';
        
		tic;
		ranks(i) = 0;
        
        if k <= 2        
            L1s = ek_struct(L1, false);
            L2s = ek_struct(L2, false);
            qsranks(i) = max(hmrank(L1), hmrank(L2));
        else
            pp1 = ones(n, 1); % pp(t', beta1);
            pm1 = ones(n, 1); % pm(t', beta1);
            
            D1 = .5 * tau1 * spdiags(pp1, 0, n, n);
            D2 = .5 * tau1 * spdiags(pm1, 0, n, n);
            
            if beta1 < 1.5  && false % Disabled because P2 is more efficient
                B = spdiags(ones(n,1) * [ 1 -1 ], 0 : 1, n, n);
            else
                B = spdiags(ones(n,1) * [ -1 2 -1 ], -1 : 1, n, n);
            end
            
            [LL1, UU1] = lu(D1 * B + D2 * B' + .5 * speye(n));
            L1s = ek_gmres_struct(@(x) nrm \ mat_mul1D(am1 * tau1, ap1 * tau1, pp1, pm1, .5, x), ...
                @(x) UU1 \ (LL1 \ x), norm(L1));
            
            qp1 = ones(n, 1); % qp(t', beta2);
            qm1 = ones(n, 1); % qm(t', beta2);
            
            D1 = .5 * tau2 * spdiags(qp1, 0, n, n);
            D2 = .5 * tau2 * spdiags(qm1, 0, n, n);
            if beta2 < 1.5 && false % Disabled because P2 is more efficient
                B = spdiags(ones(n,1) * [ 1 -1 ], 0 : 1, n, n);
            else
                B = spdiags(ones(n,1) * [ -1 2 -1 ], -1 : 1, n, n);
            end
            
            [LL2, UU2] = lu(D1 * B + D2 * B' + .5 * speye(n));
            L2s = ek_gmres_struct(@(x) nrm \ mat_mul1D(am2 * tau2, ap2 * tau2, qp1, qm1, .5, x), ...
                @(x) UU2 \ (LL2 \ x), norm(L2));       
        end
        
        Xu = zeros(n, 0); Xv = zeros(n, 0);
        for j = 1 : 8
			f1t = sin(pi * t.') .* 100 * dt * sin(10*j*dt).*t' / nrm;
			f2t = t' .*(1 - t.');
            
            UU = [f1, f1t, Xu / nrm];
            VV = [f2, f2t, Xv ];
            
            [UU, VV] = compress_low_rank(UU, VV, 1e-6);

			[Xu, Xv] = ek_sylv(L1s, L2s, -UU, VV, inf, ...
                @(r, nrm) r < 1e-6 * nrm * n, false, 'fro');
            
            %
            % The following might be used to check if the residual is
            % sufficiently small, if the problem is sufficiently small
            % scale. 
            %
            % norm(L1 * Xu * Xv' + Xu * Xv' * L2' - UU * VV', 'fro') / norm(Xu * Xv', 'fro')
			ranks(i) = max(ranks(i), size(Xu, 2));
        end

		times(i) = toc;
        
		fprintf('N = %d, time = %e, rank = %d, qsrank = %d\n', n, times(i), ranks(i), qsranks(i));
	end
	
	dlmwrite(sprintf('fd-times_%d.dat', k), [ Ns ; times ; ranks ; qsranks ]', '\t');
end

V = [ dlmread('fd-times_1.dat'), dlmread('fd-times_3.dat') ];
dlmwrite('fd-times_13.dat', V, '\t');

V = [ dlmread('fd-times_2.dat'), dlmread('fd-times_4.dat') ];
dlmwrite('fd-times_24.dat', V, '\t');

end

