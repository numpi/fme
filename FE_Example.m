function FE_Example

hodlroption('threshold', 1e-8);
show_plot = false;

Ns = 2.^( 9 : 16 );

for k = 1 :  4

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
        
        [am1, ap1] = symbol_FE2(beta1, n);
        [am2, ap2] = symbol_FE2(beta2, n);
        
        L1 = (h^(-beta1) / gamma(4-beta1)) * hodlr('toeplitz', am1, ap1, n);
        L2 = (h^(-beta2) / gamma(4-beta2)) * hodlr('toeplitz', am2, ap2, n);
        
        L1 = L1 + L1';
        L2 = L2 + L2';

	% Mass matrix
        sM = spdiags(ones(n, 1) * [1/6 2/3 1/6], -1:1, n, n);
        M = hodlr('banded', sM, 1, 1);
        
        % Choose a reasonable time step for this problem
        dt = .1;
        
        L1 = dt * (M\L1) - eye(n, 'like', L1) / 2;
        L2 = dt * (M\L2) - eye(n, 'like', L2) / 2;
	
		tic;
		ranks(i) = 0;
        
        if k <= 2        
            L1s = ek_struct(L1, false);
            L2s = ek_struct(L2, false);
            qsranks(i) = max(hodlrrank(L1), hodlrrank(L2));
        else
            pp1 = ones(n, 1); % pp(t', beta1);
            pm1 = ones(n, 1); % pm(t', beta1);
            
            tau1 = dt / h^(beta1) / gamma(4-beta1);
            tau2 = dt / h^(beta2) / gamma(4-beta2);
            nrm = 1;
            
            D1 = .5 * tau1 * spdiags(pp1, 0, n, n);
            D2 = .5 * tau1 * spdiags(pm1, 0, n, n);
            
            if beta1 < .5
                B = spdiags(ones(n,1) * [ 1 -1 ], 0 : 1, n, n);
            else
                B = spdiags(ones(n,1) * [ -1 2 -1 ], -1 : 1, n, n);
            end
            
            [LL1, UU1] = lu(D1 * B + D2 * B' + .5 * speye(n));
	    % Construct symbol of the Toeplitz part of the linear coefficients
	    tam1 = am1 * tau1;
	    tam1(1) = tam1(1) - .25 * sM(1,1); tam1(2) = tam1(2) - .25 * sM(2,1);
            tap1 = ap1 * tau1; tap1(1) = tam1(1); tap1(2) = tap1(2) - sM(1,2) * .25;
            L1s = ek_gmres_struct(@(x) sM \ (nrm \ mat_mul1D(tam1, tap1, pp1, pm1, 0, x)), ...
                @(x) UU1 \ (LL1 \ x), norm(L1));
            
            qp1 = ones(n, 1); % qp(t', beta2);
            qm1 = ones(n, 1); % qm(t', beta2);
            
            D1 = .5 * tau2 * spdiags(qp1, 0, n, n);
            D2 = .5 * tau2 * spdiags(qm1, 0, n, n);
            if beta2 < .5
                B = spdiags(ones(n,1) * [ 1 -1 ], 0 : 1, n, n);
            else
                B = spdiags(ones(n,1) * [ -1 2 -1 ], -1 : 1, n, n);
            end
            
            [LL2, UU2] = lu(D1 * B + D2 * B' + .5 * speye(n));
            tam2 = am2 * tau2;
	        tam2(1) = tam2(1) - sM(1,1) * .25; tam2(2) = tam2(2) - sM(2,1) * .25;
            tap2 = ap2 * tau2; tap2(1) = tam2(1); tap2(2) = tap2(2) - sM(1,2) * .25;

            L2s = ek_gmres_struct(@(x) nrm \ mat_mul1D(tam2, tap2, qp1, qm1, 0, x), ...
                @(x) UU2 \ (LL2 \ x), norm(L2));       
        end
        
        Xu = [ zeros(3*n/8, 1) ;  .5 ; ones(n/4-2,1) ; .5 ; zeros(3*n/8, 1) ];
        Xv = Xu;
        
        f1 = Xu; f2 = Xv;
        
        timesteps = 8;
        
        if show_plot
            t = linspace(0, 1, n);

            if mod(timesteps, 2) == 0
               subplot(2, timesteps / 2, 1);
            else
                subplot(1, timesteps, 1);
            end
        
            [XX, YY] = meshgrid(t, t);
            mesh(XX, YY, Xu * Xv');
	    pause
        end
        
        for j = 1 : timesteps - 1
            % Make sure the RHS is in compressed form
            [UU, VV] = compress_low_rank([dt * f1, -Xu], [f2, Xv], 1e-6);

			[Xu, Xv] = ek_sylv(L1s, L2s, -UU, VV, inf, ...
                @(r, nrm) r < 1e-6 * nrm * n, false, 'fro');
			ranks(i) = max(ranks(i), size(Xu, 2));            
            
            if show_plot
                if mod(timesteps, 2) == 0
                    subplot(2, timesteps / 2, j+1);
                else
                    subplot(1, timesteps, j+1);
                end

                mesh(XX, YY, Xu * Xv')
            end
        end           

		times(i) = toc;

		fprintf('N = %d, time = %e, rank = %d, qsrank = %d\n', n, times(i), ranks(i), qsranks(i));
	end
	
	dlmwrite(sprintf('fe-times_%d.dat', k), [ Ns ; times ; ranks ; qsranks ]', '\t');
end

%V = [ dlmread('fe-times_1.dat'), dlmread('fe-times_3.dat') ];
%dlmwrite('fe-times_13.dat', V, '\t');

end

