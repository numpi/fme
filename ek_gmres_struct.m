function S = ek_gmres_struct(Afun, Pfun, nrm)
%EK_GMRES_STRUCT Construct a GMRES based struct for rktoolbox. 
%
% S = EK_GMRES_STRUCT(AFUN, PFUN, NRM) constructs a struct compatible with
%     RKTOOLBOX that implements the solver based on the preconditioned
%     GMRES. The function AFUN needs to implement the matrix-vector
%     multiplication, PFUN the application of the preconditioner, and NRM
%     is the norm of A. 
%
% AFUN and PFUN are meant to be passsed to MATLAB's GMRES, in the form 
%
%    gmres(@(x) Afun(x), v, 15, 1e-7, size(v, 1), Pfun);
%

S = struct(); 

S.solve = @(nu, mu, x) gmres_solve(nu, mu, Afun, Pfun, x);
S.multiply = @(rho, eta, x) rho * (Afun(x)) - eta * x;
S.nrm = nrm; 
S.isreal = true;

    function y = gmres_solve(nu, mu, Afun, Pfun, x)
        if nu > mu
            y = x;
            for j = 1 : size(y, 2)
                [y(:,j),~] = gmres(@(x) Afun(x), x(:,j), ...
                    [], 1e-7, 150, Pfun);
            end
            y = nu \ y;
        else
            y = -mu \ x;
        end
    end

end

