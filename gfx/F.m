function y=F(t,T,par,dt, sigma)
    %sigma=(4*par.K2*par.Ts*par.RTs)
    %sigma=sqrt(sigma);
    %sigma=0;
%    sigma=0.004;
    
    %sigma=1e-3/6;
    
    y=1/(par.RTs*exp(par.alpha*(T-par.Ts)))...
        *(par.Vb(t)^2+2*(dt)^(-1/2)*sigma*randn*par.Vb(t)+sigma^2)...
        +par.e*(par.Pt+par.Ps)-(2*par.A)*par.e*par.K*T^4-par.Gleg*(T-par.Ts);
    y=y/par.C;
end