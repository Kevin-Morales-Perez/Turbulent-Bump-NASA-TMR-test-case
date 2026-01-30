function [psi] = vanLeerLimiter(r)
%Van Leer Limiter

psi=(r + abs(r))/(1 + abs(r));


end