function [pZo, pWgZo, pDgZo] = cs_plsi (X,k,B,pZ,pWgZ,pDgZ)                  %#ok
%CS_PLSI sparse PLSI
%   [pZ, pWgZ, pDgZ] = cs_plsi (X,k,B,pZ,pWgZ,pDgZ) 
%
%   Example:
%       P = UFget(165);
%       A = P.A;
%       [m,n] = size(A);
%       k = 10;
%       pZ = rand(10,1); pWgZ = rand(m,k); pDgZ = rand(n,k);
%       [pZo, pWgZo, pDgZo] = cs_plsi(A,k,0.9,pZ,pWgZ,pDgZ);
%

%   cs_plsi uses CSparse, which is:
%   Copyright 2006-2007, Timothy A. Davis.
%   http://www.cise.ufl.edu/research/sparse

error ('cs_plsi mexFunction not found') ;
