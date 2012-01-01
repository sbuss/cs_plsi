function [pZ, pZc, pWgZ, pWgZc, pDgZ, pDgZc, l, lc, t, tc] = cs_plsi_demo(B,k, A)
    beta = B;

    if nargin < 3
        P = UFget(1172);
        A = P.A;
        A = A(:,sum(A,1)~=0);
        A = A(sum(A,2)~=0,:);
    end
    pZ = rand(k,1);
    pWgZ = rand(size(pZ,1), size(A,1)); 
    pDgZ = rand(size(pZ,1), size(A,2));
    pWgZ(1,33)
    size(pWgZ)
    size(pDgZ)
    size(A)
    
    figure(1);
    tic;
    [pWgZc, pDgZc, pZc, lc] = mockSparsePlsi(full(A), k, B, pWgZ', pDgZ, pZ);
    tc = toc;
    %[pWgZc, pDgZc, pZc, pZgWD, lc] = plsi(full(A), k, B, pWgZ', pDgZ, pZ);
    

    ll = zeros(200,1);
    tic;
    for i=(1:200)
        [pZo,pWgZo,pDgZo,l] = cs_plsi(A,size(pZ,1),beta, pZ, pWgZ, pDgZ);
        if l == 100 || l == -1
            break;
        end
        ll(i) = l;
        figure(2);
        pZ = pZo;
        pWgZ = pWgZo;
        pDgZ = pDgZo;
        %subplot(1,2,1);
        cspy(pWgZ' * diag(pZ) * pDgZ);
        %subplot(1,2,2);
        %plot(ll);
        %keyboard;
    end
    t = toc;
end

