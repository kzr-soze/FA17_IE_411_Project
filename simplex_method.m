function [istatus,X,eta,iB,iN,xB] = simplex_method(A,b,c,irule)

    [istatus,iB,iN,xB] = simplex_init(A,b,c);
    fprintf('\nsimplex_init complete\n');

    % Terminate if problem is infeasible
    if (istatus == 16)
        istatus = 4;
        X = [];
        eta = 0;
        iB = [];
        iN = [];
        xB = [];
        fprintf('\nProblem is infeasible \n');
        return;
    end
    istatus_2 = istatus;
    iN_2 = iN;
    iN = iN(iN < size(A,2)+1);
    %A = [A, eye(size(b,1),size(b,1))];
    if (b >= 0)
        fprintf('\nb is Positive\n')
    else
        for(zzz = 1:size(b,1))
            if b(zzz) < 0
                A(zzz,:) = -1*A(zzz,:);
                b(zzz) = -b(zzz);
            end
        end
    end
    Binv = inv(A(:,iB));
    c = [c, zeros(1,size(b,1))];
    X = zeros(1,size(iB,2)+size(iN,2));
    istatus = 0;
    while (istatus == 0)
        [istatus,iB,iN,xB,Binv] = simplex_step(A,b,c,iB,iN,xB,Binv,irule);
        X = zeros(1,size(iB,2)+size(iN,2));
        X(iB) = xB;
        eta = c(iB)*xB;
    end
    if(istatus == -1)
        istatus = 0;
    end
    if(istatus == 16)
        istatus = 32;
    end
    if(istatus_2 == 4 & xB >= 0)
        istatus = 16;
    end
    if(istatus_2 == 16)
        istatus = 4;
    end
end