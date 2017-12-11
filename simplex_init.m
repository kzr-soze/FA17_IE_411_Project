function [istatus,iB,iN,xB] = simplex_init(A,b,c)
    
    eps = 10^-10;

    % Initialize iN, iB, and a list of Artificial variables
    iN = zeros(1,size(A,2));
    for ab = 1:size(A,2)
        iN(ab) = ab;
    end
    iB = zeros(1,size(b,1));
    Artificial_Variables = zeros(1,size(b,1));
    for ac = 1:size(b,1)
        iB(ac) = size(A,2)+ac;
        Artificial_Variables(ac) = iB(ac);
    end

    % Construct a new cost vector with positive 1 values for artificial
    % variables, 0 for others, and construct a new corresponding A matrix
    c_aux = [zeros(1,size(A,2)), ones(1,size(b,1))]; %M_row
    A_aux = [A, eye(size(b,1),size(b,1))];
    
    % Check that all b values are non-negative
    if (b >= 0)
        fprintf('\nb is Positive\n');
    else
        for(zzz = 1:size(b,1))
            if b(zzz) < 0
                A_aux(zzz,:) = -1*A_aux(zzz,:);
                A(zzz,:) = -1 *A(zzz,:);
                b(zzz) = -b(zzz);
            end
        end
    end
    
    
    Binv = eye(size(b,1));
    xB = b;
    irule = 1;
    
    while(1)
        [istatus,iB,iN,xB,Binv] = simplex_step(A_aux,b,c_aux,iB,iN,xB,Binv,irule);
        
        % Check for negative reduced costs on the non-basic variables.
        if(c_aux(iN)-(c_aux(iB)*Binv)*A_aux(:,iN) >=0)
            
            istatus = 0;

                       
            % Check that any artificial variables remaining in basis have
            % value of zero.
            
            for (chk = 1:size(iB,2))
                for chkl = 1:size(Artificial_Variables,2)
                    if (iB(chk) == Artificial_Variables(chkl) )
                        if (xB(chk) > eps)
                            istatus = 16;
                            return;
                        end
                    end
                end
            end
            
            X = zeros(1,size(A,2))';
            X(iB)=xB;
            X = X(1:size(A,2));
%             disp(A);
%             disp(X);
%             disp(A*X);
%             disp(b);
            
            err = abs((A*X-b));
            dist = max(err);
            distb = min(xB);
%             disp(distb);
%             disp(dist);
            % Check that all basic variables are non-negative, and that
            % solution is a BFS
            if(distb < -eps || dist > eps)
                istatus = 4;
                return;
            end
            
            return
        end
    end
end