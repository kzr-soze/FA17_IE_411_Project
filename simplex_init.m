function [istatus,iB,iN,xB] = simplex_init(A,b,c)

iN = [];
    for ab = 1:size(A,2)
        iN(end+1) = ab;
    end
iB = [];
    for ac = size(A,2)+1:size(A,2)+size(b,1)
        iB(end+1) = ac;
    end
Artificial_Variables = [];
    for av = size(A,2)+1:size(A,2)+size(b,1)
        Artificial_Variables(end+1) = av;
    end
c = [zeros(1,size(A,2)), ones(1,size(b,1))]; %M_row
A = [A, eye(size(b,1),size(b,1))];
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
Binv = eye(size(b,1));
xB = b;
irule = 1;
while(1)
    [istatus,iB,iN,xB,Binv] = simplex_step(A,b,c,iB,iN,xB,Binv,irule);
    if(c(iN)-(c(iB)*Binv)*A(:,iN) >=0)
    if(xB >= 0)
        istatus = 0;
    end
    for(chk = 1:size(iB,2))
        for(chkl = 1:size(Artificial_Variables,2))
            if(iB(chk) == Artificial_Variables(chkl) && xB(chk) == 0)
                istatus = 4;
            end
        end
    end
    for(chklm = 1:size(iB,2))
        for(chkln = 1:size(Artificial_Variables,2))
            if(iB(chklm) == Artificial_Variables(chkln) && xB(chklm) ~= 0)
                istatus = 16;
            end
        end
    end
        return
    end
end