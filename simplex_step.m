function [istatus,iB,iN,xB,Binv] = simplex_step(A,b,c,iB,iN,xB,Binv,irule)

%Call using the following format: "[istatus,iB,iN,xB,Binv] =
%   simplex_step(A,b,c,iB,iN,xB,Binv,irule)" after assigning values to
%   required variables. 
%
%Performs a single step of the revised simplex method
%
%   Input Parameters:
%
%A-(m,n) constraint matrix
%b-(m,1) POSITIVE vector appearing in the constraint equation above
%c-(1,n) vector giving the coefficients of the objective function
%iB-(1,m) integer vector specifying the indices of the basic variables at
%         the beginning of the simplex step
%iN-(1,n-m)- integer vector specifying the indices of the nonbasic
%            variables at the beginning of the simplex step
%xB-(m,1) vector specifying the values of the basic variables at the
%         beginning of the simplex step
%Binv-(m,m) inverse matrix of the basis B
%irule-integer parameter specifying which pivot rule to use:
%   irule = 0 indicates that the smalles coefficient rule should be
%   used
%   irule = 1 indicates that Bland's rule should be used
%  
%   Output Parameters:
%
%istatus-integer parameter reporting on the progress or lake thereof made
%        by this function
%   istatus = 0 indicates normal nondegenerate simplex method step
%             completed
%   istatus = 16 indicates the program is unbounded 
%   istatus = -1 indicates an optimal feasible vector has been found
%iB-integer vector specifying the m indices of the basic variables after
%   the simplex step
%iN-integer vector specifying the n-m indices of the nonbasic variables after
%   the simplex step
%xB-vector of length m specifying the values of the basic variables after
%   the simplex step

    if (irule == 0)
        w_T = c(iB)*Binv;
        Negative_w_T = -1*(w_T);
        obj = -1*(c(iB)*Binv*b);
    
        %Prints Initial Tableau of Revised Simplex
        fprintf('\nInitial Tableau\n')
        g = sprintf('x%i ', iB);
        fprintf('z %s RHS\n',g)
        h = sprintf(' %i ', Negative_w_T);
        i = sprintf('%i ', obj);
        fprintf('-1 %s %s \n',h,i)
        for j = 1:size(iB,2)
            k = sprintf(' %i ', Binv(j,:));
            l = sprintf(' %i', xB(j));
            fprintf('0 %s %s\n',k,l)
        end
        %End of Initial Tableau of Revised Simplex
        
        %Start of 1-step Simplex
        Cost_check = [];
        for m = 1:size(iN,2)
            Cost_check(end+1) = c(iN(m))-w_T*A(:,iN(m));
        end
        
    %   index of lowest cost in Cost_check (not necessarily enter_var)
        lowest = Inf;
        for z = 1:size(iN,2)
            if(Cost_check(z) < lowest && Cost_check(z) < 0 && lowest == Inf)
                lowest = Cost_check(z);
                n = z;
            elseif(Cost_check(z) < lowest)
                lowest = Cost_check(z);
                n = z;
            else
                n = n;
            end
        end
        if(Cost_check >= 0)
            istatus = -1;
            return
        end
        enter_var = iN(n); %Read as x[enter_var]
        y = transpose(Binv*A(:,enter_var));
        Min_ratio = [];
        Unbounded_check = [];
        for delta = 1:size(iB,2)
            if(y(delta) < 0)
                Min_ratio(end+1) = Inf;
            else
                Min_ratio(end+1) = xB(delta)/y(delta);
            end
            if(Min_ratio(delta) == Inf || Min_ratio(delta) == -Inf || Min_ratio(delta) < 0 || isnan(Min_ratio(delta)) || y(delta) < 0)
                Unbounded_check(end+1) = -1;
            else
                Unbounded_check(end+1) = 1;
            end
        end
        if(Unbounded_check == -1)
            istatus = 16;
            fprintf('\nThe problem is Unbounded\n')
            return
        end
        lowest_ratio = Inf;
        o = Inf;
        for delta_2 = 1:size(iB,2)
            if(Min_ratio(delta_2) < lowest_ratio && lowest_ratio == Inf && Min_ratio(delta_2) >= 0)
                lowest_ratio = Min_ratio(delta_2);
                o = delta_2;
            elseif((Min_ratio(delta_2) < lowest_ratio) && Min_ratio(delta_2) >= 0)
                lowest_ratio = Min_ratio(delta_2);
                o = delta_2;
            elseif(lowest_ratio == Inf &&(Min_ratio(delta_2) == Inf || isnan(Min_ratio(delta_2))))
                o = Inf;
            else
                o = o;
            end
        end
        exit_var = iB(o); %Read as x[exit_var]
        iB(o) = enter_var;
        iN(n) = exit_var;
        Negative_w_T = -1*(Cost_check(n)/y(o))*Binv(o,:) + Negative_w_T;
        obj = -1*(Cost_check(n)/y(o))*xB(o) + obj;
        for p = 1:size(iB,2)
            if(p ~= o)
                Binv(p,:) = -1*(y(p)/y(o))*Binv(o,:) + Binv(p,:);
                xB(p) = -1*(y(p)/y(o))*xB(o) + xB(p);
            end
        end
        Binv(o,:) = Binv(o,:)/y(o);
        xB(o) = xB(o)/y(o);
        %   Check for Optimal Solution (if so, istatus = -1, otherwise 
        %   return istatus = 0) A check for Unboundedness should be done earlier
    
        %   End of 1-step Simplex
        %   Prints Iteration 1
        fprintf('\nIteration 1\n')
        g = sprintf('x%i ', iB);
        fprintf('z %s RHS\n',g)
        h = sprintf(' %i ', Negative_w_T);
        i = sprintf('%i ', obj);
        fprintf('-1 %s %s \n',h,i)
        for j = 1:size(iB,2)
            k = sprintf(' %i ', Binv(j,:));
            l = sprintf(' %i', xB(j));
            fprintf('0 %s %s\n',k,l)
        end
        %End of Iteration 1
        
        %Start of Status Code
        if(Cost_check >= 0)
            istatus = -1;
        else
            istatus = 0;
        end
    %End of Status Code
    %End of Lowest coefficient code   
    end


    %Beginning of Bland's Rule code
    if(irule == 1)
        w_T = c(iB)*Binv;
        Negative_w_T = -1*(w_T);
        obj = -1*(c(iB)*Binv*b);
        
        %Prints Initial Tableau of Revised Simplex   
        fprintf('\nInitial Tableau\n')
        g = sprintf('x%i ', iB);
        fprintf('z %s RHS\n',g)
        h = sprintf(' %i ', Negative_w_T);
        i = sprintf('%i ', obj);
        fprintf('-1 %s %s \n',h,i)
        for j = 1:size(iB,2)
            k = sprintf(' %i ', Binv(j,:));
            l = sprintf(' %i', xB(j));
            fprintf('0 %s %s\n',k,l)
        end
        
        %End of Initial Tableau of Revised Simplex
        %Start of 1-step Simplex
        Cost_check = [];
        for m = 1:size(iN,2)
            Cost_check(end+1) = c(iN(m))-w_T*A(:,iN(m));
        end
        
        %   index of non-basic variables with negative Cost and smallest index in Cost_check (not necessarily enter_var)
        %   index = Value of smallest index (acts as a comparison; only appears in
        %   the following code.
        %   n = Index number of smallest index in array Cost_check
        %   Example: Cost_check = [0 -12 -1 -.5 0], iN = [5 3 1 2 4]
        %   index = 1, n = 3, so iN(3) = 1 (x1)
        index = max(iN);
        n = size(iN,2);
        for z = 1:size(iN,2)
            if(Cost_check(z) < 0 && iN(z) <= index)
                index = iN(z);
                n = z;
            else
                index = index;
                n = n;
            end
        end
        if(Cost_check >= 0)
            istatus = -1;
            return
        end
        enter_var = iN(n); %Read as x[enter_var]
        y = transpose(Binv*A(:,enter_var));
        Min_ratio = [];
        Unbounded_check = [];
        for delta = 1:size(iB,2)
            if(y(delta) < 0)
                Min_ratio(end+1) = Inf;
            else
                Min_ratio(end+1) = xB(delta)/y(delta);
            end
            if(Min_ratio(delta) == Inf || Min_ratio(delta) == -Inf || Min_ratio(delta) < 0 || isnan(Min_ratio(delta)) || y(delta) < 0)
                Unbounded_check(end+1) = -1;
            else
                Unbounded_check(end+1) = 1;
            end
        end
        if(Unbounded_check == -1)
            istatus = 16;
            fprintf('\nThe problem is Unbounded\n')
            return
        end
    %   index of basic variables with >=0 Min_ratio smallest index in Min_ratio (not necessarily exit_var)
    %   index_ratio = Value of smallest index (acts as a comparison; only appears in
    %   the following code.
    %   o = Index number of smallest index in array Min_ratio with >=0
    %   Min_ratio and smallest Min_ratio; tie in Min_ratio goes to smaller
    %   index
    %   Example: Min_ratio = [-Inf Inf 0 -200 0], iB = [5 4 3 2 1]
    %   index_ratio = 1, o = 5, so iB(5) = 1 (x1)
    %   Example 2: Min_ratio = [Inf -Inf 0 0 0], iB = [1 2 5 4 3]
    %   index_ratio = 3, o = 5, so iB(5) = 3 (x3)
        lowest_ratio = Inf;
        index_ratio = max(iB);
        o = size(iB,2);
        for delta_2 = 1:size(iB,2)
            if(Min_ratio(delta_2) < lowest_ratio && lowest_ratio == Inf && Min_ratio(delta_2) >= 0 && delta_2 == 1)
                lowest_ratio = Min_ratio(delta_2);
                o = delta_2;
                index_ratio = iB(delta_2);
            elseif((Min_ratio(delta_2) <= lowest_ratio) && Min_ratio(delta_2) >= 0 && iB(delta_2) < index_ratio)
                lowest_ratio = Min_ratio(delta_2);
                o = delta_2;
                index_ratio = iB(delta_2);
            elseif((Min_ratio(delta_2) < lowest_ratio) && Min_ratio(delta_2) >= 0)
                lowest_ratio = Min_ratio(delta_2);
                o = delta_2;
                index_ratio = iB(delta_2);    
            elseif(lowest_ratio == Inf &&(Min_ratio(delta_2) == Inf || isnan(Min_ratio(delta_2))))
                o = o;
                index_ratio = index_ratio;
            else
                o = o;
                index_ratio = index_ratio;
            end
        end
        exit_var = iB(o); %Read as x[exit_var]
        iB(o) = enter_var;
        iN(n) = exit_var;
        Negative_w_T = -1*(Cost_check(n)/y(o))*Binv(o,:) + Negative_w_T;
        obj = -1*(Cost_check(n)/y(o))*xB(o) + obj;
        for p = 1:size(iB,2)
            if(p ~= o)
                Binv(p,:) = -1*(y(p)/y(o))*Binv(o,:) + Binv(p,:);
                xB(p) = -1*(y(p)/y(o))*xB(o) + xB(p);
            end
        end
        Binv(o,:) = Binv(o,:)/y(o);
        xB(o) = xB(o)/y(o);
    %   Check for Optimal Solution (if so, istatus = -1, otherwise 
    %   return istatus = 0) A check for Unboundedness should be done earlier
    %End of 1-step Simplex
    %Prints Iteration 1
        fprintf('\nIteration 1\n')
        g = sprintf('x%i ', iB);
        fprintf('z %s RHS\n',g)
        h = sprintf(' %i ', Negative_w_T);
        i = sprintf('%i ', obj);
        fprintf('-1 %s %s \n',h,i)
        for j = 1:size(iB,2)
            k = sprintf(' %i ', Binv(j,:));
            l = sprintf(' %i', xB(j));
            fprintf('0 %s %s\n',k,l)
        end
    %End of Iteration 1
    %Start of Status Code
        if(Cost_check >= 0)
            istatus = -1;
        else
            istatus = 0;
        end
    %End of Status Code
    %End of Bland's Rule code 
    end
end