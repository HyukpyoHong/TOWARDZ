function [pos_sol, pos_sol_TF] = CRN_positive_solution_optim(Aeq)
% This function find a coordinate-wise positive solution of a homogeneous
% system of linear equations, Ax = 0.
% If there is a solution, then the variable pos_sol contains the solution
% or 0 otherwise.
% pos_sol_TF represents whether there is (or is not) a solution.

[numrow, numcol] = size(Aeq);
if rank(Aeq) == numcol
    % if the stoichiometric matrix has the full-column rank, there is only
    % the trivial solution. 
    pos_sol_TF = 0;
    pos_sol = 0;
else
    %Aeq = [2 -1];
    
    fun0 = @(x) x;
    x0 = zeros(numcol,1);
    A = [];
    b = [];
    lb = -ones(numcol,1);
    ub = zeros(numcol, 1);
    beq = zeros(numrow, 1);
    nonlcon = [];
    options = optimoptions('fminimax','Display','none');
    x = -fminimax(fun0,x0,A,b,Aeq,beq,lb,ub,nonlcon, options);
    
    if isequal(x, zeros(numcol, 1))
        pos_sol_TF = 0;
        pos_sol = 0;
    else
        pos_sol_TF = 1;
        pos_sol = x;
    end
end
end
