function [R]= bound_r_LP_calculator(G,A,b)

    ng = size(A,2);
    R=zeros(ng,2); 

    optSolver.lpSolver='gurobi';
    for i = 1:ng

        s = zeros(ng,1);
        s(i) = 1;
        lbs=-ones(size(G,2),1);
        lbs(i)=-inf;
        ubs=ones(size(G,2),1);
        ubs(i)=inf;

        % [r_u, ~, ~] = solveLP(s,[],[], A, b, lbs, ubs, optSolver);
        % R(i,2) = r_u(i);
        % [r_l, ~, ~] = solveLP(-s,[],[], A, b, lbs, ubs, optSolver);
        % R(i,1) = r_l(i);


        % options = optimoptions('linprog','Display','off');
        % [x,~,~,~] = linprog(-s,[],[],A,b,lbs,ubs,options);
        % R(i,1) = x(i);
        % [x,~,~,~] = linprog(s,[],[],A,b,lbs,ubs,options);
        % R(i,2) = x(i);

        model.sense = '=';
        model.A = sparse(A);
        model.rhs = b;

        model.obj = s;
        model.lb = lbs;
        model.ub = ubs;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        result = gurobi(model,params);
        x=result.x;
        R(i,2) = x(i);

        model.obj = -s;
        model.lb = lbs;
        model.ub = ubs;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        result = gurobi(model,params);
        x=result.x;
        R(i,1) = x(i);


        %Another solver
        % lp = Opt('f',-s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        % opt = mpt_solve(lp); % Solves the LP.
        % R(i,2) = s'*opt.xopt;
        % lp = Opt('f',s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        % opt = mpt_solve(lp); % Solves the LP.
        % R(i,1) = s'*opt.xopt; 


    end
end
