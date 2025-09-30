function [E]= bound_e_LP_calculator(G,A,b)

    ng = size(A,2);
    E=zeros(ng,2);
    optSolver.lpSolver='gurobi';
    
    for i = 1:ng    
        s = zeros(ng,1);
        s(i) = 1;
        lbs=-ones(size(G,2),1);
        ubs=ones(size(G,2),1);

        % [e_u, ~, ~] = solveLP(s,[],[], A, b, lbs, ubs, optSolver);
        % E(i,2) = e_u(i);
        % 
        % [e_l, ~, ~] = solveLP(-s,[],[], A, b, lbs, ubs, optSolver);
        % E(i,1) = e_l(i);

        % options = optimoptions('linprog','Display','off');
        % [x,~,~,~] = linprog(-s,[],[],A,b,lbs,ubs,options);
        % E(i,1) = x(i);
        % [x,~,~,~] = linprog(s,[],[],A,b,lbs,ubs,options);
        % E(i,2) = x(i);

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
        E(i,2) = x(i);

        model.obj = -s;
        model.lb = lbs;
        model.ub = ubs;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        result = gurobi(model,params);
        x=result.x;
        E(i,1) = x(i);

        lp = Opt('f',-s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        E(i,2) = s'*opt.xopt;
        lp = Opt('f',s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        E(i,1) = s'*opt.xopt;    
    end
end
