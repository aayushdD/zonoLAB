function Z_reduced = Iterative_bounds_rem_function(Zi, varargin)
    p = inputParser;
    addParameter(p, 'tolR', 1e-8, @isnumeric);
    addParameter(p, 'tol', 1e-15, @isnumeric);
    addParameter(p, 'maxIter', 10000, @isnumeric);
    addParameter(p, 'plotResult', false, @islogical);
    addParameter(p, 'tolNonZero', 1e-6, @isnumeric);
    parse(p, varargin{:});
    
    % Initialize variables
    G = Zi.G;
    c = Zi.c;
    A = Zi.A;
    b = Zi.b;
    tolR = p.Results.tolR;
    tol = p.Results.tol;
    tolNonZero = p.Results.tolNonZero;
    tolE = 1e-15;
    maxIter = p.Results.maxIter;

    % [Anew,bnew, Gnew, ~] = refined_rref_with_pivoting(A, b, G);
    % A=Anew;
    % b=bnew;
    % G=Gnew;
    % 
    % Main reduction loop
    for iter = 1:maxIter
        % Remove zero columns/rows
        zeroCols = find(all(abs([G; A]) < tol, 1));
        G(:, zeroCols) = [];
        A(:, zeroCols) = [];
        
        zeroRows = find(all(abs([A, b]) < tol, 2));
        A(zeroRows, :) = [];
        b(zeroRows) = [];
        
        % Refine bounds
        [E, ~] = refine_bounds_function_con_check(A, b);
        
        % Remove factors with equal bounds
        equal_indices = find(E(:, 1) ==E(:, 2));
        if ~isempty(equal_indices)
            c = c + G(:, equal_indices) * E(equal_indices, 1);
            G(:, equal_indices) = [];
            b = b - A(:, equal_indices) * E(equal_indices, 1);
            A(:, equal_indices) = [];
            continue; % Restart loop after removing equal bounds
        end
        
        % Refine bounds again
        [~, R] = refine_bounds_function_con_check(A, b);
        
        % Check for redundant factors
        eff_r = find((R(:, 1) >= -1 + tolR) & (R(:, 2) <= 1- tolR));
        if isempty(eff_r)
            break; % No more redundancies found
        end
        
        % Remove the last redundant factor found
        remove_position = eff_r(end);
        removable_constraint = find(abs(A(:, remove_position)) > tolNonZero, 1);
        
        if isempty(removable_constraint)
            continue; % Skip if no valid constraint found
        end
        removable_constraint = removable_constraint(1);
        % Remove the generator and constraint
        [G, A, b, c] = removeGenerator(G, A, b, c, remove_position, removable_constraint);
    end
    
    % Create output zonotope
    Z_reduced = conZono(G, c, A, b);
    
end

function [G_new, A_new, b_new, c_new] = removeGenerator(G, A, b, c, remove_position, removable_constraint)
%REMOVEGENERATOR Helper function to remove generator and constraint
    [n_c, n_g] = size(A);
    
    % Step 1: Define the E matrix
    E_lam = zeros(n_g, n_c);
    E_lam(remove_position, removable_constraint) = 1;
    
    % Step 2: Compute Lambda matrices
    a_1j = A(removable_constraint, remove_position);
    Lambda_G = G * E_lam / a_1j;
    Lambda_A = A * E_lam / a_1j;
    
    % Steps 3-6: Update matrices
    G_new = G - Lambda_G * A;
    c_new = c + Lambda_G * b;
    A_new = A - Lambda_A * A;
    b_new = b - Lambda_A * b;
    
    % Step 7: Remove the generator and constraint
    G_new(:, remove_position) = [];
    A_new(removable_constraint, :) = [];
    A_new(:, remove_position) = [];
    b_new(removable_constraint) = [];
end