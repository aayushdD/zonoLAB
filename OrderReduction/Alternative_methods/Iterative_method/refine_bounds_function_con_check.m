function [E, R] = refine_bounds_function_con_check(A, b)
% REFINE_BOUNDS_FUNCTION - Iteratively compute tighter bounds for constrained zonotope generators
%
% Inputs:
%   A - Constraint matrix (nc x ng)
%   b - Constraint vector (nc x 1)
%
% Outputs:
%   E - Final bounds matrix for generators (ng x 2)
%   R - Final bounds matrix from constraints (ng x 2)

    % Dimensions
    [nc, ng] = size(A); % nc = number of constraints, ng = number of generators

    % Initialize bounds
    E = repmat([-1, 1], ng, 1);  % Initial generator bounds
    R = repmat([-inf, inf], ng, 1);  % Initial constraint-based bounds

    eq=false;

    while ~eq
        E_prev = E;  % Store bounds before this iteration

        % Iterate over constraints and generators
        for i = 1:nc
            for j = 1:ng
                ajj = A(i, j);
                if ajj == 0
                    continue;
                end

                % Contributions from other variables (k ≠ j)
                sum_rest_low = 0;
                sum_rest_high = 0;

                for k = 1:ng
                    if k ~= j
                        R_temp_low = A(i, k) / ajj * E(k, 1);
                        R_temp_high = A(i, k) / ajj * E(k, 2);
                        if R_temp_low > R_temp_high
                            tmp = R_temp_low;
                            R_temp_low = R_temp_high;
                            R_temp_high = tmp;
                        end
                        sum_rest_low = sum_rest_low + R_temp_low;
                        sum_rest_high = sum_rest_high + R_temp_high;
                    end
                end

                % Compute constraint-implied bounds on ξ_j
                Rj = b(i) / ajj - [sum_rest_low, sum_rest_high];
                Rj_lower = min(Rj);
                Rj_upper = max(Rj);

                % Update R and E bounds for generator j
                R(j, 1) = max(R(j, 1), Rj_lower);
                R(j, 2) = min(R(j, 2), Rj_upper);

                E(j, 1) = max(E(j, 1), R(j, 1));
                E(j, 2) = min(E(j, 2), R(j, 2));

                % Ensure lower ≤ upper
                if E(j, 1) > E(j, 2)
                    tmp = E(j, 1);
                    E(j, 1) = E(j, 2);
                    E(j, 2) = tmp;
                end
            end
        end

        % Check for convergence
        tol = 1e-6;
        if all(abs(E(:) - E_prev(:)) < tol)
            eq=true;
        end
     end
end
