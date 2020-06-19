function [Nprime, u, cost] = mclc_without_gcv(x_sv, v_sv, u_sv, x_leader, v_leader, x_pv, v_pv, x_b, N, tb)
    umax = 2.5;
    vmax = 13.89;
    if x_pv == -1
        Nprime = 1;
        u = min(umax, vmax - v_sv);
        cost = 0;
        return
    end
    
    A = tril(ones(N));
    
    v_pv = v_pv * ones(N, 1);
    x_pv = x_pv + A * v_pv;
    
    if x_leader == -1
        x_leader = -1 * ones(N, 1);
    else
        v_leader = v_leader * ones(N, 1);
        x_leader = x_leader + A * v_leader;
    end

    x_sv = x_sv * ones(N, 1);
    v_sv = v_sv * ones(N, 1);
    
    cost = zeros(N, 2);
    for n = 1:N
        [u, us] = mclc_cost_without_gcv(x_sv, v_sv, u_sv, x_leader, x_pv, x_b, n, N, tb);
        cost(n, :) = [u, us];
    end
    
    filter = cost(:, 2) == min(cost(:, 2));
    k = find(filter == 1);
    k = k(1);
    Nprime = k;
    u = cost(k, 1);
    cost = cost(k, 2);
end


function [u, us] = mclc_cost_without_gcv(x0_sv, v0_sv, u0_sv, x_leader, x_pv, x_b, n, N, tb)
    vmax = 13.89;
    v_length = 5;
    dumax = 1;
    dumin = -1;
    umax = 2.5;
    umin = -2.5;
    Pv = 3;
    Pa = 1;
    Px = 1;
    v_des = vmax;
    
    A = tril(ones(N));
    E = eye(N);
    % 0 <= v <= vmax
    Av = [-A; A];
    bv = [v0_sv; vmax - v0_sv];
    % dumin <= du <= dumax
    du = -diag(diag(ones(N))) + diag(diag(ones(N), 1), 1);
    du = du(1:end - 1, :);
    Au = [du; -du];
    bu = [dumax * ones(N - 1, 1); -dumin * ones(N - 1, 1)];
    Au = [Au; eye(1, N); -eye(1, N)];
    bu = [bu; dumax + u0_sv; -dumin - u0_sv];
    % x <= x_b - tb * v FOR k < n
    Ab = (A + tb * E) * A;
    bb = x_b - x0_sv - (A + tb * E) * v0_sv;
    Ab = Ab(1:n - 1, :);
    bb = bb(1:n - 1, :);
    % x <= x_leader for k < n
    if x_leader(1) == -1
        Alv = [];
        blv = [];
    else
        Alv = (A + tb * E) * A;
        blv = x_leader - v_length - x0_sv - (A + tb * E) * v0_sv;
        Alv = Alv(1:n - 1, :);
        blv = blv(1:n - 1, :);
    end
    % x <= x_pv for k >= n
    if x_pv(1) == -1
        Apv = [];
        bpv = [];
    else
        Apv = (A + tb * E) * A;
        bpv = x_pv - v_length - x0_sv - (A + tb * E) * v0_sv;
        Apv = Apv(n:end, :);
        bpv = bpv(n:end, :);
    end
    
    H = 2 * (Pv * (A' * A) + Pa);
    f = (2 * Pv * (v0_sv - v_des)' * A + Px * E(n, :) * A^2)';

    Aineq = [Av; Au; Ab; Apv; Alv];
    bineq = [bv; bu; bb; bpv; blv];
    
    lb = umin * ones(N, 1);
    ub = umax * ones(N, 1);
    
    options = optimoptions(@quadprog,'Display','off');
    [u, val, exitflag] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], options);
    if exitflag == 1 
        us = val + Pv * (v0_sv - v_des)' * (v0_sv - v_des) + Px * E(n, :) * (x0_sv + A * v0_sv);
        u=u(1);
    else
        u = 0;
        us = 1e9;
    end
end