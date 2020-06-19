function [Nprime, u, a_gcv, cost] = mclc_with_gcv(x_sv, v_sv, u_sv, x_leader, v_leader, x_gcv, v_gcv, x_pv, v_pv, x_b, N, tb)
    A = tril(ones(N));
    if x_pv == -1
        x_pv = -1 * ones(N, 1);
    else
        v_pv = v_pv * ones(N, 1);
        x_pv = x_pv + A * v_pv;
    end
    if x_leader == -1
        x_leader = -1 * ones(N, 1);
    else
        v_leader = v_leader * ones(N, 1);
        x_leader = x_leader + A * v_leader;
    end
    x_sv = x_sv * ones(N, 1);
    v_sv = v_sv * ones(N, 1);
    
    a = -2.5:0.5:2.5;
    cost = zeros(N, length(a), 3);
    for n = 1:N
        for i = 1:length(a)
            u_gcv = a(i) * ones(N, 1);
            v_gcv = v_gcv + A * u_gcv;
            x_gcv = x_gcv + A * v_gcv;
            [u, us, ug] = mclc_cost_with_gcv(x_sv, v_sv, u_sv, x_leader, x_gcv, v_gcv, x_pv, x_b, n, a(i), N, tb);
            cost(n, i, :) = [u, us, ug];
        end
    end
    
    filter = cost(:, :, 3) == max(cost(:, :, 3), [], 2);  % max ugcv
    usv = cost(:, :, 2) .* filter;
    filter = usv == max(usv, [], 2);
    usv(~filter) = 1e9;
    filter = usv == min(min(usv));
    [r, c] = find(filter == 1);
    r=r(1);
    c=c(1);
    Nprime = r;
    u = cost(r, c, 1);
    cost = cost(r, c, 2);
    a_gcv = a(c);
end

function [u, us, ug] = mclc_cost_with_gcv(x0_sv, v0_sv, u0_sv, x_leader, x_gcv, v_gcv, x_pv, x_b, n, a, N, tb)
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
    A2 = A * A;
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
    % x >= x_gcv for k >= n
    Agcv = -A * A;
    bgcv = x0_sv + A * v0_sv - x_gcv - v_length - tb * v_gcv;
    Agcv = Agcv(n:N, :);
    bgcv = bgcv(n:N, :);
    
    H = 2 * (Pv * (A' * A) + Pa);
    f = (2 * Pv * (v0_sv - v_des)' * A + Px * E(n, :) * A^2)';

    Aineq = [Av; Au; Ab; Apv; Agcv; Alv];
    bineq = [bv; bu; bb; bpv; bgcv; blv];
    
    lb = umin * ones(N, 1);
    ub = umax * ones(N, 1);
    
    options = optimoptions(@quadprog,'Display','off');
    [u, val, exitflag] = quadprog(H, f, Aineq, bineq, [], [], lb, ub, [], options);
    if exitflag == 1
        us = val + Pv * (v0_sv - v_des)' * (v0_sv - v_des) + Px * E(n, :) * (x0_sv + A * v0_sv);
        x_sv = x0_sv + A * v0_sv + A2 * u;
        ug = ugcv(x_sv, x_gcv, x_pv, a);
        u=u(1);
    else
        u = 0;
        us = 1e9;
        ug = -1;
    end
end

function ug = ugcv(x_sv, x_gcv, x_pv, a)
    beta = 0.7;
    ksp = 1.4;
    ksf = 1.4;
    kacc = 0.256;
    a_max = 4;
    
    x_sg = x_sv - x_gcv;
    x_sg = x_sg / 30;
    
    usp = max(0, min(1, 1 - x_sg));
    usf1 = min(abs(x_sg), 1);
    uacc = -(a / a_max)^2;
    
    if x_pv(1) == -1
        usf2 = ones(length(usf1), 1);
    else
        x_pg = x_pv - x_gcv;
        x_pg = abs(x_pg) / 30;
        usf2 = min(x_pg, 1);
    end
    
    ug = sum(beta * ksp * usp + (1 - beta) * ksf * (usf1 + usf2)) + kacc * uacc * length(x_sv);
end