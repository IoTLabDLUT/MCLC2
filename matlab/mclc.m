function [Nprime, u, a_gcv, cost] = mclc(x_sv, v_sv, u_sv, x_leader, v_leader, x_gcv, v_gcv, x_pv, v_pv, x_b, N, tb)
    if x_gcv == -1
        [Nprime, u, cost] = mclc_without_gcv(x_sv, v_sv, u_sv, x_leader, v_leader, x_pv, v_pv, x_b, N, tb);
        a_gcv = 0;
    else
        [Nprime, u, a_gcv, cost] = mclc_with_gcv(x_sv, v_sv, u_sv, x_leader, v_leader, x_gcv, v_gcv, x_pv, v_pv, x_b, N, tb);
    end
end