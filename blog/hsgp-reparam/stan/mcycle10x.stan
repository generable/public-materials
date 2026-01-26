// Generated with [StanBlocks.jl v0.1.5](https://github.com/nsiccha/StanBlocks.jl). 
functions {
vector normal_lpdfs(
    vector obs,
    vector loc,
    real scale
) {
    int n = dims(obs)[1];
    return jbroadcasted_normal_lpdfs(obs, loc, scale);
}
vector jbroadcasted_normal_lpdfs(
    vector x1,
    vector x2,
    real x3
) {
    int n = dims(x1)[1];
    vector[n] rv;
    for(i in 1:n) {
        rv[i] = normal_lpdfs(broadcasted_getindex(x1, i), broadcasted_getindex(x2, i), broadcasted_getindex(x3, i));
    }
    return rv;
}
real normal_lpdfs(
    real args1,
    real args2,
    real args3
) {
    return normal_lpdf(args1 | args2, args3);
}
real broadcasted_getindex(vector x, int i) {
    int m = dims(x)[1];
    return x[i];
}
real broadcasted_getindex(real x, int i) {
    return x;
}
}
data {
    int x_n;
    vector[x_n] x;
    int n_functions;
    int x_scale_scale;
    int y_scale_scale;
    int y_n;
    vector[y_n] y;
}
transformed data {
    vector[x_n] gp_loc_xi = (-1 + ((2 * (x - min(x))) / (max(x) - min(x))));
    real gp_loc_L = 1.5;
    matrix[x_n, n_functions] gp_loc_X = (
        sin(
            (
                (3.141592653589793 / (2 * gp_loc_L)) *
                (gp_loc_xi + gp_loc_L) *
                (linspaced_vector(n_functions, 1, n_functions)')
            )
        ) /
        sqrt(gp_loc_L)
    );
}
parameters {
    real gp_loc_log_x_scale;
    real gp_loc_log_y_scale;
    vector[n_functions] gp_loc_unit_weight;
    real gp_log_scale;
}
transformed parameters {
    real gp_loc_x_scale = exp(gp_loc_log_x_scale);
    vector[n_functions] gp_loc_log_scale = (
        (
            -0.25 *
            (linspaced_vector(n_functions, 1, n_functions) ^ 2) *
            (((gp_loc_x_scale * 3.141592653589793) / (2 * gp_loc_L)) ^ 2)
        ) +
        rep_vector((gp_loc_log_y_scale + 0.45946926660233633 + (0.5 * gp_loc_log_x_scale)), n_functions)
    );
    vector[x_n] gp_loc = (gp_loc_X * (exp(gp_loc_log_scale) .* gp_loc_unit_weight));
}
model {
    gp_loc_log_x_scale ~ normal(0, x_scale_scale);
    gp_loc_log_y_scale ~ normal(0, y_scale_scale);
    gp_loc_unit_weight ~ std_normal();
    gp_log_scale ~ normal(0, y_scale_scale);
    y ~ normal(gp_loc, exp(gp_log_scale));
}
generated quantities {
    real gp_loc_y_scale = exp(gp_loc_log_y_scale);
    vector[y_n] y_likelihood = normal_lpdfs(y, gp_loc, exp(gp_log_scale));
    vector[x_n] y_gen = to_vector(normal_rng(gp_loc, exp(gp_log_scale)));
}