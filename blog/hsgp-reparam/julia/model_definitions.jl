HSGP = @slic begin 
    xi = -1 + 2 * (x - min(x)) / (max(x) - min(x))
    L = 1.5
    n_functions = 10
    X = sin(pi/(2L) * (xi+L) * range(1, n_functions)')/sqrt(L)
    log_x_scale ~ normal(0, x_scale_scale)
    x_scale = exp(log_x_scale)
    log_y_scale ~ normal(0, y_scale_scale)
    y_scale = exp(log_y_scale)
    log_scale = (
        -0.25*range(1, n_functions)^2*(x_scale*pi/2L)^2 + 
        rep_vector(log_y_scale + 0.45946926660233633 + .5 * log_x_scale, n_functions)
    )
    unit_weight ~ std_normal(;n=n_functions)
    return X * (exp(log_scale) .* unit_weight)
end
mcycle = @slic begin 
    gp_loc ~ HSGP(;x, n_functions, x_scale_scale, y_scale_scale)
    gp_log_scale ~ HSGP(;x, n_functions, x_scale_scale, y_scale_scale)
    y ~ normal(gp_loc, exp(gp_log_scale))
end

adaptive_HSGP = HSGP(quote 
    adaptive_weight ~ normal(0, exp(log_scale .* adaptive_centeredness); n=n_functions)
    unit_weight = exp(-log_scale .* adaptive_centeredness) .* adaptive_weight
    return X * (exp(log_scale .* (1 - adaptive_centeredness)) .* adaptive_weight)
end)
adaptive_mcycle = @slic begin 
    gp_loc ~ adaptive_HSGP(;x, n_functions, adaptive_centeredness=gp_loc_centerednesses, x_scale_scale, y_scale_scale)
    gp_log_scale ~ adaptive_HSGP(;x, n_functions, adaptive_centeredness=gp_log_scale_centerednesses, x_scale_scale, y_scale_scale)
    y ~ normal(gp_loc, exp(gp_log_scale))
end




mcycle10 = @slic begin 
    gp_loc ~ HSGP(;x, n_functions, x_scale_scale, y_scale_scale)
    gp_log_scale ~ normal(0, y_scale_scale)
    y ~ normal(gp_loc, exp(gp_log_scale))
end
adaptive_mcycle10 = @slic begin 
    gp_loc ~ adaptive_HSGP(;x, n_functions, adaptive_centeredness=gp_loc_centerednesses, x_scale_scale, y_scale_scale)
    gp_log_scale ~ normal(0, y_scale_scale)
    y ~ normal(gp_loc, exp(gp_log_scale))
end

for model in (mcycle, adaptive_mcycle, mcycle10, adaptive_mcycle10)
    model.data[:docstring] = """
    Generated with [StanBlocks.jl v0.1.5](https://github.com/nsiccha/StanBlocks.jl).
    """
end