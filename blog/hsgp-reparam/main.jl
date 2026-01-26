using BridgeStan,
    DataFrames,
    JSON,
    LinearAlgebra,
    Plots,
    RDatasets,
    Random,
    StanBlocks,
    StanLogDensityProblems,
    Statistics,
    Term,
    WarmupHMC,
    DataFrames, 
    PrettyTables

include("julia/convenience_functions.jl")
include("julia/model_definitions.jl")

data = dataset("MASS", "mcycle")
y_scale = std(data.Accel)
y = data.Accel ./ y_scale
n_functions = 20
x_scale_scale = y_scale_scale = 4
centerednesses = range(0, 1, 101)
plot_idxs = [1, 2, n_functions-1, n_functions]
plot_idx_labels = "basis function #" .* string.(plot_idxs')
xi = range(-1.5, +1.5, 501)
didxs = "(numbers 1, 2, 19, and 20)"
dcris = "Shaded areas represent 90, 80, and 50% central credible intervals."

gp_basis_function(xi, idx; L=1.5) = sin(pi/2L * (xi+L) * idx) / sqrt(L)
gp_log_scale(log_x_scale, log_y_scale, idx; L=1.5) = (
    -0.25*(idx*exp(log_x_scale)*pi/2L)^2 .+ 
    (log_y_scale + 0.45946926660233633 + .5 * log_x_scale)
)

# start snippet loss_function
"""
Given posterior draws of

* non-centered parameters (basis function weights in our HSGP example) and
* log scales (usually the scale parameter of the prior of that parameter)

and a candidate centeredness parameter (between 0 and 1), computes the loss function, 
i.e. the KL-divergence of a normal approximation of the reparametrized posterior
from that reparametrized posterior.
"""
loss_function(noncentered_weights::AbstractVector, log_scales::AbstractVector, centerednesses::Real) = begin 
    transformed_weights = noncentered_weights .* exp.(log_scales .* centerednesses)
    log(std(transformed_weights)) - mean(log_scales .* centerednesses)
end
# end snippet loss_function
gps(df) = map(("gp_loc", "gp_log_scale")) do stem
    value = df[:, Regex("$(stem)\\.")] |> Matrix
    length(value) == 0 && return
    log_x_scale = df[:, "$(stem)_log_x_scale"] |> Vector
    log_y_scale = df[:, "$(stem)_log_y_scale"] |> Vector
    log_scale = df[:, Regex("$(stem)_log_scale\\.")] |> Matrix
    unit_weight = df[:, Regex("$(stem)_unit_weight\\.")] |> Matrix
    adaptive_weight = df[:, Regex("$(stem)_adaptive_weight\\.")] |> Matrix
    centered_weight = unit_weight .* exp.(log_scale)
    loss_functiones = loss_function.(eachcol(unit_weight), eachcol(log_scale), centerednesses')
    optimal_centerednesses = centerednesses[argmin.(eachrow(loss_functiones))]
    optimal_weight = unit_weight .* exp.(log_scale .* optimal_centerednesses')
    (;
        value,
        log_x_scale,
        log_y_scale,
        log_scale,
        unit_weight,
        adaptive_weight,
        centered_weight,
        loss_functiones,
        optimal_centerednesses,
        optimal_weight,
    )
end
problem = stan_instantiate(
    mcycle(;x=data.Times, y, n_functions, x_scale_scale, y_scale_scale); 
    path="stan/mcyclex.stan"
)
fit = @isdefined(fit) ? fit : WarmupHMC.adaptive_warmup_mcmc(
    Xoshiro(1), problem; progress=Term.ProgressBar, n_draws=10_000
)
gp_loc, gp_log_scales = gps(drawsdf(problem, fit.posterior_position))

adaptive_problem = stan_instantiate(
    adaptive_mcycle(;
        x=data.Times, y, n_functions, x_scale_scale, y_scale_scale,
        gp_loc_centerednesses=gp_loc.optimal_centerednesses,
        gp_log_scale_centerednesses=gp_log_scales.optimal_centerednesses
    ); 
    path="stan/adaptive_mcyclex.stan"
)
adaptive_fit = @isdefined(adaptive_fit) ? adaptive_fit : WarmupHMC.adaptive_warmup_mcmc(
    Xoshiro(1), adaptive_problem; progress=Term.ProgressBar, n_draws=10_000
)
agp_loc, agp_log_scales = gps(drawsdf(adaptive_problem, adaptive_fit.posterior_position))

myplot([
    vline!(
        plot(
            xi, 
            gp_basis_function.(xi, plot_idxs'); 
            color=plot_idxs', alpha=ternary.(abs.(xi) .< 1, 1, .5), xticks=([-1.5, -1, 0, +1, 1.5], ["-L", -1, 0, +1, "+L"]),
            label=plot_idx_labels, xlabel="unit x", ylabel="basis function value", linewidth=3
        ),
        [-1, +1]; color=:black, label="", linewidth=3, z_order=:back
    )
    vline!(
        plot(
            exp.(gp_log_scale.([-2 -1 0], 0, 1:n_functions)); marker=:circle, color=palette([:purple, :brown], 3) |> collect |> permutedims,
            label="GP log length scale = " .* string.([-2 -1 0]), xlabel="basis function number", ylabel="weight prior scale"
        ),
        plot_idxs, xticks=plot_idxs, label="", color=plot_idxs, linewidth=3, z_order=:back
    )
] |> permutedims) |> savefig("_figs/hsgp_viz.png"; caption="""
HSGP Vizualization. 
Left: shapes of select basis functions $didxs. The shapes of the basis functions are independent of the GP hyperparameters. 
Black vertical lines mark the end of the domain of good approximation, x-ticks at +/- L mark the discretization parameter L.
Right: scale of the zero-mean normal prior of the individual basis function weights, as a function of the basis function number (x-axis) and GP length scale (color). GPs with smaller length scales put more prior weight on higher-frequency basis functions. Vertical lines mark the basis functions visualized on the left (and throughout this post).
""")
myplot([
    qplot!(scatter(data.Times, y; color=:black, label="observations"), data.Times, gp_loc.value', qs; label="inferred mean", ),
    qplot!(plot(), data.Times, exp.(gp_log_scales.value)', qs, label="inferred scale", yscale=:log10)
] |> permutedims; xlabel="time (milliseconds after impact)", ylabel="scaled acceleration") |> savefig("_figs/hsgp_posterior.png"; caption="""
HSGP posteriors. 
$dcris 
Left: observations (black dots) and posterior of the mean function GP. 
Right: posterior of the scale function GP. 
""")

myplot([
    [
        scatter(
            exp.(gp_loc.log_x_scale), gp_loc.unit_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP marginal std.", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_loc.log_y_scale), gp_loc.unit_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP length scale", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_log_scales.log_x_scale), gp_log_scales.unit_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP marginal std.", ylabel="scale GP weight #$idx"
        )  scatter(
            exp.(gp_log_scales.log_y_scale), gp_log_scales.unit_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP length scale", ylabel="scale GP weight #$idx"
        ) 
    ]
    for idx in plot_idxs
]; link=:x, xscale=:log10, markerstrokewidth=0, label="") |> savefig("_figs/unit_scatter.png"; caption="""
Dependency of non-centered basis function weights on GP hyperparameters.
Posterior pair plots of select $didxs **non-centered** basis function weights against the GP hyperparameters.
All quantities are plotted such that they represent the sampler geometry.
Row-wise (from top to bottom): Different basis functions $didxs.
Col-wise (from left to right): Marginal standard deviation and length scale of the mean function GP, and marginal standard deviation and length scale of the scale function GP.
""")

myplot([
    [
        scatter(
            exp.(gp_loc.log_x_scale), gp_loc.centered_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP marginal std.", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_loc.log_y_scale), gp_loc.centered_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP length scale", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_log_scales.log_x_scale), gp_log_scales.centered_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP marginal std.", ylabel="scale GP weight #$idx"
        )  scatter(
            exp.(gp_log_scales.log_y_scale), gp_log_scales.centered_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP length scale", ylabel="scale GP weight #$idx"
        ) 
    ]
    for idx in plot_idxs
]; link=:x, xscale=:log10, markerstrokewidth=0, label="") |> savefig("_figs/centered_scatter.png"; caption="""
Dependency of centered basis function weights on GP hyperparameters.
Posterior pair plots of select $didxs **centered** basis function weights against the GP hyperparameters.
All quantities are plotted such that they represent the sampler geometry.
Row-wise (from top to bottom): Different basis functions $didxs.
Col-wise (from left to right): Marginal standard deviation and length scale of the mean function GP, and marginal standard deviation and length scale of the scale function GP.
""")

myplot([
    [
        scatter(
            exp.(gp_loc.log_x_scale), gp_loc.optimal_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP marginal std.", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_loc.log_y_scale), gp_loc.optimal_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP length scale", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(gp_log_scales.log_x_scale), gp_log_scales.optimal_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP marginal std.", ylabel="scale GP weight #$idx"
        )  scatter(
            exp.(gp_log_scales.log_y_scale), gp_log_scales.optimal_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP length scale", ylabel="scale GP weight #$idx"
        ) 
    ]
    for idx in plot_idxs
]; link=:x, xscale=:log10, markerstrokewidth=0, label="") |> savefig("_figs/optimal_scatter.png"; caption="""
Dependency of optimally parametrized basis function weights on GP hyperparameters.
Posterior pair plots of select $didxs **optimally parametrized** basis function weights against the GP hyperparameters.
All quantities are plotted such that they represent the sampler geometry.
Row-wise (from top to bottom): Different basis functions $didxs.
Col-wise (from left to right): Marginal standard deviation and length scale of the mean function GP, and marginal standard deviation and length scale of the scale function GP.
""")

myplot([
    plot(
        centerednesses, scaledto.(eachrow(gp_loc.loss_functiones)[plot_idxs]);
        label=plot_idx_labels, color=plot_idxs', xlabel="partial centeredness (mean function GP weights)", linewidth=3
    ) plot(
        centerednesses, scaledto.(eachrow(gp_log_scales.loss_functiones)[plot_idxs]);
        label=plot_idx_labels, color=plot_idxs', xlabel="partial centeredness (scale function GP weights)", linewidth=3
    )
]; ylabel="rescaled loss", link=:y, ylim=[0,1]) |> savefig("_figs/loss_profile.png"; caption="""
Loss profiles. 
Profiles of select $didxs basis function specific loss functions, rescaled to [0, 1].
Left: basis function weights of mean function GP.
Right: basis function weights of scale function GP.

""")

myplot([
    vline!(
        plot(
            [gp_loc.optimal_centerednesses, gp_log_scales.optimal_centerednesses];
            label=["mean function GP" "scale function GP"], 
            xlabel="basis function number", ylabel="optimal centeredness", ylim=[0,1], linewidth=3, marker=:circle
        ),
        plot_idxs, xticks=plot_idxs, label="", color=plot_idxs, linewidth=3, z_order=:back
    )
]) |> savefig("_figs/wave_profile.png"; caption="""
Optimal centerdnesses profiles. The optimal centerednesses of the basis function weights for the two GPs (color) as a function of the basis function number (x-axis).
""")


myplot([
    qplot!(scatter(data.Times, y; color=:black, label="observations"), data.Times, agp_loc.value', qs; label="inferred mean", ),
    qplot!(plot(), data.Times, exp.(agp_log_scales.value)', qs, label="inferred scale", yscale=:log10)
] |> permutedims; xlabel="time (milliseconds after impact)", ylabel="scaled acceleration") |> savefig("_figs/ahsgp_posterior.png"; caption="""
Optimally resampled HSGP posteriors. 
$dcris 
Left: observations (black dots) and posterior of the mean function GP. 
Right: posterior of the scale function GP. 
""")

myplot([
    [
        scatter(
            exp.(agp_loc.log_x_scale), agp_loc.adaptive_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP marginal std.", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(agp_loc.log_y_scale), agp_loc.adaptive_weight[:, idx];
            alpha=.1, color=idx, xlabel="mean GP length scale", ylabel="mean GP weight #$idx"
        ) scatter(
            exp.(agp_log_scales.log_x_scale), agp_log_scales.adaptive_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP marginal std.", ylabel="scale GP weight #$idx"
        )  scatter(
            exp.(agp_log_scales.log_y_scale), agp_log_scales.adaptive_weight[:, idx];
            alpha=.1, color=idx, xlabel="scale GP length scale", ylabel="scale GP weight #$idx"
        ) 
    ]
    for idx in plot_idxs
]; link=:x, xscale=:log10, markerstrokewidth=0, label="") |> savefig("_figs/adaptive_scatter.png"; caption="""
Dependency of optimally resampled basis function weights on GP hyperparameters.
Posterior pair plots of select $didxs **optimally resampled** basis function weights against the GP hyperparameters.
All quantities are plotted such that they represent the sampler geometry.
Row-wise (from top to bottom): Different basis functions $didxs.
Col-wise (from left to right): Marginal standard deviation and length scale of the mean function GP, and marginal standard deviation and length scale of the scale function GP.
""")



begin



homoscedastic_problem = stan_instantiate(
    mcycle10(;x=data.Times, y, n_functions, x_scale_scale, y_scale_scale); 
    path="stan/mcycle10x.stan"
)
homoscedastic_fit = @isdefined(homoscedastic_fit) ? homoscedastic_fit : WarmupHMC.adaptive_warmup_mcmc(
    Xoshiro(1), homoscedastic_problem; progress=Term.ProgressBar, n_draws=10_000
)
homoscedastic_gp_loc, _ = gps(drawsdf(homoscedastic_problem, homoscedastic_fit.posterior_position))

homoscedastic_adaptive_problem = stan_instantiate(
    adaptive_mcycle10(;
        x=data.Times, y, n_functions, x_scale_scale, y_scale_scale,
        gp_loc_centerednesses=homoscedastic_gp_loc.optimal_centerednesses,
    ); 
    path="stan/adaptive_mcycle10x.stan"
)

target_acceptance_rates = [.8, .9, .99]
n_chains = 40
for (key, (ncp, acp)) in [
    :heteroscedastic => (problem, adaptive_problem), 
    :homoscedastic => (homoscedastic_problem, homoscedastic_adaptive_problem),
]
    tmp = [
            WarmupHMC.adaptive_warmup_mcmc(
                Xoshiro.(1:n_chains), p; 
                progress=Term.ProgressBar, n_draws=1_000, target_acceptance_rate,
                ntries=100
            )
            for target_acceptance_rate in target_acceptance_rates, p in (ncp, acp)
        ]
 
    hcat(
        DataFrame("target acceptance rate"=>target_acceptance_rates),
        map(1:2) do i 
            prefix = ["non centered", "partially centered"][i]
            map(tmp[:, i]) do fits
                ess = WarmupHMC.MCMCDiagnosticTools.ess(permutedims(stack(getproperty.(fits, :posterior_position)), (2, 3, 1))) |> (x->trunc(Int, minimum(x)))
                n_evaluations = sum(fit->getproperty(fit, :total_evaluation_counter), fits)
                d = Dict()
                for fit in fits
                    d[fit.n_divergent_samples] = get(d, fit.n_divergent_samples, 0) + 1
                end
                n_divergent_samples = string(sum(fit->getproperty(fit, :n_divergent_samples), fits)) * " (" * join([
                    "$(v)x$k"
                    for (k, v) in sort(pairs(d))
                ], "+") * ")"
                (;ess, n_evaluations, n_divergent_samples)
            end |> (x->rename(DataFrame(x), 
                :ess=>"$prefix ESS", 
                :n_evaluations=>"$prefix #evals", 
                :n_divergent_samples=>"$prefix #divergent"
            ))
        end...
    ) |> savetable("_figs/$key-sampling.md"; formatters=[(v, i, j)->WarmupHMC.short_string(v)], caption="""
Sampling statistics for the $key model. $n_chains independent chains with 1000 post-warmup draws each are used to sample from the posterior for varying target acceptance rates (left-most column, 0.8, 0.9 and 0.99). We report the minimal ESS (2nd and 5th columns), the total number of gradient evaluations (3rd and 6th columns), and the number of divergences (as "total (histogram over chains)", 4th and 7th columns) for the non-centered parametrizations (columns 2, 3, and 4) and the partially centered parametrizations (columns 5, 6, and 7). While the individual chains are independent, all partially centered chains use the same optimized parametrization due to the initial long chain from the main text (if apllicable) with 10,000 post warm-up draws.
""")
    end
end