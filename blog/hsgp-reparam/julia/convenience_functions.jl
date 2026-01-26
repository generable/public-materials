drawsdf(posterior, draws; seed=0, include_tp=true, include_gq=true, nan_on_error=false) = begin
    stan_rng = StanRNG(posterior.model, seed)
    names = param_names(posterior.model; include_tp, include_gq)
    d, n_draws = size(draws)
    tmp_draw = zeros(d)
    tmp_rv = zeros(length(names))
    rv_matrix = zeros((length(names), n_draws))
    for (i, draw) in enumerate(eachcol(draws))
        tmp_draw .= draw
        try
            param_constrain!(posterior.model, tmp_draw, tmp_rv; include_tp, include_gq, rng=stan_rng)
        catch
            nan_on_error || rethrow()
            tmp_rv .= NaN
        end
        rv_matrix[:, i] .= tmp_rv
    end
    DataFrame(rv_matrix', names)
end
const qs = [.05, .1, .25, .5]
qplot!(p, x::AbstractVector, Y::AbstractVector{<:AbstractVector}; q=.025, medianalpha=1., fillalpha=.25, label="", kwargs...) = begin 
    plot!(p, x, quantile.(Y, .5); fillrange=quantile.(Y, .5), alpha=medianalpha, fillalpha=medianalpha, label, z_order=:back, kwargs...)
    plot!(
        p, 
        x, quantile.(Y, q); fillrange=quantile.(Y, 1-q), linewidth=0, fillalpha, label="", z_order=:back, kwargs...
    )
end
qplot!(p, x::AbstractVector, Y::AbstractMatrix; alpha=.25, n=0, kwargs...) = begin 
    n > 0 && plot!(p, x, rand(eachcol(Y), n); alpha, label="", kwargs...)
    qplot!(p, x, eachrow(Y); kwargs...)
end
qplot!(p, x, Y, qs; n=0, kwargs...) = begin 
    for q in qs
        qplot!(p, x, Y; q, n, color="#0B7BEC", kwargs...)
        kwargs = merge((;kwargs...), (;label=""))
    end
    p
end
myplot(plots::AbstractVector{<:AbstractVector}; kwargs...) = myplot(reduce(hcat, plots); kwargs...)
myplot(plots::AbstractVector{<:AbstractMatrix}; kwargs...) = myplot(reduce(vcat, plots); kwargs...)
myplot(plots::AbstractVector; kwargs...) = myplot(reshape(plots, (:, 1)); kwargs...)
myplot(plots::AbstractMatrix; 
    scale=(.75, 3 / size(plots, 2)),
    col_width=400, row_height=400, 
    heights=fill(1, size(plots, 1)), widths=fill(1, size(plots, 2)),
    margin=tuple.((2, 2, 6, 6), :mm), 
    kwargs...
) = plot(
    permutedims(plots)...; 
    layout=Plots.grid(
        size(plots)..., 
        heights=length(heights) == 1 ? nothing : normalize(heights, 1), 
        widths=length(widths) == 1 ? nothing : normalize(widths, 1)
    ),
    size=reverse(scale .* (row_height * sum(heights), col_width * sum(widths))), 
    top_margin=margin[1], right_margin=margin[2], bottom_margin=margin[3], left_margin=margin[4],
    kwargs...
)

doandpass(f, x) = (f(x); x)
doandpass(f) = Base.Fix1(doandpass, f)
savefig(f::Function; path, force=true, caption=nothing, kwargs...) = if force || !isfile(path)
    mkpath(dirname(path))
    @info "Writing to $path."
    if !isnothing(caption)
        md_path = replace(path, ".png"=>".md")
        label = "fig-" * replace(basename(path), ".png"=>"")
        @info "Writing to $md_path."
        plot_title = split(caption, ".")[1]
        kwargs = (;plot_title, kwargs...)
        caption = replace(caption, plot_title=>"**$plot_title**")
        open(md_path, "w") do fd 
            write(fd, """
            ::: {#$label}

            ![]($path)

            $caption

            :::
            """)
        end
    end
    Plots.savefig(plot!(f(); kwargs...), path)
end 
savefig(f; kwargs...) = savefig(()->f; kwargs...)
savefig(path::AbstractString; kwargs...) = f->savefig(f; path, kwargs...)
savetable(f::Function; path, force=true, caption, kwargs...) = if force || !isfile(path)
    mkpath(dirname(path))
    @info "Writing to $path."
    label = "tbl-" * replace(basename(path), ".md"=>"")
    plot_title = split(caption, ".")[1]
    caption = replace(caption, plot_title=>"**$plot_title**")
    df = f()
    open(path, "w") do fd 
        write(fd, """
        ::: {#$label}

        """)
        pretty_table(fd, df; backend=:markdown, show_first_column_label_only=true, kwargs...)
        write(fd, """

        $caption

        :::
        """)
    end
end 
savetable(f; kwargs...) = savetable(()->f; kwargs...)
savetable(path::AbstractString; kwargs...) = f->savetable(f; path, kwargs...)

scaledto(x; lo=0, hi=1) = lo .+ (hi .- lo) .* (x .- minimum(x)) ./ (maximum(x) - minimum(x))
scaledto(x, lo, hi) = scaledto(x; lo, hi)
ternary(c, t, f) = c ? t : f