module MixedFidelityModelSelection

import CairoMakie: Figure, Axis, density!, lines!, text!
using ColorSchemes
using Distributions
using LaTeXStrings
using LinearAlgebra
using Markdown
using Measures
using Plots; default(fontfamily="Computer Modern", framestyle=:box)
using POMDPs
using Parameters
using Random
using Statistics
using StatsBase

export
    JobParameters,
    Configuration,
    Trial,
    Results
include("parallel.jl")

export
    configurations,
    initialize,
    evaluate,
    save
include("interface.jl")

export
    count_true_positives,
    count_true_negatives,
    count_false_positives,
    count_false_negatives,
    precision,
    recall,
    accuracy,
    confusion_matrix,
    get_true_decisions,
    confusion_table,
    get_colorbar_accuracy,
    plot_confusion,
    plot_accuracies,
    plot_accuracy
include("confusion.jl")

export
    regret,
    get_colorbar_regret,
    plot_sweep_regret,
    plot_regret,
    plot_cumulative_regret
include("regret.jl")

export
    goal_programming,
    plot_pareto,
    plot_pareto_3d
include("pareto.jl")

export
    initialize_betas,
    update_betas,
    betavalue
include("betas.jl")

export
    plot_ore_mass_distributions,
    plot_relative_errors,
    plot_rel_error_aggregate,
    draw_zero,
    get_data_xy,
    get_colorbar,
    get_colorbar_bounded,
    bounded_values_multi,
    get_colorbar_bounded_multi
include("visualizations.jl")

export
    sample_massive,
    filltomatch!,
    mismatch_mean,
    mismatch_std
include("utils.jl")

export
    BeliefMDP
include("belief_mdp.jl")

end # module
