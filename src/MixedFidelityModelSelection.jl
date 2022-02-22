module MixedFidelityModelSelection

import CairoMakie: Figure, Axis, density!, lines!, text!
using ColorSchemes
using LaTeXStrings
using Plots; default(fontfamily="Computer Modern", framestyle=:box)
using POMDPs
using Parameters
using Random
using Statistics

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
    plot_confusion
include("confusion.jl")

export
    plot_ore_mass_distributions,
    plot_relative_errors,
    plot_rel_error_aggregate,
    draw_zero
include("visualizations.jl")

export
    sample_massive,
    filltomatch!,
    mismatch_mean,
    mismatch_std
include("utils.jl")

end # module
