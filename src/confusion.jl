count_true_positives(decisions, true_decisions) = sum(decisions .== :mine .&& true_decisions .== :mine)
count_true_negatives(decisions, true_decisions) = sum(decisions .== :abandon .&& true_decisions .== :abandon)
count_false_positives(decisions, true_decisions) = sum(decisions .== :mine .&& true_decisions .== :abandon)
count_false_negatives(decisions, true_decisions) = sum(decisions .== :abandon .&& true_decisions .== :mine)


function Base.precision(decisions::Vector, true_decisions::Vector)
    tp = count_true_positives(decisions, true_decisions)
    fp = count_false_positives(decisions, true_decisions)
    return tp / (tp+fp)
end


function recall(decisions, true_decisions)
    tp = count_true_positives(decisions, true_decisions)
    fn = count_false_negatives(decisions, true_decisions)
    return tp / (tp+fn)
end


function get_true_decisions(results::Dict; extraction_cost=150)
    return map(r_massive->r_massive > extraction_cost ? :mine : :abandon, results[:r_massive])
end


function accuracy(results::Dict; extraction_cost=150)
    decisions = results[:last_action]
    true_decisions = get_true_decisions(results; extraction_cost=extraction_cost)
    return accuracy(decisions, true_decisions)
end


function accuracy(decisions, true_decisions)
    tp = count_true_positives(decisions, true_decisions)
    tn = count_true_negatives(decisions, true_decisions)
    fp = count_false_positives(decisions, true_decisions)
    fn = count_false_negatives(decisions, true_decisions)
    return (tp+tn) / (tp+tn+fp+fn)
end

function confusion_matrix(decisions, true_decisions)
    N = length(decisions)
    tpr = count_true_positives(decisions, true_decisions) / N
    tnr = count_true_negatives(decisions, true_decisions) / N
    fpr = count_false_positives(decisions, true_decisions) / N
    fnr = count_false_negatives(decisions, true_decisions) / N
    return [tpr fpr;
            fnr tnr]
end

function plot_confusion(cm)
    heatmap(rotr90(cm)', c=:viridis, ratio=1, xlims=(0.5, 2.5), ylims=(0.5, 2.5), clims=(0,1))
end