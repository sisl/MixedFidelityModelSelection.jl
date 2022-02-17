function configurations(enumerations)::Configuration end
function initialize(config::Configuration)::Trial end
function evaluate(trial::Trial)::Results end
function save(results::Results) end
