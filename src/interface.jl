# (::Type{<:Configuration})::Vector{<:Configuration}
function configurations end

# (config::Configuration)::Trial
function initialize end

# (trial::Trial)::Results
function evaluate end

# (results::Results)
function save end
