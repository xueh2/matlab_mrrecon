
function priors = normalizePriors(priors, notincluded)
% normalize the probability



priors(:) = priors(:) ./ sum(priors);
return