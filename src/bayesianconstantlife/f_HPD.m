function probability=f_HPD(lprior,CI)
%for the prior distribution find the Highest Prosterior Density (HPD) probability for a given comfort
%interval. This is calculated as interval with the probability of interest 
%(e.g., 95% (CI=0.95)) of the distribution mass around the center of the distribution
%Inputs: lprior 2d matrix, CI 
normlprior=exp(lprior);
ordnormlprior=sort(normlprior(:),'descend');
int=cumsum(ordnormlprior);

[~,id]=min(abs(int-CI));
probability=ordnormlprior(id);
end
