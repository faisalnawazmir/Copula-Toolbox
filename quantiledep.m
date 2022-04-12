function [out1,out2,out3] = quantiledep(X,Y,q,bootreps)
% Computes the sample estimator of the quantile dependence measure ('near tail' dependence)
%
% INPUTS:	X, a Tx1 vector of data
%			Y, a Tx1 vector of data
%			q, a kx1 vector of quantiles to evaluate the dependence measure at
%
% OUTPUT:	out1, a kx1 vector of the measure at each quantile (output is k+1 if one of the q values is exactly 0.5)
%           out2, a kxk matrix, the bootstrap covariance matrix of the estimated quantile dep measuures
%
% Andrew Patton
%
% Tuesday, 24 April, 2001
%
% REVISED 24 March 2011 to include boostrap std errors, based on theory in Remillard (2010, wp) and Oh and Patton (2011, wp)

if nargin<4 || isempty(bootreps) || bootreps==0
    bootreps=0;
    out3 = nan;
end

U = empiricalCDF(X);
V = empiricalCDF(Y);

T = size(X,1);
k = size(q,1);

if sum(q==0.5)>0
    k=k+1;
end

out1 = nan(k,1);
out2 = nan(k,k);
out3 = nan(k,bootreps);

out1(:,1) = quantiledep_calc(U,V,q);

if bootreps>0
    bootdates = randi(1,T,T);  % a set of dates to use in bootstrap (data are assumed iid so just shuffle dates without using blocks
    
    outboot = nan(bootreps,k);
    for bb=1:bootreps
        outboot(bb,:) = quantiledep_calc(U(bootdates(:,bb)),V(bootdates(:,bb)),q);
    end
    out1(:,2) = std(outboot)';
    out2 = cov(outboot);
    out3 = outboot';
end

% helper function: computes the quantile dependence
    function qdep = quantiledep_calc(U,V,q)
        
        qdep= nan(length(q)+(sum(q==0.5)>0),1);
        counter=0;
        for jj = 1:length(q)
            if q(jj)<0.5    	% then this is a 'lower' quantile dependence measure
                qdep(jj) = mean( (U<=q(jj)).*(V<=q(jj)) )/q(jj);
            elseif q(jj)==0.5  % then it's right in the middle: return both lower and upper
                counter=1;
                qdep(jj) = mean( (U<=q(jj)).*(V<=q(jj)) )/q(jj);
                qdep(jj+counter) = ( 1-2*q(jj)+mean( (U<=q(jj)).*(V<=q(jj)) ) ) /(1-q(jj));
            else
                qdep(jj+counter) = ( 1-2*q(jj)+mean( (U<=q(jj)).*(V<=q(jj)) ) ) /(1-q(jj));
            end
        end
    end


end