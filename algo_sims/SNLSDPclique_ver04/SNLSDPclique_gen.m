function [Xorig,A,Dpartial,stats] = SNLSDPclique_gen(n,m,r,R,nf,randstate)
% SNLSDPclique_gen  Generate a random instance of the Sensor Network 
%                   Localization problem in dimension r.
%
% [Xorig,A,Dpartial] = SNLSDPclique_gen(n,m,r,R,nf,randstate) generates 
% the original sensor positions Xorig, the anchor positions A, and the 
% partial distance (squared) matrix Dpartial where entry (i,j) is 
% ||pi-pj||^2 if the distance between sensor i and senor j, ||pi-pj||, 
% is less than the radio range R; otherwise Dpartial(i,j) = 0;
% Note that the matrix [Xorig; A] has n rows and r columns.
%
% The sensors and anchors are uniformly distributed in the [0,1]^r region.
%
% Optional parameters:
%
%        nf :   the noise factor applied to the exact distances to obtain 
%               a noisy partial distance matrix.  The noise is not applied 
%               to the distances between the anchors
%
% randstate :   an integer between 0 and 2^32 which sets the default
%               random number stream
% 
% =====================
% Citation information:
% =====================
% Nathan Krislock and Henry Wolkowicz. Explicit sensor network localization
% using semidefinite representations and facial reductions. SIAM Journal on
% Optimization, 20(5):2679-2708, 2010.

% SNLSDPclique, version 0.4
% Copyright (C) 2009 by
% Nathan Krislock, Henry Wolkowicz
% Last Modified June 20, 2011



%% Process inputs and do error-checking

if (nargout > 4)
    error('SNLSDPclique_gen:TooManyOutputs', 'Too many output arguments.')
elseif (nargout < 4)
    error('SNLSDPclique_gen:TooManyOutputs', 'Too few output arguments.')
end

if (nargin < 4)
    error('SNLSDPclique_gen:TooFewInputs','Too few input arguments.');
elseif (nargin == 4)
    nf = 0;
    randstate = [];
elseif (nargin == 5)
    randstate = [];
elseif (nargin > 6)
    error('SNLSDPclique_gen:TooManyInputs','Too many input arguments.');
end


%% generate random points in dim r

if ~isempty(randstate)
    
    % For MATLAB 7.7 and later
%     RandStream.setDefaultStream(RandStream('mt19937ar','seed',randstate));
      RandStream.setGlobalStream(RandStream('mt19937ar','seed',randstate));
    
    % For MATLAB versions prior to 7.7
    %rand('twister',randstate);
    
end

Porig = rand(n,r);
scatter(Porig(1:end,1),Porig(1:end,2))
Xorig = Porig(1:n-m,:);
A = Porig(n-m+1:end,:);


%% Compute Dpartial

II = cell(n,1);
JJ = cell(n,1);
DD = cell(n,1);

% This memory saving technique for computing the partial distance matrix
% is due to Franz Rendl of the University of Klagenfurt.
for ip = 1:n-m
    I = false(n,1);
    %Pdiff = bsxfun(@minus,Porig(ip+1:end,:),Porig(ip,:));
    Pdiff = Porig(ip+1:end,:) - repmat(Porig(ip,:),n-ip,1);
    I(ip+1:end) = all(abs(Pdiff) < R,2);
    %Pdiff = bsxfun(@minus,Porig(I,:),Porig(ip,:));
    Pdiff = Porig(I,:) - repmat(Porig(ip,:),nnz(I),1);
    Pdist = sum(Pdiff.^2,2);
    small_dist = Pdist < R^2;
    Pdist = Pdist(small_dist);
    I(I) = small_dist;
    II{ip} = ip*ones(sum(I),1);
    JJ{ip} = find(I);
    DD{ip} = Pdist;
end

Ainds = n-m+1:n;
AJJ = repmat(Ainds,m,1);
AII = AJJ';
ADD = triu(K(A*A'));

II = [cell2mat(II); AII(:)];
JJ = [cell2mat(JJ); AJJ(:)];

DDnoisy = cell2mat(DD);
DDnoisy = DDnoisy.*( (1+nf*randn(length(DDnoisy),1)).^2 );
DD = [DDnoisy; ADD(:)];

Dpartial = sparse(II,JJ,DD,n,n);
Dpartial = Dpartial + Dpartial';

v = full(sum(Dpartial>0));

stats.mindeg = min(v);
stats.maxdeg = max(v);
stats.avgdeg = sum(v)/2/n;


end


%% ========= External functions ==========

%% K
function D = K(B)
% D = K(B)
%
% Operator K of B  - forms EDM if B is psd

D = zeros(size(B));
for ii=1:length(B)
    for jj=1:length(B)
        D(ii,jj) = B(ii,ii) + B(jj,jj) - 2*B(ii,jj);
    end
end

% v = diag(B);
% vt = v';
% D = bsxfun(@plus,v,vt);
% D = D - 2*B;
end
