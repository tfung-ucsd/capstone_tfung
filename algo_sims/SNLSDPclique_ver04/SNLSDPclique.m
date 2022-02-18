function [X,exitflag,OUTPUT] = SNLSDPclique(Dpartial,A,r,OPTS)
% SNLSDPclique      Solve the Sensor Network Localization problem using 
%                   the SDP facial reduction method
%
% X = SNLSDPclique(Dpartial,A,r) solves the sensor network problem with
% partial distance (squared) matrix Dpartial, and anchor positions A, in
% dimension r.  
%
% Set A=[] to solve the anchorless (Euclidean distance matrix completion) 
% problem in dimension r.
% 
% NOTE:  When A~=[], the distances between the anchors must be in the
% bottom right corner of the matrix Dpartial, and A must have r columns.
%
% X = SNLSDPclique(Dpartial,A,OPTS) solves the sensor network problem
% with the default parameters replaced by those in the structure OPTS with:
%
%   OPTS.level      =   {1,[2],3} where:
%                           1 : perform only Rigid Clique Union
%                           2 : perform 1 and Rigid Node Absorption
%                           3 : perform 2 and Non-Rigid Clique Union
%
%   NOTE:   When OPTS.level > 2, this program requires the "fminsearch"
%           function from the MATLAB Optimization Toolbox.
%
%   OPTS.verbose    =   {[0],1,2,3} where 0 is the default no output,
%                       1 prints minimal output of the progress of the
%                       algorithm, 2 prints full output, and 3 provides
%                       graphical output
%
% [X,exitflag] = SNLSDPclique(...) returns an exitflag with:
%
%   exitflag = 1 if some sensors were localized
%   exitflag = 0 if no sensors were localized
%
% [X,exitflag,OUTPUT] = SNLSDPclique(...) returns a structure OUTPUT with:
%
%   OUTPUT.cliques              =   entry i of this sparse vector contains
%                                   the size of clique i
%
%   OUTPUT.cliquematrix         =   column i of this sparse matrix is the
%                                   indicator vector of clique i
%
%   OUTPUT.cliqueintersections  =   (i,j) entry of this sparse matrix
%                                   gives the size of the intersection of
%                                   clique i and clique j
%
%   OUTPUT.cliquefaces          =   the {i} entry of this cell array
%                                   contains the matrix describing the face
%                                   of the remaining clique i
%
% =====================
% Citation information:
% =====================
% Nathan Krislock and Henry Wolkowicz. Explicit sensor network localization
% using semidefinite representations and facial reductions. SIAM Journal on
% Optimization, 20(5):2679-2708, 2010.
%
% Special thanks to Ting Kei Pong for his feedback on the code and for his
% help locating the source of a particular bug.  Due to his help, the code is 
% now much more compatible with different versions of MATLAB.

% SNLSDPclique, version 0.4
% Copyright (C) 2009 by
% Nathan Krislock, Henry Wolkowicz
% Last Modified July 4, 2011


%% Initialize parameters

condtolerscaling = 1e2;  % number of digits that are possibly lost
                         % from scaling in the subspace intersection
maxscaling = 1e8;  % maximum increase of the condition scaling


%% Process inputs and do error-checking

if (nargout > 3)
    error('SNLSDPclique:TooManyOutputs', 'Too many output arguments.')
end

if (nargin < 3)
    error('SNLSDPclique:TooFewInputs','Too few input arguments.');
elseif (nargin == 3)
    OPTS.level = 2;
    OPTS.verbose = 0;
elseif (nargin > 4)
    error('SNLSDPclique:TooManyInputs','Too many input arguments.');
end

n = length(Dpartial);
m = size(A,1);

% Check that all anchor distances exist
if ( nnz(Dpartial(end-m+1:end,end-m+1:end)) + m ) ~= m^2,
    error('SNLSDPclique:MissingAnchorDistances',...
        'Some distances for anchors are missing.');
end

if OPTS.verbose >= 1
    disp(' ');
    disp('=================== Start of SNLSDPclique ===================');
end



%% Grow cliques & Initialize
eigvs = '';
GrowCliques

intcliques = Dcq'*Dcq;
csizes = sparse(diag(intcliques));
intcliques = triu(intcliques,1);

% needed at the end and in node absorption step
Dcqinit = Dcq;
csizesinit = full(csizes);



%% Initialize indicators
Co = find(Cp)';

grow = 1;
paramchange = 0;
mainwiters = 0;          % # iterations through while loop

if OPTS.verbose >= 1
    disp('Finished initialization -- starting WHILE loop');
end


%% MAIN WHILE LOOP:
%  while number of cliques is reduced or a clique grows or a parameter
%  changes paramchange = 1 if flagred = 1 occurs. This signals that an
%  increase in condtolerscaling should be tried.

while length(Co)>1 && (grow || paramchange),
    
    %% Reset indicators
    grow = 0;
    paramchange = 0;
    mainwiters = mainwiters + 1;
    
    if OPTS.verbose >= 1
        disp('______________________________________');
        disp(['# cliques remaining is ',num2str(full(sum(Cp)))]);
    end
    
    
    %% Rigid Clique Union
    
    if OPTS.verbose >= 1
        disp('   Rigid Clique Union');
    end
    
    for ict = Co  % for each clique ic
        
        if Cp(ict) % if ic is still a clique
            
            if sum(Dcq(:,ict)) >= r+1
                
                % Find the cliques intersecting clique ic
                icintcliques = Dcq'*Dcq(:,ict);
                icintcliques(ict) = 0;
                Copj = find(icintcliques)';
                
                for jc = Copj
                    
                    flagred = RigidCliqueUnion(ict,jc);
                    
                    if flagred == 1,
                        paramchange = 1; % try lower precision
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
    %% Rigid Node Absorption
    
    if ~grow && OPTS.level >= 2
        
        if OPTS.verbose >= 1
            disp('      Rigid Node Absorption');
            
        end
        
        for ict = Co
            
            p = Dcq(:,ict)>0;
            
            if sum(p) >= r+1  % clique ic >= r+1,
                
                % Find nodes connected to at least r+1 nodes in clique ic
                NeighbourDegrees = sum(Dpartial(:,p)>0,2) .* ~p;
                [sND,inds] = sort(full(NeighbourDegrees),'descend');
                icConnectedNodes = inds(sND >= r+1)';
                
                numnodes = min(length(icConnectedNodes),r+1);
                
                for jc = icConnectedNodes(1:numnodes)
                    
                    flagred = RigidNodeAbsorption(ict,jc);
                    
                    if flagred == 1,
                        paramchange = 1; % try lower precision
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
    %% Non-Rigid Clique Union
    
    if ~grow && OPTS.level >= 3
        
        if OPTS.verbose >= 1
            disp('         Non-Rigid Clique Union');
        end
        
        for ict = Co  % for each clique ic
            
            if Cp(ict) % if ic is still a clique
                
                p = Dcq(:,ict)>0;
                
                if nnz(p) >= r+1
                    
                    % Find the cliques intersecting clique ic
                    icintcliques = Dcq'*p;
                    icintcliques(ict) = 0;
                    Copj = find(icintcliques==r)';
                    
                    for jc = Copj
                        
                        flagred = NonRigidCliqueUnion(ict,jc);
                        
                        if flagred == 1,
                            paramchange = 1; % try lower precision
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
    
    %% Reset everything
    
    Cp  = logical(spones(Cp));  % save memory
    Dcq = double(spones(Dcq));
    
    % Compute sizes of clique intersections and clique sizes
    intcliques = Dcq'*Dcq;
    csizes = sparse(diag(intcliques));
    intcliques = triu(intcliques,1);
    
    Co = find(Cp)';
    
    if OPTS.verbose >= 3
        figure(2)
        spy(Dcq);
        pause(.000001);
    end
    
    
    %% Scaling
    if ~grow && paramchange
        if  condtolerscaling < maxscaling
            condtolerscaling = 10*condtolerscaling;
            if OPTS.verbose >= 1
                disp(['condtolerscaling increased to ',...
                    num2str(condtolerscaling)]);
            end
        else
            paramchange = 0;
        end
    end
    
    
end



%% Print Output

if OPTS.verbose >= 1
    disp('______________________________________');
    disp('End of WHILE loop');
    disp(['number of remaining cliques is: ', ...
        num2str(full(sum(Cp)))])
    disp(['sizes of remaining cliques is: ', ...
        num2str(full(csizes(csizes>0)'))]);
end
clear NeighbourDegrees e22 e33 e11 p a1 a2


%% Compute sensor positions


if isempty(A)
    
    % let ic be the clique of maximum size
    [dummy,ic] = max(csizes);
    
    [P,flagred] = CompleteClique(ic);
    
    if flagred
        X = [];
    else
        p = Dcq(:,ic) > 0;
        X = P(p,:);
    end
    
else
    
    % let ic be the clique of maximum size containing the anchors
    anchors = sparse(false(n,1));  anchors(n-m+1:end) = true;
    acliques = (anchors'*Dcq == m)';
    [dummy,ic] = max(csizes.*acliques);
    
    if csizes(ic) >= r+1
        [P,flagred] = CompleteClique(ic);
    else
        flagred = 1;
    end
    
    if flagred
        if OPTS.verbose >= 1
            disp('completeclique FLAGRED during final completion')
        end
        X = [];
    else
        
        icsensors = Dcq(:,ic) > anchors;
        
        if nnz(icsensors) %all(ismember(n-m+1:n,inds))
            % if clique ic contains more than just anchors
            
            if OPTS.verbose >= 2
                disp(['Final completion:  ',...
                    'largest clique containing anchors ',...
                    'also contains sensors']);
            end
            
            % Translate P so that P2 is centred at zero
            P2 = P(n-m+1:n,:);
            mp2 = sum(P2)/m;  % centre of P2
            P = P - repmat(mp2,n,1);  % translate P
            
            % Translate A so that A is centred at zero
            ma = sum(A)/m;  % centre of A
            A = A - repmat(ma,m,1);  % centre A at zero
            
            % solve min ||A-P2*Q|| s.t. Q'*Q = I
            % this is the orthogonal Procrustes problem
            %   -- Golub, Van Loan, Algorithm 12.4.1
            C = P2'*A;
            [Uc,dummy,Vc] = svd(C);
            Q = Uc*Vc';
            
            % Compute final X
            P = P*Q;  % rotate P to align P2 with anchors A
            
            P = P + repmat(ma,n,1);  % translate P back
            A = A + repmat(ma,m,1);  % translate A back
            
            P2 = P(n-m+1:n,:);
            
            if OPTS.verbose >= 1
                disp(['Anchor agreement error = ',...
                    num2str(norm(A-P2,'fro'))]);
            end
            
            X = P(icsensors,:);  % return X = P1
            
        else
            
            X = [];
            
        end
        
    end
    
end


%% Set exitflag and OUTPUT

if isempty(X)
    exitflag = 0;
else
    exitflag = 1;
end

OUTPUT.cliques              =   csizes;
OUTPUT.cliquematrix         =   Dcq;
OUTPUT.cliqueintersections  =   intcliques;
OUTPUT.cliquefaces          =   eigvs;

if OPTS.verbose >= 1
    disp('=================== End of SNLSDPclique ===================');
    disp(' ');
end



%% ========= Nested functions ==========

%% GrowCliques
    function GrowCliques
        
        if OPTS.verbose >= 1
            disp('Try to grow cliques');
        end
        
        MaxCliqueSize = 3*(r+1);
        
        tt = tic;
        II = zeros(n*MaxCliqueSize+m,1);
        JJ = zeros(n*MaxCliqueSize+m,1);
        k = 0;
        for jct = 1:n % For each clique jc
            
            % Add node jc to clique jc
            k = k + 1;  II(k) = jct;  JJ(k) = jct;
            NodesToAdd = MaxCliqueSize - 1;
            
            q = (Dpartial(:,jct)>0);
            while NodesToAdd && nnz(q)
                ic = find(q);
                % Add node ic to clique jc
                k = k + 1;  II(k) = ic(end);  JJ(k) = jct;
                NodesToAdd = NodesToAdd - 1;
                q = q & Dpartial(:,ic(end));
            end
            
        end
        tt = toc(tt);
        
        if OPTS.verbose >= 1
            disp(['Time to grow cliques is  ', num2str(tt)]);
        end
        
        
        if isempty(A)
            Dcq = sparse(II(1:k),JJ(1:k),1,n,n);  clear ii jj
            Cp = sparse(1:n,1,true);
            
            eigvs = cell(n,1);  % Cell array used to store eigenvectors
            %for iq = 1:n
            %    eigvs{iq} = sparse([],[],[],n,r);
            %end
        else
            % Add anchor clique to column n+1
            II(k+1:k+m) = n-m+1:n;
            JJ(k+1:k+m) = n+1;
            k = k + m;
            
            Dcq = sparse(II(1:k),JJ(1:k),1,n,n+1);  clear II JJ
            Cp = sparse(1:n+1,1,true);
            
            eigvs = cell(n+1,1);  % Cell array used to store eigenvectors
            %for iq = 1:n+1
            %    eigvs{iq} = sparse([],[],[],n,r);
            %end
        end
                        
        if OPTS.verbose >= 3
            figure(2)
            spy(Dcq);
            pause(.000001);
        end
        
    end


%% RigidCliqueUnion
    function flagred = RigidCliqueUnion(ic,jc)
        
        flagred = 0;
        
        e22 = Dcq(:,ic) & Dcq(:,jc);  % intersect cliques ic and jc
        
        e11 = (Dcq(:,ic) - e22)>0;   e33 = (Dcq(:,jc) - e22)>0;
        
        inds = [find(e11); find(e22); find(e33)];
        
        if ~any(e11)   % nothing in e11
            % clique ic is a subset of clique jc,
            % so copy clique jc onto clique ic and delete clique jc
            
            Cp(jc) = 0;
            Dcq(:,ic) = Dcq(:,jc);
            Dcq(:,jc) = 0;
            eigvs{ic} = eigvs{jc};
            eigvs{jc} = [];
            grow = 1;
            
        elseif ~any(e33)  % nothing in e33
            % clique jc is a subset of clique ic, so delete clique jc
            
            Cp(jc) = 0;
            Dcq(:,jc) = 0;
            eigvs{jc} = [];
            grow = 1;
            
        elseif sum(e22) >= r+1
            
            nvec = full([sum(e11) sum(e22) sum(e33)]);
            
            a1 = 1:sum(nvec(1:2));   a2 = nvec(1)+1:sum(nvec);
            a1inds = inds(a1);       a2inds = inds(a2);
            
            % Find Ub1
            if isempty(eigvs{ic})
                Dbar = full(Dpartial(a1inds,a1inds));
                B = Kdag(Dbar);
                [Ub,dummy] = eig(B); 
                Ub = Ub(:,end-r+1:end);
                k = length(a1);  e = ones(k,1);
                Ub1 = [ Ub, e/sqrt(k) ];
            else
                Ub1 = full(eigvs{ic}(a1inds,:));
            end
            
            % Find Ub2
            if isempty(eigvs{jc})
                Dbar = full(Dpartial(a2inds,a2inds));
                B = Kdag(Dbar);
                [Ub,dummy] = eig(B);
                Ub = Ub(:,end-r+1:end);
                k = length(a2);  e = ones(k,1);
                Ub2 = [ Ub, e/sqrt(k) ];
            else
                Ub2 = full(eigvs{jc}(a2inds,:));
            end
            
            % Find U
            [U,flagred] = SubspaceIntersection(nvec,Ub1,Ub2);
            
            if ~flagred
                
                % Store U
                ii = repmat(inds,1,r+1);
                jj = repmat(1:r+1,length(inds),1);
                eigvs{ic} = sparse(ii(:),jj(:),U(:),n,r+1);
                
                % Update Dcq
                Cp(jc) = 0;
                Dcq(:,ic) = Dcq(:,ic) | Dcq(:,jc);
                Dcq(:,jc) = 0;
                eigvs{jc} = [];
                grow = 1;
                
                if OPTS.verbose >= 2
                    disp(['   [ic jc |e22|] = [',...
                        num2str([ic jc nvec(2)]),']']);
                end
                
                if OPTS.verbose >= 3
                    Dcq = double(spones(Dcq));
                    figure(2);
                    spy(Dcq);
                    pause(.000001);
                end
                
            end
            
        end
        
    end


%% RigidNodeAbsorption
    function flagred = RigidNodeAbsorption(ic,jc)
        
        flagred = 0;
        
        e22 = Dpartial(:,jc) & Dcq(:,ic);  % e22 = nhbrs of jc in clique ic
        ne22 = sum(e22);
        
        if ne22 == sum(Dcq(:,ic)) && isempty(eigvs{ic})
            % node jc is connected to every node in clique ic
            % and there are no eigvs for clique ic,
            % so clique ic is an original clique --
            % simply add jc to clique ic
            
            % Update Dcq
            Dcq(jc,ic) = 1;  % absorb jc into clique ic
            grow = 1;
            
            % Output
            if OPTS.verbose >= 2
                disp(['      [ic jc |e22|] = [',...
                    num2str([ic jc sum(e22)]),'] -- simple add']);
            end
        else
            
            IntersectionComplete = (nnz(Dpartial(e22,e22))==ne22*(ne22-1));
            
            % Complete clique ic, if necessary
            if ~IntersectionComplete
                
                % Complete clique ic
                [P,flagred] = CompleteClique(ic);
                
            end
            
            % If complete clique was successful, perform node absorption
            if flagred
                if OPTS.verbose >= 2
                    disp('      completeclique FLAGRED');
                end
            else
                
                e11 = (Dcq(:,ic) - e22)>0;
                e33 = sparse(jc,1,true,n,1);
                
                inds = [find(e11); find(e22); find(e33)];
                nvec = full([sum(e11) sum(e22) sum(e33)]);
                
                a1 = 1:sum(nvec(1:2));   a2 = nvec(1)+1:sum(nvec);
                a1inds = inds(a1);       a2inds = inds(a2);
                
                % Find Ub1
                if isempty(eigvs{ic})
                    if IntersectionComplete
                        B = Kdag(full(Dpartial(a1inds,a1inds)));
                    else
                        B = Kdag(K(full(P(a1inds,:)*P(a1inds,:)')));
                    end
                    [Ub,dummy] = eig(B);
                    Ub = Ub(:,end-r+1:end);
                    k = length(a1);  e = ones(k,1);
                    Ub1 = [ Ub, e/sqrt(k) ];
                else
                    Ub1 = full(eigvs{ic}(a1inds,:));
                end
                
                % Find Ub2
                if IntersectionComplete
                    B = Kdag(full(Dpartial(a2inds,a2inds)));
                else
                    v = full(Dpartial(e22,jc));
                    B = Kdag([K(full(P(e22,:)*P(e22,:)')) v; v' 0]);
                end
                [Ub,dummy] = eig(B);
                Ub = Ub(:,end-r+1:end);
                k = length(a2);  e = ones(k,1);
                Ub2 = [ Ub, e/sqrt(k) ];
                
                % Find U
                [U,flagred] = SubspaceIntersection(nvec,Ub1,Ub2);
                
                if ~flagred
                    
                    % Store U
                    II = repmat(inds,1,r+1);
                    JJ = repmat(1:r+1,length(inds),1);
                    eigvs{ic} = sparse(II(:),JJ(:),U(:),n,r+1);
                    
                    % Update Dcq
                    Dcq(jc,ic) = 1;  % absorb jc into clique ic
                    grow = 1;
                    
                    % Output
                    if OPTS.verbose >= 2
                        disp(['      [ic jc |e22|] = [',...
                            num2str([ic jc ne22]),']']);
                    end
                    
                end
                
            end
            
        end
        
        if OPTS.verbose >= 3
            Dcq = double(spones(Dcq));
            figure(2);
            spy(Dcq);
            pause(.000001);
        end
        
    end


%% NonRigidCliqueUnion
    function flagred = NonRigidCliqueUnion(ic,jc)
        
        flagred = 0;
        
        e22 = Dcq(:,ic) & Dcq(:,jc);  % intersect cliques ic and jc
        
        e11 = (Dcq(:,ic) - e22)>0;  e33 = (Dcq(:,jc) - e22)>0;
        
        inds = [find(e11); find(e22); find(e33)];
        
        if ~any(e11)   % nothing in e11
            % clique ic is a subset of clique jc,
            % so copy clique jc onto clique ic and delete clique jc
            
            Cp(jc) = 0;
            Dcq(:,ic) = Dcq(:,jc);
            Dcq(:,jc) = 0;
            eigvs{ic} = eigvs{jc};
            eigvs{jc} = [];
            grow = 1;
            
        elseif ~any(e33)  % nothing in e33
            % clique jc is a subset of clique ic, so delete clique jc
            
            Cp(jc) = 0;
            Dcq(:,jc) = 0;
            eigvs{jc} = [];
            grow = 1;
            
        elseif sum(e22) == r
            
            nvec = full([sum(e11) sum(e22) sum(e33)]);
            
            a1 = 1:sum(nvec(1:2));  a2 = nvec(1)+1:sum(nvec);
            a1inds = inds(a1);      a2inds = inds(a2);
            k1 = length(a1);        k2 = length(a2);
            
            nota1inds = inds(sum(nvec(1:2))+1:end);
            nota2inds = inds(1:nvec(1));
            
            [II,JJ] = find(Dpartial(nota1inds,nota2inds));
            II = nota1inds(II);  JJ = nota2inds(JJ);
            
            if isempty(II)  % no other distances to compare with
                flagred = 1;
            else
                % Find Ub1
                if isempty(eigvs{ic})
                    Dbar = full(Dpartial(a1inds,a1inds));
                    B = Kdag(Dbar);
                    [Ub,dummy] = eig(B); %#ok<SETNU>
                    Ub = Ub(:,end-r+1:end);
                    e = ones(k1,1);
                    Ub1 = [ Ub, e/sqrt(k1) ];
                else
                    Ub1 = full(eigvs{ic}(a1inds,:));
                end
                
                % Find Ub2
                if isempty(eigvs{jc})
                    Dbar = full(Dpartial(a2inds,a2inds));
                    B = Kdag(Dbar);
                    [Ub,dummy] = eig(B);
                    Ub = Ub(:,end-r+1:end);
                    e = ones(k2,1);
                    Ub2 = [ Ub, e/sqrt(k2) ];
                else
                    Ub2 = full(eigvs{jc}(a2inds,:));
                end
                
                % Find U
                [U,flagred] = SubspaceIntersection(nvec,Ub1,Ub2);
                
                if ~flagred
                    % Find V such that V'*U'*e = 0
                    UTe = sum(U)';
                    [V,Rtemp] = qr(UTe);
                    if norm(V*Rtemp-UTe) > .000001
                        disp('||V*Rtemp-U^T*e|| > .000001');
                        keyboard
                    end
                    V = V(:,2:end);  UV = U*V;
                    
                    % Form A1 and A2
                    Ja1 = eye(k1) - ones(k1)/k1;
                    Ua1 = U(a1,:);  A1 = Ja1*Ua1*V;
                    
                    [U1,S1,V1] = svd(A1,0);
                    S1 = S1(1:end-1,1:end-1);
                    n1 = V1(:,end);
                    U1 = U1(:,1:end-1);
                    V1 = V1(:,1:end-1);
                    
                    Ja2 = eye(k2) - ones(k2)/k2;
                    Ua2 = U(a2,:);  A2 = Ja2*Ua2*V;
                    
                    [U2,S2,V2] = svd(A2,0);
                    S2 = S2(1:end-1,1:end-1);
                    n2 = V2(:,end);
                    U2 = U2(:,1:end-1);  V2 = V2(:,1:end-1);
                    
                    % Complete cliques ic and jc
                    [P1,flagredi] = CompleteClique(ic);
                    [P2,flagredj] = CompleteClique(jc);
                    
                    if flagredi || flagredj
                        if OPTS.verbose >= 1
                            disp('         flagredi || flagredj')
                        end
                        flagred = 1;
                    else
                        % Compute dZ
                        dZ = n1*n2'+n2*n1';  % Ai*dZ*Ai'=0, for i=1,2
                        
                        % Remove the nullspace from CC
                        [Qc,dummy] = qr(svec(dZ));
                        Q = Qc(:,2:end);
                        CC = [skron(V1',V1'); skron(V2',V2')];
                        CC = CC*Q;
                        
                        % Compute dd
                        R1 = S1\(U1'*P1(a1inds,:));
                        R2 = S2\(U2'*P2(a2inds,:));
                        dd = [svec(R1*R1'); svec(R2*R2')];
                        
                        % Find a particular solution Z
                        z = CC\dd;
                        Z = sMat(Q*z);
                        
                        % Find a pos definite particular solution Zbar
                        t = fminsearch(@(t) -min(eig(Z+t*dZ)), 0);
                        Zbar = Z + t*dZ;
                        
                        % Find optimal solution Dstar
                        d = eig(-dZ,Zbar,'chol');
                        indst = find(abs(d) > .000001);
                        ConstraintSatisfied = zeros(size(d));
                        Pstar = [];
                        
                        for ii = indst'
                            Zi = Zbar + 1/d(ii)*dZ;
                            
                            [Uz,Dz] = eig(Zi);
                            Uz = Uz(:,end-r+1:end);
                            Dz = Dz(end-r+1:end,end-r+1:end);
                            PZi = Uz*sqrt(Dz);
                            
                            Pi = UV*PZi;

                            err = zeros(length(II),1);
                            
                            for iq = 1:length(II)
                                pi = inds==II(iq);
                                pj = inds==JJ(iq);
                                err(iq) = sum((Pi(pi,:)-Pi(pj,:)).^2) - ...
                                    Dpartial(II(iq),JJ(iq));
                            end
                            
                            if norm(err) < 1e-14
                                ConstraintSatisfied(ii) = 1;
                                Pstar = Pi;
                            else
                                ConstraintSatisfied(ii) = 0;
                            end
                        end
                        
                        if sum(ConstraintSatisfied) ~= 1
                            if OPTS.verbose >= 2
                                disp(['         [ic jc |e22|] = [',...
                                    num2str([ic jc nvec(2)]),']  ',...
                                    'Constraints not sat. uniquely.']);
                            end
                            flagred = 1;
                        else
                            % Compute U
                            [Up,dummy1,dummy2] = svd(Pstar,0); %#ok<ASGLU,NASGU>
                            Up = Up(:,1:r);
                            k = sum(nvec);  e = ones(k,1);
                            U = [ Up, e/sqrt(k) ];
                            
                            % Store U
                            ii = repmat(inds,1,r+1);
                            jj = repmat(1:r+1,length(inds),1);
                            eigvs{ic} = sparse(ii(:),jj(:),U(:),n,r+1);
                            
                            % Update Dcq
                            Cp(jc) = 0;
                            Dcq(:,ic) = Dcq(:,ic) | Dcq(:,jc);
                            Dcq(:,jc) = 0;
                            eigvs{jc} = [];
                            grow = 1;
                            
                            if OPTS.verbose >= 2
                                disp(['         [ic jc |e22|] = [',...
                                    num2str([ic jc nvec(2)]),']']);
                            end
                            
                            if OPTS.verbose >= 3
                                Dcq = double(spones(Dcq));
                                figure(2);
                                spy(Dcq);
                                pause(.000001);
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end


%% SubspaceIntersection
    function [U,flagred] = SubspaceIntersection(nvec,Ub1,Ub2)
        
        flagred = 0;
        
        kb1 = nvec(1);
        kb2 = nvec(2);
        
        U1p  = Ub1(1:kb1,:);
        U1pp = Ub1(kb1+1:end,:);
        
        U2pp = Ub2(1:kb2,:);
        U2p  = Ub2(kb2+1:end,:);
                    
        cond1 = cond(U1pp);
        cond2 = cond(U2pp);
        
        if min(cond1,cond2) > condtolerscaling*kb2,
            flagred = 1;
        else
            if cond1 < cond2,
                Ubar = [U1p*(U1pp\U2pp); Ub2];
                if kb2 >= r+1
                    U = Ubar;
                elseif kb2 == r
                    u1 = null(U1pp);
                    U = [Ubar [U1p*u1; zeros(size(Ub2,1),1)] ];
                else
                    flagred = 1;
                end
            else
                Ubar = [Ub1; U2p*(U2pp\U1pp)];
                if kb2 >= r+1
                    U = Ubar;
                elseif kb2 == r
                    u2 = null(U2pp);
                    U = [Ubar [zeros(size(Ub1,1),1); U2p*u2] ];
                else
                    flagred = 1;
                end
            end
        end
        
        if flagred,
            if OPTS.verbose >= 2
                disp('FLAGRED in SubspaceIntersection: large cond number');
            end
            U = [];
        end
        
    end


%% CompleteClique
    function [P,flagred] = CompleteClique(ic)
        % Given a clique ic, compute all the sensor positions P;
        % the clique has the eigenvectors or all the distances specified
        
        % Complete the EDM of the clique ic
        icinds = find(Dcq(:,ic));
        
        if length(icinds) <= r
            
            disp('length(inds) <= r in completeclique');
            
            P = [];
            flagred = 1;
            return
            
        else
            
            if isempty(eigvs{ic})
                % clique ic must be a clique in the original graph
                
                % Check that we have a complete EDM
                p = Dcq(:,ic)>0;
                np = sum(p);
                if nnz(Dpartial(p,p)) ~= np*(np-1)
                    disp('EDM not complete in completecliques');
                    keyboard
                end
                
                D = full(Dpartial(icinds,icinds));
                B = Kdag(D);
                [Ubs,Dbs] = eig(B);
                Ubs = Ubs(:,end-r+1:end);
                Dbs = Dbs(end-r+1:end,end-r+1:end);
                
                k = length(icinds);  e = ones(k,1);
                Ub = [ Ubs, e/sqrt(k) ];
                
                II = repmat(icinds,1,r+1);
                JJ = repmat(1:r+1,length(icinds),1);
                eigvs{ic} = sparse(II(:),JJ(:),Ub(:),n,r+1);
                
                % Compute P
                P = Ubs*sqrt(Dbs);
                II = repmat(icinds,1,r);
                JJ = repmat(1:r,length(icinds),1);
                P = sparse(II(:),JJ(:),P(:),n,r);
                
                flagred = 0;
                return
                
            else
                
                % Find a good original clique inside clique ic
                if isempty(A)
                    indsPlusAnchorClique = icinds;
                else
                    indsPlusAnchorClique = [icinds; n+1];
                end
                [temps,tempi] = sort(...
                    csizesinit(indsPlusAnchorClique),'descend');
                inds3 = sum( temps >= r+1 );  % # of orig cliques >= r+1
                
                % cliques with sizes >= r+1
                cliqueinds = indsPlusAnchorClique(tempi(1:inds3))';
                
                for im = cliqueinds
                    % for orig cliques im -- sorted, size >= r+1
                    
                    % determine if orig clique im is inside clique ic
                    e22 = Dcqinit(:,im) & Dcq(:,ic);
                    
                    if sum(e22) == csizesinit(im)
                        % clique im is inside clique ic
                        
                        b = find(Dcqinit(:,im));
                        mm = length(b);
                        
                        % Find U and V
                        U = eigvs{ic};
                        
                        UTe = sum(U)';
                        [V,Rtemp] = qr(UTe);
                        if norm(V*Rtemp-UTe) > .000001
                            disp('||V*Rtemp-U^T*e|| > .000001');
                            keyboard
                        end
                        V = V(:,2:end);  UV = U*V;
                        
                        % Compute AA and B
                        Jb = eye(mm) - ones(mm)/mm;
                        Ub = U(b,:);
                        AA = Jb*Ub*V;
                        
                        B  = Kdag(full(Dpartial(b,b)));
                        [Vb,Db] = eig(B);
                        Vb = Vb(:,end-r+1:end);
                        Db = Db(end-r+1:end,end-r+1:end);
                        
                        condB = cond(Db);
                        
                        if condB > 50
                            if OPTS.verbose >= 2
                                disp(['****** completeclique condB = ',...
                                    num2str(condB)]);
                            end
                            P = [];
                            flagred = 1;
                            return
                        else
                            % Compute P
                            Zb = linsolve(AA,Vb*sqrt(Db));
                            P = UV(icinds,:)*Zb;
                            
                            II = repmat(icinds,1,r);
                            JJ = repmat(1:r,length(icinds),1);
                            P = sparse(II(:),JJ(:),P(:),n,r);
                            
                            flagred = 0;
                            return
                        end
                        
                    end
                    
                end
                
                P = [];
                flagred = 1;
                return
                
            end
            
        end
        
    end

end


%% =============== External functions ================

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


%% Kdag
function B = Kdag(D)
% B = Kdag(D)
%
% Pseudoinverse of K - forms centred Gram matrix B if D is an EDM

nn = length(D);

D = 0.5*D;
vn = sum(D,2)/nn;
evnn = sum(vn)/nn;

%B = bsxfun(@plus,vn,vn');
B = zeros(nn);
for ii=1:nn
    for jj=1:nn
        B(ii,jj) = vn(ii) + vn(jj);
    end
end

B = B - D;
B = B - evnn;
end


%% sMat
function X = sMat(x)
% X = sMat(x)

tn=length(x);
n=floor(sqrt(2*tn));
E=true(n);
indsu= triu(E,0);
indsuu=triu(E,1);
X=zeros(n);
X(indsu)=x;
v=diag(X);
X(1:n+1:end)=zeros(n,1);
X(indsuu)=X(indsuu)/sqrt(2);
X=X+X';
X=X+diag(v);

end


%% svec
function x = svec(X)
% x = svec(X)

Utri = triu(true(size(X)));
Ustri = triu(Utri,1);

x = X(:);
x(Ustri) = sqrt(2) * x(Ustri);
x = x(Utri);

end


%% skron (Copyright J.V. Burke, A.S. Lewis, M.L. Overton, 2001.)
function S = skron(A,B)
% symmetrized Kronecker product
%
% skron is a matrix representation for the linear operator from
% S_n to S_n:
%        K    ->    1/2 (BKA' + AKB')
% so that
%      skron(A,B) svec(K) = 1/2 svec(BKA' + AKB').
% where svec maps a symmetric matrix to a vector of length
% n(n+1)/2, multiplying off diagonal components by sqrt(2),
% so svec(K)'svec(L) = tr KL.
% This last property implies that svec is a unitary operator.
% kron is a matrix representation for the linear operator from
% M_n to M_n:
%        K    ->    BKA'
% so that
%      kron(A,B) vec(K) = vec(BKA')
% and therefore
%      (kron(A,B) + kron(B,A)) vec(K) = vec(BKA' + AKB').
% where vec maps a square matrix to a vector of length n^2.
% We make use of the identity
%
%      skron(A,B) = svec mat 1/2 (kron(A,B) + kron(B,A)) vec smat
%                 = svec mat 1/2 (kron(A,B)
%                                   + kron(B,A)) mat^* svec^*
%                 = svec sym mat 1/2 (kron(A,B)
%                                   + kron(B,A)) mat^* sym^* svec^*
% where mat is the inverse of vec (also the adjoint of vec) and
%      smat is the inverse of svec (also the adjoint of svec).
%       sym symmetrizes a matrix
%                   (sym^* embeds a symmetric matrix in M_n)
%
% Thus to compute skron(A,B), we start with
% 1/2 (kron(A,B) + kron(B,A)) (which is provided by the Matlab
% built-in function kron), then apply the operator "svec sym mat"
% to each of its columns, and then apply the same operator to the
% transposes of each of the resulting rows, and transpose the
% result.
% Copyright J.V. Burke, A.S. Lewis, M.L. Overton, 2001.

% Modified slightly by N. Krislock, 2010.

[m,n] = size(A);
K = 0.5*(kron(A,B) + kron(B,A));
L = zeros(m*(m+1)/2,n^2);  % Added by N. Krislock
for k = 1: n^2
    col = K(:,k);
    M = reshape(col,m,m);
    M = 0.5*(M + M');
    L(:,k) = svec(M);
    %L(:,k) = svec(M,n);
end
S = zeros(m*(m+1)/2,n*(n+1)/2);  % Added by N. Krislock
for k = 1: m*(m+1)/2
    row = L(k,:);
    M = reshape(row',n,n);
    M = 0.5*(M + M');
    %S(k,:) = (svec(M,n))';
    S(k,:) = (svec(M))';
end

end
