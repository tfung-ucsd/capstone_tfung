%% Runfile for SNLSDPclique
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
% Last Modified July 1, 2011

disp(' ');
disp('SNLSDPclique, version 0.3');
disp('Copyright (C) 2009 by Nathan Krislock and Henry Wolkowicz');
disp('This program comes with ABSOLUTELY NO WARRANTY. This is free');
disp('software, and you are welcome to redistribute it under certain');
disp('conditions; see Copyright_SNLSDPclique.txt and glp.txt for');
disp('details.');

disp(' ');
disp('Enter values for n (# of sensors), r (dim.), m (# of anchors),');
disp('                 R (radio range), and nf (noise factor):');
disp('(e.g., n = 10000, r = 2, m = 4, R = .04, nf = 0):');
% n  = input('n  = ');
% r  = input('r  = ');
% m  = input('m  = ');
% R  = input('R  = ');
% nf = input('nf = ');

n  = 2;
r  = 2;
m  = 2;
R  = 0.5;
nf = 0;

% Select the number of random problems you wish to solve
NumProbs = 5;

disp(' ');
% fprintf('Enter values for OPTS.level {1,2,3} ');
% fprintf('and OPTS.verbose {0,1,2,3}:\n');
% OPTS.level = input('OPTS.level = ');
% OPTS.verbose = input('OPTS.verbose = ');

% see 'help SNLSDPclique' for a description
%      of the OPTS parameter
OPTS.level = 30;
OPTS.verbose = 3;

NN = n;
RR = R;
NF = nf;

NN = NN + m;

RandStates = 1:NumProbs;

t = clock;

filename = ['results_',...
    int2str(t(1)),'y-',...
    int2str(t(2)),'m-',...
    int2str(t(3)),'d_',...
    int2str(t(4)),'h-',...
    int2str(t(5)),'m-',...
    int2str(t(6)),'s.txt'];

for n = NN
    
    for R = RR
        
        for nf = NF
            
            AvgAvgDeg = 0;
            AvgPositioned = 0;
            AvgCPU = 0;
            AvgMaxError = 0;
            AvgRMSD = 0;
            
            NumSuccess = 0;
            
            for prob = 1:NumProbs
                
                fprintf('\nRunning problem = %d / %d\n',prob,NumProbs);
                
                randstate = RandStates(prob);
                
                tt = tic;
                disp(' ');
                fprintf('Running [Xorig,A,Dpartial,stats] = ');
                fprintf('SNLSDPclique_gen(n,m,r,R,nf,randstate)...\n');
                %profile on
                [Xorig,A,Dpartial,stats] = ...
                    SNLSDPclique_gen(n,m,r,R,nf,randstate);
                %profile report; profile off
                tt = toc(tt);
                
                fprintf('min degree = %g\n',stats.mindeg);
                fprintf('max degree = %g\n',stats.maxdeg);
                fprintf('avg degree = %.1f\n',stats.avgdeg);
                fprintf('Elapsed time is %.1f seconds.\n',tt);
                
                tt = tic;
                disp(' ');
                fprintf('Running [X,exitflag,OUTPUT] = ');
                fprintf('SNLSDPclique(Dpartial,A,r,OPTS)...\n');
                %profile on
                [X,exitflag,OUTPUT] = SNLSDPclique(Dpartial,A,r,OPTS);
                %profile report; profile off
                tt = toc(tt);
                
                fprintf('Elapsed time is %.1f seconds.\n',tt);
                
                if exitflag == 1                    
                    Dcq = OUTPUT.cliquematrix;
                    csizes = OUTPUT.cliques;
            
                    % let ic be the clique of maximum size 
                    % containing the anchors
                    anchors = sparse(false(n,1));  
                    anchors(n-m+1:end) = true;
                    acliques = (anchors'*Dcq == m)';
                    [dummy,ic] = max(csizes.*acliques);
                    inds = find(Dcq(:,ic));
                    inds = setdiff(inds,n-m+1:n);
                    
                    disp(' ');
                    disp('==============================================');
                    disp([num2str(length(inds)),'/',num2str(n-m),...
                        ' sensors localized']);
                    
                    XmXorig        = X - Xorig(inds,:);
                    norm2XmXorig   = sum(XmXorig.^2,2);
                    max_sensor_err = sqrt(max(norm2XmXorig));
                    rmsd_error     = ...
                        sqrt(sum(norm2XmXorig)/length(inds));
                    
                    disp(' ');
                    fprintf('max_sensor_err = %.1d\n', max_sensor_err);
                    fprintf('rmsd_error     = %.1d\n', rmsd_error);
                    
                    AvgPositioned = AvgPositioned + length(inds);
                    AvgMaxError   = AvgMaxError   + max_sensor_err;
                    AvgRMSD       = AvgRMSD       + rmsd_error;
                    
                    NumSuccess = NumSuccess + 1;
                else
                    disp(' ');
                    disp(['0/',num2str(n-m),' sensors localized']);
                    disp(' -- try increasing OPTS.level');
                end
                
                AvgCPU    = AvgCPU    + tt;
                AvgAvgDeg = AvgAvgDeg + stats.avgdeg;
                
                disp(' ');
                fprintf('r  = %d\n',r);
                fprintf('R  = %g\n',R);
                fprintf('nf = %.0e\n',nf);
                disp('==============================================');
                fprintf('NumSuccess = %d / %d\n\n', NumSuccess,prob);
                %fprintf('\nPress any key.\n');  pause
                
            end
            
            AvgAvgDeg     = AvgAvgDeg/NumProbs;
            AvgPositioned = AvgPositioned/NumProbs;
            AvgCPU        = AvgCPU/NumProbs;
            
            AvgMaxError = AvgMaxError/NumSuccess;
            AvgRMSD     = AvgRMSD/NumSuccess;
            
            diary(filename);
            
            fprintf('n  = %g\n',n-m);
            fprintf('r  = %g\n',r);
            fprintf('m  = %g\n',m);
            fprintf('R  = %g\n',R);
            fprintf('nf = %.0e\n\n',nf);
            
            fprintf('OPTS.level    = %g\n',OPTS.level);
            fprintf('NumSuccess    = %d / %d\n\n', NumSuccess,NumProbs);
            
            fprintf('AvgAvgDeg         = %.1f\n', AvgAvgDeg);
            fprintf('AvgPositioned     = %.1f\n', AvgPositioned);
            fprintf('AvgCPU            = %.1f\n', AvgCPU);
            fprintf('AvgMaxError(Succ) = %.0e\n', AvgMaxError);
            fprintf('AvgRMSD(Succ)     = %.0e\n', AvgRMSD);
            
            disp('==================================================');
            
            diary off;
                        
        end
        
    end
    
end

fprintf('\nResults stored in file %s\n\n',filename);