function safe_stat_contrasts(SPM, opt_con, run)
% 
% SAFE_STAT_CONTRASTS creates the contrasts for one subject. It checks for
% the estimability of contrasts before creating them.
% 
% 
%USAGE
%-----
% safe_stat_contrasts(spm, use_rp, run)
% 
% 
%INPUT
%-----
% - SPM    : SPM.mat file after Model Estimation
% - OPT_CON: structure with the following fields:
%   - SID: subject identification
% - RUN    : 'run' or 'norun'
% 
% 
%OUTPUT
%------
% - SID-contrasts.mat: batch file in the same folder as SPM.mat
% - If RUN is 'run', SPM.mat will be updated.
% 
% See also SAFE_PREPROC, SAFE_STAT_DESIGN, SAFE_STAT_RESULTS
% 
%__________________________________________________________________________
% Copyright (C) 2012-2014 Guilherme Coco Beltramini

% Guilherme Coco Beltramini - 2013-Jan-25


% "Run" option
%--------------------------------------------------------------------------
switch run
    case 'run'
        run = 1;
    case 'norun'
        run = 0;
    otherwise
        error('Invalid "run" option')
end


sID = opt_con.sID;


tmp = clock;
disp('-------------------------------------------------------------------')
fprintf('(%.2d:%.2d:%02.0f) Running %s.m...\n',...
    tmp(4),tmp(5),tmp(6),mfilename)


% SPM.mat
%--------------------------------------------------------------------------
fpath = fileparts(SPM);
if isempty(fpath)
    fpath = pwd;
    SPM = fullfile(fpath,'SPM.mat');
end


matlabbatch{1}.spm.stats.con.spmmat = {SPM};
matlabbatch{1}.spm.stats.con.delete = 1;


% Load SPM.mat
%--------------------------------------------------------------------------
try
    load(SPM)
catch
    warning('SPM.mat not found. You must specify and estimate the model first')
    SPM = struct('xX',struct('X',0),'Sess',[]);
end


% Total number of columns and sessions in the design matrix
%--------------------------------------------------------------------------
C     = size(SPM.xX.X,2); % = size(SPM.Vbeta,2) = ...
Nsess = size(SPM.Sess,2);

contr = 0; % number of contrasts
for ss = 1:Nsess % loop over sessions
    
    for cc = 1:size(SPM.Sess(1,ss).Fc,2) % loop over conditions for each session
        name      = SPM.Sess(1,ss).Fc(1,cc).name;
        same_cond = 0;
        for curr_contr=1:contr % compare current condition name with others
            if strcmp(consess{curr_contr}.tcon.name,name) % same condition
                same_cond = 1;
                break
            end
        end
        
        cols = SPM.Sess(1,ss).Fc(1,cc).i; % columns with reference to the session
        cols = SPM.Sess(1,ss).col(cols);  % columns with reference to the whole scan
        
        % Condition exists already
        %------------------------------------------------------------------
        if same_cond
            consess{curr_contr}.tcon.convec(cols(1)) = 1;
            if size(cols,2)>1 % more than 1 basis function
                for row=1:size(cols,2)
                    consess{curr_contr}.fcon.convec{1}(row,cols(row)) = 1;
                end
            end
        
        % New condition
        %------------------------------------------------------------------
        else
            contr = contr + 1;
            
            % Insert in alphabetic order
            %---------------------------
            aux_names{contr,1} = name;
            tmp = contr;
            while ~issorted(aux_names)
                aux_names([tmp-1 tmp],1) = aux_names([tmp tmp-1],1);
                tmp = tmp - 1;
            end
            if contr>1
                consess = [consess(1:tmp-1) struct('tcon',struct('name',name,'convec',zeros(1,C))) consess(tmp:end)];
            else
                consess{1}.tcon.name   = name;
                consess{1}.tcon.convec = zeros(1,C);
            end
            
            consess{tmp}.tcon.convec(cols(1)) = 1;
            if size(cols,2)>1 % more than 1 basis function
                consess{tmp}.fcon.name   = [sID ' - ' name ' <> baseline (all basis functions) - avg'];
                consess{tmp}.fcon.convec = {zeros(size(cols,2),C)};
                for row=1:size(cols,2)
                    consess{tmp}.fcon.convec{1}(row,cols(row)) = 1;
                end
            end
            
        end
        
    end
    
end


if contr>0 % there are contrasts

    % Organize contrasts and create negative BOLD contrasts
    %======================================================================
    if isfield(consess{1},'fcon')

        % More than 1 basis function
        %------------------------------------------------------------------

        for tmp=contr:-1:1

            % Negative BOLD (t-contrast)
            %--------------------------------------------------------------
            consess{3*tmp  }.tcon.name    = [sID ' - ' consess{tmp}.tcon.name ' < baseline (avg)'];
            consess{3*tmp  }.tcon.convec  = -consess{tmp}.tcon.convec;
            consess{3*tmp  }.tcon.sessrep = 'none';

            % Positive BOLD (t-contrast)
            %--------------------------------------------------------------
            consess{3*tmp-1}.tcon.name    = [sID ' - ' consess{tmp}.tcon.name ' > baseline (avg)'];
            consess{3*tmp-1}.tcon.convec  = consess{tmp}.tcon.convec;
            consess{3*tmp-1}.tcon.sessrep = 'none';

            % F-contrast
            %--------------------------------------------------------------
            consess{3*tmp-2}.fcon         = consess{tmp}.fcon;
            consess{3*tmp-2}.fcon.sessrep = 'none';

            if tmp~=1, consess{tmp} = {};
            else consess{1} = rmfield(consess{1},'tcon');
            end
        end

    else

        % Only 1 basis function
        %------------------------------------------------------------------

        for tmp=contr:-1:1

            % Negative BOLD (t-contrast)
            %--------------------------------------------------------------
            consess{2*tmp}.tcon.name      = [sID ' - ' consess{tmp}.tcon.name ' < baseline (HRF) - avg'];
            consess{2*tmp}.tcon.convec    = -consess{tmp}.tcon.convec;
            consess{2*tmp}.tcon.sessrep   = 'none';

            % Positive BOLD (t-contrast)
            %--------------------------------------------------------------
            consess{2*tmp-1}.tcon.name    = [sID ' - ' consess{tmp}.tcon.name ' > baseline (HRF) - avg'];
            consess{2*tmp-1}.tcon.convec  = consess{tmp}.tcon.convec;
            consess{2*tmp-1}.tcon.sessrep = 'none';

        end

    end


    % Contrast for movement parameters as the last contrast
    %----------------------------------------------------------------------
    if ~isempty(SPM.Sess(1,1).C.C)
        % no real. pars.  : SPM.Sess(1,1).C.C=[]  & SPM.Sess(1,1).C.name=cell(1,0)
        % with real. pars.: SPM.Sess(1,1).C.C=Nx6 & SPM.Sess(1,1).C.name={'R1' 'R2' 'R3' 'R4' 'R5' 'R6'}
        % check only the first session because the others should also have
        contr                       = size(consess,2) + 1;
        consess{contr}.fcon.name    = [sID ' - Movement'];
        consess{contr}.fcon.convec  = {zeros(6,C)};
        consess{contr}.fcon.sessrep = 'none';
        for s=1:Nsess
            % For every session, the realignment parameters are in these columns:
            % SPM.Sess(1,s).col(end-5:end)
            tmp = SPM.Sess(1,s).col(end-5);
            consess{contr}.fcon.convec{1}(1,tmp  ) = 1;
            consess{contr}.fcon.convec{1}(2,tmp+1) = 1;
            consess{contr}.fcon.convec{1}(3,tmp+2) = 1;
            consess{contr}.fcon.convec{1}(4,tmp+3) = 1;
            consess{contr}.fcon.convec{1}(5,tmp+4) = 1;
            consess{contr}.fcon.convec{1}(6,tmp+5) = 1;
        end
    end


    % Delete invalid contrasts
    %======================================================================
    if isfield(SPM.xX,'xKXs')
        est = spm_SpUtil('IsCon',SPM.xX.xKXs); % parameter estimability
        % 0 where parameter is not good
        c = 0;
        while c<size(consess,2) %c<contr
            c = c + 1;
            if isfield(consess{c},'tcon') % T contrast
                ctr = logical(consess{c}.tcon.convec);
                tmp = 'T';
            else % F contrast
                ctr = logical(sum(consess{c}.fcon.convec{1},1));
                tmp = 'F';
            end
            if ~isequal(and(est,ctr),ctr) % 0 in "est" coincides with 1 in "ctr"
                switch tmp
                    case 'T'
                        ctr = est.*ctr;
                        % same as and(est,ctr), but gives double instead of logical
                        if any(ctr) % at least one "1" in "ctr"
                            consess{c}.tcon.convec = ctr;
                        else % "ctr" contains only zeros
                            consess(c) = [];
                            %contr      = contr - 1;
                            c          = c - 1;
                        end
                    case 'F'
                        for ii=size(consess{c}.fcon.convec{1},1):-1:1
                            ctr = consess{c}.fcon.convec{1}(ii,:).*est;
                            if any(ctr)
                                consess{c}.fcon.convec{1}(ii,:) = ctr;
                            else
                                consess{c}.fcon.convec{1}(ii,:) = [];
                            end
                        end
                        if size(consess{c}.fcon.convec{1},1)==0 % all rows were deleted
                            consess(c) = [];
                            %contr      = contr - 1;
                            c          = c - 1;
                        end
                end
            end
        end
        
        
        matlabbatch{1}.spm.stats.con.consess = consess;
        
        
        % Save and run
        %==================================================================
        save(fullfile(fpath,[sID '-contrasts.mat']),'matlabbatch')
        if isempty(consess) % removed all the contrasts
            warning('%s.m was skipped',mfilename)
        elseif run
            %spm_jobman('initcfg');
            spm_jobman('run',matlabbatch)
        end
        
    else % probably there was a problem with the model estimation
        warning('%s.m was skipped',mfilename)
    end


else % no contrasts (contr=0)
    warning('%s.m was skipped',mfilename)
end


tmp = clock;
fprintf('(%.2d:%.2d:%02.0f) %s.m done!\n',...
    tmp(4),tmp(5),tmp(6),mfilename)
disp('-------------------------------------------------------------------')