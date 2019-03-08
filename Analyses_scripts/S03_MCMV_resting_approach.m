%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script estimates connectivity using a scalar LCMV and a muti-source 
% (MCMV) BFs. Based on a resting state approach where there is only a 
% signal covariance matrix (no noise cov.)
% The approach is either:
%       a PairWise MCMV: for each pair of sources, time series are 
%           reconstructed and connectivity estimated 
%       an Augmented PairWise MCMV: for each pair extra sources are added 
%           to further reduce leakage from neighboring sources
%
% The choice of connectivity is envelope correlations, and these are
% estimated with an LCMV, PW-MCMV, APW-MCMV and from the symmetrically
% orthogonalized LCMV time series (Colclough 2015, Neuroimage). In here,
% the data is surrogated, thus all significant connectivity is spurious!
%  
% Based from the study: 
%       https://www.biorxiv.org/content/10.1101/567768v1
%
%
% DEPENDENCIES:
% To preprocess MEG data it requires:
%       Fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip
%
% To calculate envelope correlations and plot brains it requires:
%       OHBA toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip
%
% DATA:
% HCP data can be obtained from:
%       MEG resting: https://www.humanconnectome.org/study/hcp-young-adult
%
% Adonay Nunes, SFU, Vancouver, March 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes/MultiSource_Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% add necessary paths

addpath(genpath('~/Documents/MATLAB/OHBA-analysis'))
ft_pat = dir('~/Documents/MATLAB/fieldtrip*');
addpath (['~/Documents/MATLAB/',ft_pat.name])
ft_defaults
addpath('functions')

%% set dirs
outdir              = 'Res_resting_analysis/';
if ~exist(outdir, 'dir'), mkdir(outdir), end

subdir              = '/Volumes/4TB_drive/HCP/MEG_dat/';
subject             = dir(subdir);
not_fls             = ismember({subject.name}, { '.', '..', '.DS_Store'});
subject(not_fls)    = [];


%% make AAL template

temp          = load (sprintf('%s%s/MEG/anatomy/%s_MEG_anatomy_sourcemodel_3d6mm.mat', subdir. subject(1).name, subject(1).name));  
template_grid = temp.sourcemodel3d.cfg.grid.template;
template_grid = ft_convert_units(template_grid, 'mm');

coord_aal   = importdata('functions/aal116_COG.txt');
lab_aal     = importdata('functions/aal116_LABELS.txt');
rois        = [1:36, 43:70,79:90]; % only cortical
coord_aal   = coord_aal(rois,:);

A           = template_grid.pos;
B           = coord_aal;
[IDX, dis]  = knnsearch(A,B); % find closest grid point

template_grid.inside        = zeros(size(template_grid.pos,1),1);
template_grid.inside(IDX)   = 1;
template_grid.inside        = logical(template_grid.inside);
 

%% loop for each subject
subject([3, 13, 18, 35,45,67]) =[]; % out subj. without resting

for s = 1:numel(subject)
    
    
    load(sprintf('%s%s/MEG/anatomy/%s_MEG_anatomy_headmodel.mat', subdir, subject(s).name, subject(s).name))
    load(sprintf('%s%s/MEG/anatomy/%s_MEG_anatomy_sourcemodel_3d6mm.mat', subdir, subject(s).name, subject(s).name));
    
    headmodel     = ft_convert_units(headmodel, 'mm');
    sourcemodel3d = ft_convert_units(sourcemodel3d, 'mm');
    
    sourcemodel3d.inside = template_grid.inside;
    load(sprintf('%s%s/MEG/Restin/rmegpreproc/%s_MEG_3-Restin_rmegpreproc.mat',subdir,subject(s).name, subject(s).name));
    

    %% band pass filter  signal
    cfg             = [];
    cfg.channel     = 'MEG';
    cfg.bpfilter    = 'yes';
    cfg.bpfreq      = [8 12];
    cfg.padding     =  6;
    data            = ft_preprocessing(cfg,data);
    
    data_sen        = cat(3, data.trial{:});
    data_sen_c      = genSurro(reshape(data_sen, size(data_sen,1), []));

    %% make lf
    cfg                 = [];
    cfg.headmodel       = headmodel;
    cfg.grad            = data.grad;
    cfg.grid            = sourcemodel3d;
    cfg.reducerank      = 3;
    cfg.backproject     = 'no';
    lf                  = ft_prepare_leadfield(cfg, data);
    
    LF = cat(3,lf.leadfield{IDX});
    LF = permute(LF,[3 2 1]);
    
    %% make cov
    cfg                     = [];
    cfg.trials              = 'all';
    cfg.covariance          = 'yes';
    cfg.covariancewindow    = 'all';
    cfg.keeptrials          = 'yes';
    cfg.removemean          = 'yes';
    timelock                = ft_timelockanalysis(cfg, data);
    
    Cov_S =squeeze(mean(timelock.cov));
    
    clear sIn
    lambdamat = 0.01 * trace(Cov_S)/size(Cov_S,1);
    sIn.Rm1 = inv(Cov_S + lambdamat * eye(size(Cov_S)));
    
    %% env correlation settings
    Settings                                = struct();
    Settings.fs                             = data.fsample ;
    Settings.EnvelopeParams.windowLength    = 2;
    Settings.EnvelopeParams.overlap         = 0.6;
    Settings.useFilter                      = 1;
    Settings.Regularize.do                  = 0;
    Settings.leakageCorrectionMethod        = 'sym';
    
    %% Calculate BF weights and connectivity
    %%%%%%%%%%%%%%%%%%
    %%%% LCMV      %%%
    %%%%%%%%%%%%%%%%%%
    
    clear W1
    for p = 1:size(LF,1)
        sIn.arrH = LF( p,:,:); 
        sOut1b   = constructMCMVWeights22sMER(sIn);
        W1(p,:)  = sOut1b.arrW';
    end
    
    % LCMV source data
    data_src_c    = W1 * data_sen_c;
    % SO source data
    data_src_orth = ROInets.remove_source_leakage(data_src_c,'symmetric');
    
    % envelope correlations
    [CorrMats_sym, CorrMats_LCMV] =  LCMV_compEnv(data_src_orth,data_src_c,Settings);
    save([outdir,subject(s).name, '_sym_raw.mat'],'CorrMats_sym' , 'CorrMats_LCMV')
    
    
    %%%%%%%%%%%%%%%%%%
    %%%% PW-MCMV   %%%
    %%%%%%%%%%%%%%%%%%

    [W ]= run_2sBF(LF, Cov_S,0.01);
    
    
    Settings.ARmodel    = CorrMats_sym.ARmodel;
    Settings.sigma      = CorrMats_sym.H0Sigma;
    
    [CorrMats_PW_MCMV]     =  MCMV_compEnv(W,data_sen_c, Settings);
    save([outdir,subject(s).name, '_MCMV.mat'],'CorrMats_PW_MCMV' )
    
    %%%%%%%%%%%%%%%%%%%
    %%%% APW-MCMV   %%%
    %%%%%%%%%%%%%%%%%%%
    
    h_m            = ROInets.false_discovery_rate(CorrMats_PW_MCMV.envCorrelation_z, 0.05);
    pos            = sourcemodel3d.pos(IDX,:);
    
    W              = run_APW_MCMV(LF, Cov_S, pos, 0.01, [],h_m); 

    CorrMats_APW_MCMV = MCMV_compEnv(W,data_sen_c, Settings);
    save([outdir,subject(s).name,'_APW_MCMV.mat'],'CorrMats_APW_MCMV')
    
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group Stats
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
corrMat= struct;
for s = 1:numel(subject)
   
        load([outdir,subject(s).name, '_MCMV.mat'])
        corrMat.PW_MCMV.raw(:,:,s)       =  CorrMats_PW_MCMV.envCorrelation;
        corrMat.PW_MCMV.raw_z(:,:,s)     =  CorrMats_PW_MCMV.envCorrelation_z;
         
        load([outdir,subject(s).name, '_sym_raw.mat'])
        corrMat.LCMV.raw(:,:,s)         =  CorrMats_LCMV.envCorrelation;
        corrMat.LCMV.raw_z(:,:,s)       =  CorrMats_LCMV.envCorrelation_z;
        
        corrMat.SO.raw(:,:,s)           =  CorrMats_SO.envCorrelation;
        corrMat.SO.raw_z(:,:,s)         =  CorrMats_SO.envCorrelation_z;
        
        load([outdir,subject(s).name, '_APW_MCMV.mat'])
        corrMat.APW_MCMV.raw(:,:,s)     =  CorrMats_APW_MCMV.envCorrelation;
        corrMat.APW_MCMV.raw_z(:,:,s)   =  CorrMats_APW_MCMV.envCorrelation_z;
end


bf_str = {'LCMV','PW_MCMV','APW_MCMV','SO'};
mn = numel(bf_str);
%% 'fixed-effects'
% this should be the first method to try, the z-score env. correlations are
% are zero mean and scaled to the STD of the null ( leakage-free env. corr
% distributions generated with an AR model. Significant connections will be
% different from the mean of the network connenctions and normalized by the
% distribution of the connectivity under the null. 

for m =1:mn
    mtx                         = corrMat.(bf_str{m}).raw_z;
    corrMat.(bf_str{m}).avg_z   =  mean(mtx, 3); 
end

cgram = clustergram(squareform(pdist(coord_aal))); % distance clustering
lr    = cell2mat(cellfun(@str2num,  cgram.ColumnLabels, 'un', 0));

figure('name' , 'Average Z scores'),
lbs = lr;
cr  = [floor(sqrt(mn)), ceil(sqrt(mn))];
for m =1:mn
    mt_pot  = corrMat.(bf_str{m}).avg_z(lbs,lbs);
    ttl     = [bf_str{m}, ' Z avg'];
    min_max = [percentile(mt_pot,5) percentile(mt_pot,95)];
    subplot(cr(1),cr(2), m),imagesc(mt_pot),title(ttl), colorbar;caxis(min_max),  axis square; set(gca,'fontsize', 18)
end
colormap(bluewhitered(256))

%%%%%%%%%%%%%%%%%%%%%
%%% FDR correction
%%%%%%%%%%%%%%%%%%%%%

for m =1:mn
h.(bf_str{m})    = ROInets.false_discovery_rate(corrMat.(bf_str{m}).avg_z, 0.05);
end

%%%%%%%%%%%%%%%%%%%%%
%%% Brain plot
%%%%%%%%%%%%%%%%%%%%%

nSrcs       = size(coord_aal,1);
sphereCols  = repmat([30 144 255]/255, nSrcs, 1);
sphereSiz   = repmat(0.5,nSrcs,1);
for m =1:mn
    
    graph           = corrMat.(bf_str{m}).avg_z .* h.(bf_str{m});
    graph(graph==0) = NaN;
    
    th_d      = min(graph(~isnan(graph))) - (min(graph(~isnan(graph)))/10);
    th_u      = percentile(graph(~isnan(graph)), 99);
    colorLims = [th_d th_u];
    
    subplot(2,mn,m);    osl_braingraph(graph, colorLims, sphereSiz, [0 1], coord_aal, [], 0, sphereCols,[3 6]); view([0 90]),
    subplot(2,mn,m+mn); osl_braingraph(graph, colorLims, sphereSiz, [0 1], coord_aal, [], 0, sphereCols,[3 6]); view([-90 0]), colorbar off
    
end


%% 'mixed-effects'
% This should be used only if the fixed-effects model did not work as a
% work around. It will detect connections that deviate from the mean of all
% the network connections, considered to be noise from BF numerical 
% inaccuracies. The p vals are Z transformerd and in here the saturation 
% threshold used is 8. 

nSrcs    = size(coord_aal,1);
for m =1:mn
    
    mtx     = corrMat.(bf_str{m}).raw;
    
    corrMat.(bf_str{m}).avg =  mean(corrMat.(bf_str{m}).raw, 3); 
    
    nullMean                = mean( corrMat.(bf_str{m}).avg(triu(ones(nSrcs),1) ~=0));
    [~, p]                  = ttest(shiftdim(corrMat.(bf_str{m}).raw, 2), nullMean, .5, 'right' );
    signZ                   = sign(corrMat.(bf_str{m}).avg);
    corrMat.(bf_str{m}).z   = ROInets.p_to_z_two_tailed(shiftdim(p), signZ);
end


figure
cr = [floor(sqrt(mn)), ceil(sqrt(mn))];
for m =1:mn
    mt_pot  = corrMat.(bf_str{m}).z ;
    ttl     = [bf_str{m}, ' avg '];
    min_max = [percentile(mt_pot,5) percentile(mt_pot,95)];
    subplot(cr(1),cr(2), m),imagesc(mt_pot),title(ttl), colorbar;  axis square; set(gca,'fontsize', 18)
end
colormap(bluewhitered(256))


%%%%%%%%%%%%%%%%%%%%%
%%% FDR correction
%%%%%%%%%%%%%%%%%%%%%

for m =1:mn
    h.(bf_str{m})    = ROInets.false_discovery_rate(corrMat.(bf_str{m}).z, 0.05);
end


%%%%%%%%%%%%%%%%%%%%%
%%% Plot
%%%%%%%%%%%%%%%%%%%%%

nSrcs       = size(coord_aal,1);
sphereCols  = repmat([30 144 255]/255, nSrcs, 1);
sphereSiz   = repmat(0.5,nSrcs,1);

sat_thr = 8;
upp_thr = 99;
figure
for m =1:mn
    
    graph=abs(corrMat.(bf_str{m}).z) .* h.(bf_str{m});
    graph(graph<sat_thr) = NaN;
    th_d = min(graph(~isnan(graph))) - (min(graph(~isnan(graph)))/10);
    th_u = percentile(graph(~isnan(graph)), upp_thr);
    colorLims = [th_d th_u];
    
    subplot(2,mn,m);    osl_braingraph(graph, colorLims, sphereSiz, [0 1], coord_aal, [], 0, sphereCols,[1 3]); view([0 90]), colorbar off
    subplot(2,mn,m+mn); osl_braingraph(graph, colorLims, sphereSiz, [0 1], coord_aal, [], 0, sphereCols,[1 3]); view([-90 0]), colorbar off
    
end


