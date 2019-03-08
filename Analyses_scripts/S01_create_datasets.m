%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script creates three datasets to be used as simulated signals and
% noise from real MEG resting data. 
% estimate connectivity using a scalar LVMV and a muti-source (MCMV) BFs.
%  
% Based from the study: 
%       https://www.biorxiv.org/content/10.1101/567768v1
%
% DEPENDENCIES:
% To preprocess MEG data it requires:
%       Fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip
%
%
% DATA:
% Data can be obtained from:
%       ECoG:        https://exhibits.stanford.edu/data/catalog/zk881ps0522
%       MEG resting: https://www.humanconnectome.org/study/hcp-young-adult
%
% Adonay Nunes, SFU, Vancouver, March 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes/MultiSource_Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ft_pat = dir('~/Documents/MATLAB/fieldtrip*');
addpath (['~/Documents/MATLAB/',ft_pat.name])
ft_defaults

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Ecog data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the subject
subj        = 'bp';
Ecog_folder = '/Users/adonay/Desktop/projects/MCMV/simulations/';

load([Ecog_folder,'intra-eeg/dg/data/', subj,'/',subj,'_fingerflex.mat']), 
load([Ecog_folder,'intra-eeg/dg/data/', subj,'/',subj,'_stim.mat'])
 
data = detrend(car(data)); % rereference to avg
data = data(1:end-40,:);
stim = stim(1:end-40,:);
flex = flex(1:end-40,:);
 
% take thumb condition only
Stim_4use               = stim;
Stim_4use(Stim_4use~=1) = -1; 
Stim_4use(stim ==0)     = 0;
 
flex_4use               = zeros(size(flex,1),1) - 1;
flex_4use(Stim_4use==1) = flex(Stim_4use==1,1);
flex_4use(Stim_4use ==0)= 0;

flex_4use = (flex_4use - mean(flex_4use))/std(flex_4use);

% get peaks above STD
STD                         = 2.75;
peaks_vec                   = ones(size(flex,1),1)*STD;
peaks_vec(flex_4use >STD)   = flex_4use(flex_4use >STD);
[~,plocs]                   = findpeaks(peaks_vec);

% reject  peaks too close
dst     = 200;
locdis  =[];
for p1 = 1:numel(plocs)-1; 
  locdis(p1) =  plocs(p1+1)- plocs(p1);
end
   
%figure, hold on, plot(Stim_4use), plot(flex_4use), hline([STD STD])
%plot(plocs, peaks_vec(plocs), '*')
%plot(plocs(locdis<dst), peaks_vec(plocs(locdis<dst)), 'ro')
    
rm_pk = (find (locdis<dst))+1; % remove peaks too close
plocs(rm_pk) = [];

% peak locations to use
pks = [20 14 2 5 ];
pks = [12 14 10 29 ];

% notch filter data  
data_notch  = data(:,pks);
f2notch     = [60 120 180 ];
fs          = 1000;      % sampling rate
fn          = fs/2;      % Nyquist frequency

 for fk = 1:numel(f2notch)
    f0              = f2notch(fk);                     
    freqRatio       = f0/fn;                    
    notchWidthRatio = freqRatio/35;% width of the notch

    [b,a]           = iirnotch(freqRatio,notchWidthRatio);

    data_notch      = filter(b,a,data_notch);
end 

data_notch = data_notch';
data_trial = [];
for d = 1:60
    data_trial(:,:,d ) = data_notch(:,plocs(d)- 500: plocs(d)+499);
end

data_ecog = data_trial;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% resting MEG data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HCP_subj_folder = '/Volumes/4TB_drive/HCP/MEG_dat/100307/';

load([HCP_subj_folder,'MEG/anatomy/100307_MEG_anatomy_headmodel.mat'])
load([HCP_subj_folder,'MEG/anatomy/100307_MEG_anatomy_sourcemodel_3d6mm.mat'])

sourcemodel3d = ft_convert_units(sourcemodel3d, 'mm');

RS_Data = load([HCP_subj_folder,'/MEG/Restin/rmegpreproc/100307_MEG_3-Restin_rmegpreproc.mat']);
RS_Data = RS_Data.data;

cfg             = [];
cfg.channel     = 'MEG';
cfg.bpfilter    = 'yes';
cfg.bpfreq      = [1 55];
cfg.padding     =  4;
RS_Data         = ft_preprocessing(cfg,RS_Data);
    
cfg             = [];
cfg.resamplefs  = 1000;
RS_Data         = ft_resampledata(cfg, RS_Data);

cfg             = [];
cfg.trials      = 1:size(data_ecog,3);
cfg.latency     = [RS_Data.time{1}(1003) RS_Data.time{1}(2002)];
RS_Data         = ft_selectdata(cfg,RS_Data);

data_rsMEG      = cat(3,RS_Data.trial{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% get ecog LF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

template_grid   = sourcemodel3d.cfg.grid.template;
template_grid   = ft_convert_units(template_grid, 'mm');
[i, d]          = knnsearch(template_grid.pos, locs(pks,:));

LF_pos  = sourcemodel3d.pos(i,:);
LF_dist = squareform(pdist(LF_pos));

evalc(fileread([HCP_subj_folder,'MEG/anatomy/100307_MEG_anatomy_transform.txt']));

T                       = transform.vox07mm2bti/transform.vox07mm2spm;
template_grid           =  ft_transform_geometry(T, template_grid);
template_headmodel      =  ft_transform_geometry(T, template_headmodel);

template_grid.inside    = false(size(template_grid.inside));
template_grid.inside(i) = true;
sourcemodel3d.inside    = template_grid.inside;

cfg             = [];
cfg.headmodel   = headmodel;
cfg.grad        = RS_Data.grad;
cfg.grid        = sourcemodel3d;
cfg.normalize   = 'no';
cfg.reducerank  = 'no';
cfg.backproject = 'no';
lf              = ft_prepare_leadfield(cfg, RS_Data);

LF = cat(3,lf.leadfield{i});
LF = permute(LF, [3 2 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sinusoidal signals dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

winlen  = 1000;
dt      = 1/winlen;
t       = 0:dt:(60000*dt)-dt;
f1      = 10;
f2      = 35;

x1 = sin(2*pi*f1*t) +sin(2*pi*f2*t)  +randn(size(t));
x2 = sin(2*pi*f1*t) +sin(2*pi*f2*t)  +randn(size(t));
x3 =                +sin(2*pi*f2*t)  +randn(size(t));
x4 = sin(2*pi*f1*t)                  +randn(size(t));

% scale to same rms
x1 = x1* 1e-8;
x2 = x2 * (rms(x1)/rms(x2));
x3 = x3 * (rms(x1)/rms(x3));
x4 = x4 * (rms(x1)/rms(x4));

data_sinSig = [x1;x2;x3;x4];
data_sinSig = reshape(data_sinSig, size(data_sinSig,1) , size(data_ecog,2),size(data_ecog,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save work done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save( 'simulation_datasets.mat', 'data_ecog', 'data_rsMEG', 'LF', 'LF_pos', 'LF_dist', 'data_sinSig')
 

