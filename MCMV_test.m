%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script uses a set of simple sinusoidal and gaussian signals to 
% estimate connectivity using a scalar LVMV and a muti-source (MCMV) BFs.
%  
% To compute connectivity measures it requires:
%       Fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip
%
% Adonay Nunes, SFU, Vancouver, March 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes/MultiSource_Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ft_pat = dir('~/Documents/MATLAB/fieldtrip*');
% addpath (['~/Documents/MATLAB/',ft_pat.name])
ft_defaults
 
load('LF_4Src.mat') % LF is nSrcs x 3 D x nSens
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sinusoidal signals for simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs      = 1000;
dt      = 1/fs;
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
ntrials     = 60;
data_sinSig = reshape(data_sinSig, size(data_sinSig,1), [], ntrials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Send to sensors 

data                        = data_sinSig;
[nsrcs, nsampls, ntrials]   = size(data);
nsens                       = size(LF,3);
Sensor_signal               = squeeze(LF(:,1,:))'* reshape(data, nsrcs,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add noise to sensors 
data_noise_wgn  = wgn(nsens,nsampls*ntrials,1);

Sensor_noise    = data_noise_wgn * ((rms(Sensor_signal(:))/rms(data_noise_wgn(:)))/2.5);

Sensor_data     = Sensor_signal + Sensor_noise;
Sensor_data_d   = reshape(Sensor_data, nsens, nsampls, ntrials);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Do signal & noise cov

SCov = cov(Sensor_data' );
NCov = cov(Sensor_noise');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Beamforming

% MCMV
sIn.iR     = inv(SCov);
sIn.arrH    = LF;    
sIn.arrN    = NCov;
sOut2       = MCMV_BF(sIn);   
W_MCMV      = sOut2.arrW'; 


% LCMV
for p = 1:size(LF,1)
    sIn.arrH    = LF( p ,:,:);  
    sOut1b      = MCMV_BF(sIn);
    W_LCMV(p,:) = sOut1b.arrW';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make FT dataset

variable                = data;
source_data1            = struct();
source_data1.fsample    = fs;
source_data1.time       = {linspace(0,size(variable,2)/fs, size(variable,2))};
source_data1.label      = cellfun( @num2str, num2cell([1:size(variable,1)]'), 'un', 0);
source_data2            = source_data1;
source_data3            = source_data1;

for d = 1:size(variable,3)
    
    source_data1.trial{d} = W_LCMV * Sensor_data_d(:,:,d); 
    source_data2.trial{d} = W_MCMV * Sensor_data_d(:,:,d); 
    source_data3.trial{d} = variable(:, :, d) ; 
    
    source_data1.time{d} = source_data1.time{1};
    source_data2.time{d} = source_data2.time{1};    
    source_data3.time{d} = source_data3.time{1};
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot reconstructed sources
 
s1 = ft_timelockanalysis([],source_data1);
s2 = ft_timelockanalysis([],source_data2);
s3 = ft_timelockanalysis([],source_data3);

t = s1.time - .5;

figure,
for s = 1:nsrcs
    
    subplot(2,2,s), hold on,
    if corr(s1.avg(s,:)',s3.avg(s,:)') <0;   s1.avg(s,:) = -s1.avg(s,:);  end
    if corr(s2.avg(s,:)',s3.avg(s,:)') <0;   s2.avg(s,:) = -s2.avg(s,:);  end
    
    plot(t,(s1.avg(s,:) - mean(s1.avg(s,:))) *1000,'LineWidth',2, 'Color',  [0      0.447 0.741])
    plot(t,(s3.avg(s,:) - mean(s3.avg(s,:))) *1000,'LineWidth',3, 'Color',  [0.929  0.694 0.125 0.7])
    plot(t,(s2.avg(s,:) - mean(s2.avg(s,:))) *1000,'LineWidth',2, 'Color',  [0.85   0.325 0.098])
    
    set(gca, 'fontsize',12, 'LineWidth',1)
end
set(gcf,'color','w');
set(gcf, 'Position', [864   822   917   235])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Connectivity

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foilim    = [0 50];
cfg.pad       = 'nextpow2';
cfg.tapsmofrq = 2;
freq1         = ft_freqanalysis(cfg, source_data1);
freq2         = ft_freqanalysis(cfg, source_data2);
freq3         = ft_freqanalysis(cfg, source_data3);

conn    = { 'amplcorr' 'coh' 'csd' 'dtf' 'granger' 'pdc' 'plv' 'powcorr' 'powcorr_ortho' 'ppc'  'wpli' 'wpli_debiased' 'wppc'};
connout = {'amplcorrspctrm' 'cohspctrm' 'crsspctrm' 'dtfspctrm' 'grangerspctrm' 'pdcspctrm' 'plvspctrm' 'powcorrspctrm' ...
    'powcorrspctrm' 'ppcspctrm' 'wplispctrm' 'wpli_debiasedspctrm' 'wppcspctrm'  };

for c = [2,7]%1:numel(conn);
    
    cfg            = [];
    cfg.method     = conn{c};
    
    cohm1          = ft_connectivityanalysis(cfg, freq1);
    cohm2          = ft_connectivityanalysis(cfg, freq2);
    cohm3          = ft_connectivityanalysis(cfg, freq3);
    
    eval(['c1 = cohm1.',connout{c},';' ])
    eval(['c2 = cohm2.',connout{c},';' ])
    eval(['c3 = cohm3.',connout{c},';' ])
    
    thr = freq1.freq(end);
    
    figure
    for row=1:nsrcs-1
        for col=row:nsrcs
            if row ~= col
                hold off
                subplot(nsrcs,nsrcs,(row-1)*nsrcs+col);
                
                h3 = plot(cohm3.freq, squeeze(c3(row,col,:)));
                hold on
                h1 = plot(cohm1.freq, squeeze(c1(row,col,:)));
                
                set(h1,'LineWidth',2); h1.Color = [0         0.4470    0.7410];
                set(h3,'LineWidth',4); h3.Color = [0.9290    0.6940    0.1250 0.7];
                
                vline([f1 f1]),vline([f2 f2])
                xlim([ 0 thr]),% ylim( [ 0 1])
                set(gca, 'fontsize',15, 'LineWidth', 1.5)
            end
        end
    end
    
    for row=1:nsrcs
        for col=1:row
            if row ~= col
                subplot(nsrcs,nsrcs,(row-1)*nsrcs+col);
                
                h3 = plot(cohm3.freq, squeeze(c3(row,col,:)));
                hold on
                h2 = plot(cohm2.freq, squeeze(c2(row,col,:)));
                
                set( h2 ,'LineWidth',2); h2.Color = [0.8500    0.3250    0.0980];
                set( h3 ,'LineWidth',4); h3.Color = [0.9290    0.6940    0.1250 0.7];
                
                vline([f1 f1]),vline([f2 f2])
                xlim([ 0 thr]),% ylim( [ 0 1])
                set(gca, 'fontsize',15, 'LineWidth', 1.5)
            end
        end
    end
    
    suptitle([conn{c} ])
    set(gcf, 'color','w');
    set(gcf, 'Position', [2320         532         637         458])
end




