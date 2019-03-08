%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script estimates connectivity using a scalar LCMV and a muti-source 
% (MCMV) BFs. Based on a task approach where there are a signal and a noise
% covariance matrices. The noise cov. is used to find optimal and strongest
% orientation of the sources. 
%  
% Based from the study: 
%       https://www.biorxiv.org/content/10.1101/567768v1
%
%
% DEPENDENCIES:
% To get MEG data it requires:
%       Fieldtrip toolbox: ftp://ftp.fieldtriptoolbox.org/pub/fieldtrip
%
% DATA:
% Data obtained from script S01_create_datasets.m
%
% Adonay Nunes, SFU, Vancouver, March 2019
% adonay.s.nunes@gmail.com
% from github: AdoNunes/MultiSource_Beamformer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ft_pat = dir('/Users/adonay/Documents/MATLAB/fieldtrip*');
addpath (['/Users/adonay/Documents/MATLAB/',ft_pat.name])
ft_defaults


load('simulation_datasets.mat')
 
fs = 1000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make sensor dataset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data_ecog;% data_sinSig;%
 
nsrcs  = size(data, 1);
nchans = size(LF,   3);
nsmpls = size(data, 2);
ntrls  = size(data, 3);


data_Noise_o    = genSurrogateData(data_rsMEG); 
data_Noise_c_o  = reshape(data_Noise_o, nchans, []);
data_ori_c      = squeeze(LF(:,1,:))'* reshape(data, nsrcs,[]);
data_ori_c      = data_ori_c * (rms(data_Noise_c_o(:))/rms(data_ori_c(:)));


%% Set noise matrix and signal + noise sensor dataset

SNR          = 2.5; % amplitude SNR.^2 = power SNR

data_Noise_c = data_Noise_c_o / (SNR.^2);

data_ori     =  reshape(data_ori_c, nchans, nsmpls, ntrls);
data_Noise   = reshape(data_Noise_c, nchans, nsmpls, ntrls);

Sensor_data  = data_ori + data_Noise;

%figure, pwelch_d(Sensor_data, [],[],1:0.5:50,fs)

%% Do signal & noisecov

SCov = cov(reshape(Sensor_data, size(Sensor_data,1), size(Sensor_data,2)*size(Sensor_data,3))');
NCov = cov(reshape(data_Noise, size(data_Noise,1), size(data_Noise,2)*size(data_Noise,3))');


%% Beamforming
 
Bavg  = mean(Sensor_data,3);
C_avg = (Bavg * Bavg')/size(Sensor_data,3);


% MPZ
clear sIn  
lambdamat   = 0.001 * trace(SCov)/size(SCov,1); % .1% regularization
sIn.Rm1     = inv(SCov + lambdamat * eye(size(SCov)));
sIn.arrH    = LF;    
sIn.arrN    = NCov;
sIn.Cavg    = C_avg;
sOut2       = constructMCMVWeights22sMER(sIn);   % Calculate the MCMV weights
W2          = sOut2.arrW'; 
Pseudoz2    = sOut2.beamSNR/(size(sIn.arrH,1)*10);

clear W1
Pseudoz1= [];
for p = 1:size(LF,1)
    sIn.arrH    = LF( p ,:,:);    % Get lead fields for our voxels only
    sOut1b      = constructMCMVWeights22sMER(sIn);
    U3D1(:,p)   = sOut1b.U3D;
    W1(p,:)     = sOut1b.arrW';
    Pseudoz1(p) = sOut1b.beamSNR;
end


%% make FT dataset

variable                = data;
source_data1            = struct();
source_data1.fsample    = fs;
source_data1.time       = {linspace(0,size(variable,2)/fs, size(variable,2))};
source_data1.label      = cellfun( @num2str, num2cell([1:size(variable,1)]'), 'un', 0);
nsrc                    = numel(source_data1.label);
source_data2            = source_data1;
source_data3            = source_data1;

for d = 1:size(variable,3)
    
    source_data1.trial{d} = W1 * Sensor_data(:,:,d); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
    source_data2.trial{d} = W2 * Sensor_data(:,:,d); % data_ori(:,:,d);%  Noise_dat_o(:,:,d);%
    source_data3.trial{d} = variable(:, :, d) ; 
    
    source_data1.time{d}  = source_data1.time{1};
    source_data2.time{d}  = source_data2.time{1};    
    source_data3.time{d}  = source_data3.time{1};
    
end

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
    
    plot(t,(s1.avg(s,:) - mean(s1.avg(s,:))) *1000,'LineWidth',2, 'Color',  [0 0.447 0.741])
    plot(t,(s3.avg(s,:) - mean(s3.avg(s,:))) *1000,'LineWidth',3, 'Color',[0.929 0.694 0.125 0.7])  
    plot(t,(s2.avg(s,:) - mean(s2.avg(s,:))) *1000,'LineWidth',2, 'Color',  [0.85 0.325 0.098]) 
    
    set(gca, 'fontsize',12, 'LineWidth',1)
    ylim( [ -20 20])
end
set(gcf,'color','w');
set(gcf, 'Position', [864   822   917   235])
 
 %% Connectivity

cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.foilim    = [1 51];
cfg.pad       =  'nextpow2';
cfg.tapsmofrq = 2;
freq1         = ft_freqanalysis(cfg, source_data1);
freq2         = ft_freqanalysis(cfg, source_data2);
freq3         = ft_freqanalysis(cfg, source_data3);

conn    = {  'coh'  'plv' };
connout = { 'cohspctrm'  'plvspctrm' };

for c = 1
    
    cfg            = [];
    cfg.method     = conn{c};
    
    cohm1          = ft_connectivityanalysis(cfg, freq1);
    cohm2          = ft_connectivityanalysis(cfg, freq2);
    cohm3          = ft_connectivityanalysis(cfg, freq3);
    
    eval(['c1 = cohm1.',connout{c},';' ])
    eval(['c2 = cohm2.',connout{c},';' ])
    eval(['c3 = cohm3.',connout{c},';' ])
    
    fth = 50;
    th  = find(cohm1.freq > fth,1);
   
    figure
    for row=1:nsrc-1
        for col=row:nsrc
            if row ~= col
                hold off
                subplot(nsrc,nsrc,(row-1)*nsrc+col);
                h3 = plot(cohm3.freq, squeeze(c3(row,col,:)));
                hold on
                h1 = plot(cohm1.freq, squeeze(c1(row,col,:)));
                 
                set(h1,'LineWidth',2); h1.Color = [0         0.4470    0.7410];
                set(h3,'LineWidth',4); h3.Color = [0.9290    0.6940    0.1250 0.7];                
                xlim([ 0 fth]), ylim( [ 0 1])
                
                set(gca, 'fontsize',15, 'LineWidth', 1.5)
            end
        end
    end
    
    for row=1:nsrc
        for col=1:row
            if row ~= col
                subplot(nsrc,nsrc,(row-1)*nsrc+col);
                h3 = plot(cohm3.freq, squeeze(c3(row,col,:)));
                hold on
                h2 = plot(cohm2.freq, squeeze(c2(row,col,:)));
                
                set( h2 ,'LineWidth',2); h2.Color = [0.8500    0.3250    0.0980];
                set( h3 ,'LineWidth',4); h3.Color = [0.9290    0.6940    0.1250 0.7];
                
                xlim([ 0 fth]), ylim( [ 0 1])
                set(gca, 'fontsize',15, 'LineWidth', 1.5)
            end
        end
    end
    
  %suptitle([conn{c} , ' SNR ', num2str(SNR)])
    set(gcf,'color','w');
    set(gcf, 'Position', [493   596   637   458])
    
end


 
