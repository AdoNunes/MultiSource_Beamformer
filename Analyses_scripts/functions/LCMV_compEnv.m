function [CorrMats_ort, CorrMats_raw] = LCMV_compEnv(data,data_raw,Settings)
% data: is the source time series leakage corrected
% data raw: the uncorrected source time series
% Setting: has to specify envelope parametersm and regulrization
%
% uses ROInets toolbox from OHBA: https://github.com/OHBA-analysis
%
% this code is an adaptation of some ROInets toolbox scripts

nH0Iter = 10;
ARorder= 1;
if nH0Iter,
    [ARmodel.coeffs,ARmodel.varianceEstimate, ARmodel.partial_coeffs  ] = ...
        ROInets.estimate_AR_coeffs(data_raw, ARorder);
else
    ARmodel = [];
end%if

time = (1:size(data,2))/Settings.fs;


if strcmpi(Settings.leakageCorrectionMethod, 'pairwise'),
    Settings.Regularize.do = 0;
    CorrMats_ort = ROInets.do_pairwise_calculation(data, time, Settings.EnvelopeParams);
    
    nodeEnv_raw = ROInets.envelope_data(data_raw, time,  Settings.EnvelopeParams);
    
    nodeEnv_raw(:,isnan(sum(nodeEnv_raw))) = [];
    CorrMats_raw = ROInets.run_correlation_analysis_noPartial(data_raw,nodeEnv_raw, Settings.Regularize);
    
else
    %% other
    
    nodeEnv     = ROInets.envelope_data(data, time,  Settings.EnvelopeParams);
    nodeEnv_raw = ROInets.envelope_data(data_raw, time,  Settings.EnvelopeParams);
    
    nodeEnv(:,isnan(sum(nodeEnv))) = [];
    CorrMats_ort = ROInets.run_correlation_analysis_noPartial(data,nodeEnv, Settings.Regularize);

    nodeEnv_raw(:,isnan(sum(nodeEnv_raw))) = [];
    CorrMats_raw = ROInets.run_correlation_analysis_noPartial(data_raw,nodeEnv_raw, Settings.Regularize);

end
%%
nNodes              = size( data,1);
nTimeSamples        = size( data,2);
Filter              = struct('band', []);
Filter.band         = [8 12];
Filter.order        = 4;
Filter.type         = 'but';
Filter.direction    = 'twopass';

sigma = ROInets.find_empirical_H0_distribution_width(nH0Iter, nNodes, nTimeSamples,  ...
    Settings,                  ...
    CorrMats_ort.Regularization, ...
    ARmodel,                   ...
    Settings.fs, Filter,                ...
    Settings.EnvelopeParams);

CorrMats_ort.ARmodel    = ARmodel;
CorrMats_ort.H0Sigma    = sigma;
CorrMats_ort.timeWindow = Settings.EnvelopeParams;

CorrMats_raw.ARmodel    = ARmodel;
CorrMats_raw.H0Sigma    = sigma;
CorrMats_raw.timeWindow = Settings.EnvelopeParams;

%% conversion of correlations to z-stats
% nullHypothesisCor = mean(corrMats.envCorrelation(triu(true(size(corrMats.envCorrelation)),1)));
% Values are corrected by the mean and normalized by the STD of the null

doRegularize = Settings.Regularize.do;
CorrMats_ort = ROInets.convert_correlations_to_normal_variables(CorrMats_ort,sigma, doRegularize);
CorrMats_raw = ROInets.convert_correlations_to_normal_variables(CorrMats_raw,sigma, doRegularize);

% average correlations from pairwise analysis - best to average after
% conversion to z-stats.
if strcmpi(Settings.leakageCorrectionMethod, 'pairwise'),
    
    ff = fieldnames(CorrMats_ort);
    for iff = 1:length(ff),
        x = CorrMats_ort.(ff{iff});
        
        if ismatrix(x) && ~isvector(x),
            CorrMats_ort.(ff{iff}) = (x +x')./2.0;
        else
            CorrMats_ort.(ff{iff}) = x;
        end
    end
    
end

%%% correct diagonal
ff = fieldnames(CorrMats_ort);
for iff = 1:length(ff),
    x = CorrMats_ort.(ff{iff});
    if ismatrix(x) && ~isvector(x),
        CorrMats_ort.(ff{iff})(logical(eye(size( CorrMats_ort.(ff{iff}))))) = 0;
    end
end

ff = fieldnames(CorrMats_raw);
for iff = 1:length(ff),
    x = CorrMats_raw.(ff{iff});
    if ismatrix(x) && ~isvector(x),
        CorrMats_raw.(ff{iff})(logical(eye(size( CorrMats_raw.(ff{iff}))))) = 0;
    end
end

end