function CorrMats_MCMV = MCMV_compEnv(W,data_sensor, Settings, data_raw)
% W: are the beamformer weights from MCMV (n pairs x M channels x 2)
% data raw: arethe uncorrected source time series
% Setting: has to specify envelope parametersm and regulrization
%
% uses ROInets toolbox from OHBA: https://github.com/OHBA-analysis

% this code is an adaptation of some ROInets toolbox scripts


if nargin==4
    compute_sigma = 1;
else
    compute_sigma = 0;
end

for pr = 1:size(W,1);
    data_src= squeeze(W(pr, :,:))'*data_sensor;
    time = (1:size(data_src,2))/Settings.fs;
    
    %% other
    [nodeEnv, time_ds] = ROInets.envelope_data(data_src, time,  Settings.EnvelopeParams);
    
    nodeEnv(:,isnan(sum(nodeEnv))) = [];
    CorrMats = ROInets.run_correlation_analysis_noPartial(data_src,nodeEnv, Settings.Regularize);
    
    CorrMats_all.correlation(pr)   = CorrMats.correlation(1,2);
    CorrMats_all.envCovariance(pr) = CorrMats.envCovariance(1,2);
    CorrMats_all.envCorrelation(pr)= CorrMats.envCorrelation(1,2);
    if isnan( CorrMats.envCorrelation(1,2))
        error (['wrong Corr analysis in pair ', num2str(pr)] )% NaNs are not possible
    end
end


ff = fieldnames(CorrMats_all);
for iff = 1:length(ff),
    x = CorrMats_all.(ff{iff});
    if ismatrix(x) && isvector(x) && numel(x) >1
        CorrMats_all.(ff{iff}) = squareform((x));
    else
        CorrMats_all.(ff{iff}) = x;
    end
end

CorrMats_all.nSamples = CorrMats.nSamples;
CorrMats_all.Regularization = CorrMats.Regularization;

CorrMats_all.ARmodel = Settings.ARmodel;
CorrMats_all.H0Sigma = Settings.sigma;

if compute_sigma
    %% ar
    
    ARorder= 1;
    [ARmodel.coeffs,ARmodel.varianceEstimate, ARmodel.partial_coeffs  ] = ...
        ROInets.estimate_AR_coeffs(data_raw, ARorder);
    
    nNodes              = size(   data,1);
    nTimeSamples        = size(   data,2);
    Filter              = struct('band', []);
    Filter.band         = [8 12 ];
    Filter.order        = 4;
    Filter.type         = 'but';
    Filter.direction    = 'twopass';

    nH0Iter= 10;
    sigma = ROInets.find_empirical_H0_distribution_width(nH0Iter, nNodes, ...
        nTimeSamples, Settings, CorrMats_all.Regularization, ...
        ARmodel,Settings.fs, Filter, Settings.EnvelopeParams);
    
    CorrMats_all.ARmodel                = ARmodel;
    CorrMats_all.H0Sigma                = sigma;
    CorrMats_all.envPartialCorrelation  = eye(nNodes);
end

CorrMats_MCMV = ROInets.convert_correlations_to_normal_variables(CorrMats_all,Settings.sigma, 0);

ff = fieldnames(CorrMats_MCMV);
for iff = 1:length(ff),
    x = CorrMats_MCMV.(ff{iff});
    if ismatrix(x) && ~isvector(x),
        CorrMats_MCMV.(ff{iff})(logical(eye(size( CorrMats_MCMV.(ff{iff}))))) = 0;
    end
end


