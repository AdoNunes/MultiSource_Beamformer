function [ W, U ] = run_2sBF( arrH, arrR, lambda, arrN )
%
% Computes multisource beamformer Weights
%
%  arrH: (nSrc x 3 X M sensors), the forward models for the sources
%  arrR: (M x M) sensor covariance
%  lambda: regularization parameters, for rank deficient arrR
%  arrN: (M x M) noise covariance

if nargin<3
    lambda = 0;
    arrN = [];
end

if exist( 'arrN', 'var')
    sIn.arrN = arrN;
else
    arrN = [];
end
n_source  = 2;  
n_chan   = size(arrH,3);
source_pair = nchoosek(1:size(arrH,1), n_source);% very slow 

lambdamat   = lambda * trace(arrR)/size(arrR,1);
sIn.Rm1     = inv(arrR + lambdamat * eye(size(arrR))); % Invert the covariance

W           = zeros(size(source_pair,1),n_chan,n_source);
U           = zeros(size(source_pair,1),3,n_source);
for ii = 1: size(source_pair,1)
    
    lstVoxNums = source_pair(ii,:);
    sIn.arrH = arrH(lstVoxNums,:,:);    % Get lead fields for the source pairs only

    sOut = constructMCMVWeights22sMER(sIn);   % Calculate the MCMV weights
    
%     % TEST: Verify that MCMV constraints W' * H = I are satisfied
%     M = size(sOut.arrW, 1);
%     H = zeros(M, nSrc);
%     
%     for i=1:nSrc
%         H(:,i) = squeeze(sIn.arrH(i,:,:))' * squeeze(sOut.U3D(:,i));
%     end

    W(ii,:,:) = sOut.arrW;
    U(ii,:,:) = sOut.U3D;
end

end
