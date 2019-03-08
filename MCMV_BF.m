function sOut = MCMV_BF(sIn)
%
%
% This function generates a set of MCMV (multi-source) scalar beamformer
% weights for a given set of N sources (nSrc) forward solutions (FS). 
% First for each source, saliences are computed independently, then from
% strongest to weakest contribution of the strongest source is removed to 
% find the orientation of the next strongest source.
% This iteration procedure is described in :
%     Moiseev et al., Neuroimage, 2011, v.58, p.481-496.
%
% The function produces a matrix W of beamformer weights (nSrc x nSensors) 
% and a matrix U of sources orientations (3 x nSrc). Given the constraints
% of unit gain within sources and zero gain across sources:
%     W * (FS*U') = Identity;
%
% This function also provides a Signal to Noise Ratio (SNR) between the 
% signal source covariance and the full (signal + noise) source covariance.
% The SNR can be used to localize the strongest sources in a grid or surface.
% In this case, run this function for each source to get a one source SNR.
%
%
% Input:
%   sIn         - a structure with input arguments, with the fields:
%                 
%                 : arrH    (nSrc x 3 X nSens)  the leadfields (or FS) for
%                           each of the sources. The FS are oriented in the
%                           x,y,z dimensions 
%                 : iR     (nSens x nSens) the inverse of the sensor 
%                           covariance matrix (the task period)
%                           
%                 Optional:
%                 : arrN    (nSens x nSens) noise covariance (the intertial
%                           or rest period). If not supplied, a diagonal 
%                           noise cov. will be used
%                 : Cavg    (nSens x nSens) Matrix of 2nd moments of the
%                           epoch-averaged sensor signals (Bavg*Bavg') 
%
%
% Output:
%   sOut        - a structure with the output results, with the fields:
%
%                 : U3D     (3 x nSrc) the found orientations of the sources
%                 : arrW    (nSens x nSrc) the MCMV weights of the sources
%                 : order   (1 x nSrc) the order of strongest sources
%                 : beamSNR (double) the beamformer joint power "SNR"
%                           value for the set of nSrcs. 
%
% Originally written by A.Moiseev, BCNI, October 2017


arrH = sIn.arrH;

[nSrc, ~, nSens] = size(arrH);


if ~isfield(sIn, 'arrN')
    % Noise covariance not specified - generate a diagonal one
    arrN = nSens*eye(nSens)/trace(sIn.iR); % sets noise power ~ mean EV of iR
else
    arrN = sIn.arrN;  
end

 
if ~isfield(sIn, 'Cavg')% SNR = MPZ
    iR          = sIn.iR;
    Is_Evoked   = 0;
else  % Event related BF, SNR = MER
    iR_C_iR     = sIn.iR * sIn.Cavg * sIn.iR;
    iR          = iR_C_iR;
    Is_Evoked   = 1;
end

iR_N_iR =  sIn.iR * arrN * sIn.iR;

%% First iterate independently to get single source SNR
for iSrc = 1:nSrc % get U3D for each source
    h = squeeze(arrH(iSrc, :, :))';
     
    T = h' * iR_N_iR * h;
    S = h' * iR * h;
    
    [V , Ev]    = eig(S, T); %% => S*V = T*V*Ev, bigger EV == bigger S
    Evals       = diag(Ev);
    [~ , idx]   = sort(Evals,'descend');
    
    u = V  ( : , idx ( 1 ) ); 
    u = real (u / norm ( u ));
    
    lstU3D(iSrc,:) = u';
    
    hu = h * u ;  
    
    iTu         = invSPD(hu' * iR_N_iR * hu);
    Su          = hu' * iR * hu;  
    SNR(iSrc)   = trace(Su*iTu);
    
    Pow(iSrc)   = trace(hu' * invSPD(sIn.Rm1) * hu); 
    Pnoise(iSrc)= trace(hu' * arrN * hu); 
    
end

% if only one source
nRef    = 0; 
order   = 1;

%% Then iterate taking one source out at a time

done = [];
for iIter = 2:nSrc
    
    nRef        = iIter -1;
    [~, idxMax] = sort(SNR(:), 'descend');  
    order       = idxMax(1:nRef);
    
    Hur = [];
    for o = 1:numel(order)
        hu_r    = squeeze(arrH(order(o), :, :))'*lstU3D(order(o),:)';
        Hur     = cat(2, Hur, hu_r);
    end
    
    iR_N_iR_Hur = iR_N_iR * Hur;
    iR_Hur      = iR      * Hur;
    
    Tr          = Hur' * iR_N_iR_Hur;
    Sr          = Hur' * iR_Hur;
    
    iTr         = invSPD(Tr);
   
    iTr_Sr_iTr  = iTr * Sr * iTr;
 
    for iSrc = 1:nSrc
        if any(iSrc == order) && ~any(done == iSrc)
            done      = cat(1, done, iSrc); 
            SNR(iSrc) = 999.9- iIter/10; % just a big number to keep order
            
        elseif iSrc ~= order 
            
            h   = squeeze(arrH(iSrc, :, :))';
 
            T   = h' * iR_N_iR * h;
            S   = h' * iR      * h;
            
            Tsr = h' * iR_N_iR_Hur;
            Ssr = h' * iR_Hur;
            
           
            D   = Tsr * iTr_Sr_iTr * Tsr' - Tsr * iTr * Ssr' - Ssr * iTr * Tsr' + S;
            F   = T - Tsr * iTr * Tsr';
            
            [V ,Ev]  = eig(D, F);
            Evals    = diag(Ev);
            [~, idx] = sort(Evals,'descend');
            
            u   = V(:, idx(1) ); 
            u   = real (u / norm(u));
           
            lstU3D(iSrc,:) = u';
            
            hu  = h * u ; 
            Hu  = horzcat(Hur, hu);
          
            iTu       = invSPD(Hu' * iR_N_iR * Hu); 
            Su        = Hu' * iR * Hu;  
            SNR(iSrc) = trace(iTu * Su); 
       

        end
    end
end


%% Compute the BF Weights
% All but the last are references
  
refU   = lstU3D(order,:);
notRef = setdiff(1:nSrc, order);

HrefT = [];
for iRef = 1:nRef
    hT      = refU(iRef, :) * squeeze(arrH(order(iRef),:,:));
    HrefT   = [HrefT; hT];    
end

if ~isempty(HrefT)
    hT = lstU3D(notRef,:) * squeeze(arrH(notRef,:,:)); 
    Ht = [HrefT; hT];
    wT = invSPD(Ht * sIn.Rm1 * Ht') * Ht * sIn.Rm1;
else
    % Single source case
    U = lstU3D;
    hT = U * squeeze(arrH); 
    wT = hT * sIn.Rm1 / (hT * sIn.Rm1 * hT');
end    

%% set outputs

sOut.U3D     =  lstU3D';
sOut.order   = [order; notRef]';
sOut.arrW(:,[order;notRef]) = wT(:,:)';

if ~Is_Evoked; SNR  = SNR -nSrc; end
 
sOut.beamSNR = SNR(notRef);

end

