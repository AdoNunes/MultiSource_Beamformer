function sOut = MCMV_BF(sIn)
%
% SYNTAX:
%   sOut = MCMV_BF(sIn)
%
% This function generates a set of MCMV (multi-source) scalar beamformer
% weights for a given set of nSrc (= number of sources) forward solution
% triplets H_i, i = 1:nSrc. Each triplet H_i(r_i) = [h_x, h_y, h_z] is a (M
% x 3) matrix, where M is a number of sensors in a sensor array and h_? is
% a (M x 1) lead field vector of a dipole located at point r_i and oriented
% in x, y, z directions, respectively. The function produces (M x nSrc)
% matrix W of beamformer weights, and (3 x nSrc) matrix U of source
% orientations. These orientations are approximate and are found using
% iteration procedure described in an article Moiseev et al., Neuroimage,
% 2011, v.58, p.481?496.
% For some sources, orientations may be known in advance. Those sources are
% referred to as "reference sources" or "references" here. IT IS ASSUMED
% THAT IN THE LIST OF TRIPLETS H_i REFERENCES (IF ANY) ALWAYS GO FIRST. If
% all orientations are known, then beamformer weights are calculated using
% a simple algebraic expression and no iterations are needed.
%
% Technically, this function is an adapter to a more general
% iterateMultiSrcBeamA() function, for a simple case when source locations
% are fixed and known in advance.
%
% Input:
%   sIn         - a structure with input arguments, with the following
%                 fields:
%                 : arrH    (nSrc x 3 X M)) (3-dimensional array: "Vector"
%                           lead fields for the sources. For each target
%                           "source" i, the 3 x M leadfield triplet is =
%                           arrH(i,:,:) = [h_x, h_y, h_z]', i = 1,...,nSrc,
%                           h_? being a lead field of a dipole oriented in
%                           ?-direction. nRef sources with known
%                           orientations should be listed first; nRef may
%                           be 0
%                 : Rm1     (M x M) an inverse of the full (= noise +src)
%                           covariance matrix
%                 Optional:
%                 : sBeamType (string). Beamformemr type, one of MPZ, MAI,
%                           MER, rMER. These types are described in the paper
%                           referenced above. If not supplied, 'MPZ'
%                           (Multi-source Pseudo-Z) is assumed by default
%                 : arrN    (M x M) noise covariance. If not supplied, a
%                           diagonal noise will be used
%                 : Cavg    (M x M) Matrix of 2nd moments of the
%                           epoch-averaged sensor signals. MUST be supplied
%                           if evoked beamformer (MER or rMER) is used; not
%                           needed for power beamformerms (MPZ, MAI)
%                 : Uref    (3 x nRef) Orientations of reference sources.
%                           Each column of Uref should be a vector of unit
%                           length. nRef should be <= nSrc
% Output:
%   sOut        - a structure with the results, with the following fields:
%                 : U3D     (3 x nSrc) - the found orientations of the sources
%                 : arrW    (M x nSrc) - a set of MCMV weights
%                 : beamSNR (double) - the beamformer joint power "SNR"
%                           value for a set. Exact mathematical expression
%                           for this quantity depends on the beamformer
%                           type, but the physical meaning is the same
%
% A.Moiseev, BCNI, October 2017

[nSrc, ~, M] = size(sIn.arrH);



if ~isfield(sIn, 'arrN')
    % Noise covariance was not specified - generate a diagonal one
    arrN = M*eye(M)/trace(sIn.Rm1);   % This sets noise power ~ mean EV of R
else
   arrN = sIn.arrN;  
end

arrH = sIn.arrH;
 
if ~isfield(sIn, 'Cavg')
    iR =  sIn.Rm1;
    Is_Evoked = 0;
else  % Event related BF
    iR_C_iR = sIn.Rm1*sIn.Cavg*sIn.Rm1;
    iR =  iR_C_iR;
    Is_Evoked = 1;
end

iR_N_iR =  sIn.Rm1*arrN* sIn.Rm1;

for iSrc = 1:nSrc % get U3D for each source
    h = squeeze(arrH(iSrc, :, :))';
     
    
    T = h' * iR_N_iR * h;
    S = h' * iR * h;
    
    [V , Ev] = eig(S, T); %% => S*V = T*V*Ev, bigger EV == bigger S
    Evals = diag(Ev);
    [Ev , idx] = sort(Evals,'descend');
    
    u = V  ( : , idx ( 1 ) ); 
    u = real (u / norm ( u ));
    
    lstU3D(iSrc,:) = u';
    
    hu = h * u ;  
    
    iTu = invSPD(hu' * iR_N_iR * hu);  %
    Su = hu' * iR * hu;  % Full G-matrix, shouldn t be iR_N_iR?
    SNR(iSrc) = trace(Su*iTu);  % Pmai = tr(S* inv(T))-n, T = Y
    
    Pow(iSrc) = trace(hu' * invSPD(sIn.Rm1) * hu); 
    Pnoise(iSrc) = trace(hu' * arrN * hu); 
    
end


% if only one source
[SNR_max idxMax] = sort(SNR(:), 'descend'); 

nRef = 0; 
order = 1;

% else MCMV
Hur = [];
done = [];
for iIter = 2:nSrc
    nRef = iIter -1;
    
    [SNR_max idxMax] = sort(SNR(:), 'descend');  
    
    order = idxMax(1:nRef);
    
    Hur = [];
    for o = 1:numel(order)
        hu_r = squeeze(arrH(order(o), :, :))'*lstU3D(order(o),:)';
        Hur = cat(2, Hur, hu_r);
    end
    
    iR_N_iR_Hur = iR_N_iR * Hur;
    iR_Hur = iR * Hur;
    
    Tr = Hur' * iR_N_iR_Hur;
    Sr = Hur' * iR_Hur;
    
    iTr = invSPD(Tr);
   
    iTr_Sr_iTr = iTr * Sr * iTr;
 
    for iSrc = 1:nSrc
        if any(iSrc == order) && ~any(done == iSrc)
            done = cat(1, done, iSrc);
              
        SNR(iSrc) = 999.9- iIter/10;
            
        elseif iSrc ~= order 
            
            h = squeeze(arrH(iSrc, :, :))';
 
            T = h' * iR_N_iR * h;
            S = h' * iR      * h;
            
            Tsr = h' * iR_N_iR_Hur;
            Ssr = h' * iR_Hur;
            
            % Calculate matrices D and F
            D = Tsr * iTr_Sr_iTr * Tsr' - Tsr * iTr * Ssr' - Ssr * iTr * Tsr' + S;
            F = T - Tsr * iTr * Tsr';
            
            [V ,Ev] = eig(D, F);
            Evals = diag(Ev);
            [~, idx] = sort(Evals,'descend');
            
            u = V  ( : , idx ( 1 ) ); 
            u = real (u / norm ( u ));
           
            lstU3D(iSrc,:) = u';
            
            hu = h * u ; 
            Hu = horzcat(Hur, hu);
          
            iTu = invSPD(Hu' * iR_N_iR * Hu);  % Reconstructed source covariance
            Su = Hu' * iR * Hu;  
            SNR(iSrc) = trace(iTu*Su); 
           

        end
    end
end


%WEIGHTS
% takes the last one, previous as references
  
refU = lstU3D(order,:);

noRef = setdiff(1:nSrc, order);


HrefT = [];
for iRef = 1:nRef
    hT = refU(iRef, :) * squeeze(arrH(order(iRef),:,:));
    HrefT = [HrefT; hT];    
end

if ~isempty(HrefT)
    hT = lstU3D(noRef,:)*squeeze(arrH(noRef,:,:)); 
    Ht = [HrefT; hT];
    wT = invSPD(Ht * sIn.Rm1 *Ht') * Ht * sIn.Rm1;
else
    % Single source case
    U = lstU3D;
    hT = U*squeeze(arrH); 
    wT = hT * sIn.Rm1 / (hT * sIn.Rm1 * hT');
end    

sOut.U3D =  lstU3D';
sOut.order = [order; setdiff(1:nSrc, order)]';
sOut.arrW(:,[order;noRef]) = wT(:,:)';

 if ~Is_Evoked; SNR  = SNR -nSrc; end
 
sOut.beamSNR = SNR(idxMax(end));

end

