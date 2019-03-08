function W = run_APW_MCMV( arrH, arrR, pos, lambda, arrN, Pmat )
 
%
% Computes multisource beamformer Weights using an 'Augmented PairWise' 
% approach (APW-MCMV). For significant source pairs, it looks sources
% within 4 cm diameters that are significant. It includes these sources as
% nulls to improve signal leakage supression. 
%
% Input:
%  arrH: (nSrc x 3 X M sensors), the forward models for the sources.
%  arrR: (M x M) sensor covariance.
%  lambda: regularization parameters, for rank deficient arrR.
%  arrN: (M x M) noise covariance.
%  Pmat: (nSrc x nSrc) is a binary matrix indicating if a source pair is
%        significantly connected. Wherever Pmat(i,j)==1, APW-MCMV approach 
%        will be applied, otherwise simple PairWise MCMV will be used.
%
%  Returns the Weights for the source pair



if nargin<4
    lambda  = 0;
    arrN    = [];
end

if exist( 'arrN', 'var') && ~isempty(arrN)
    sIn.arrN = arrN;
else
    arrN = [];
end
 
if ~exist( 'Pmat', 'var')
    Pmat  = zeros(size(arrH,1));
else
    
n_source    = 2; % pairwise 
n_chann     = size(arrH,3);
source_pair = nchoosek(1:size(arrH,1), n_source); % make pair-wise list

W = zeros(size(source_pair,1),n_chann,n_source);

lambdamat = lambda * trace(arrR)/size(arrR,1);
sIn.Rm1 = inv(arrR + lambdamat * eye(size(arrR)));    % Invert the covariance

for ii = 1: size(source_pair,1)
    
    lstVoxNums = source_pair(ii,:); % Src pair
    
    if Pmat(lstVoxNums(1),lstVoxNums(2)) % if do multisource
        % get Src-pair pZ
        sIn.arrH =  arrH(lstVoxNums,:,:);
        sOut = constructMCMVWeights22sMER(sIn);
        Pz_ori = sOut.beamSNR;
        
        
        % Get distance Src pair with all other sources
        d1 = pdist2(pos(lstVoxNums(1),:),pos)';
        d2 = pdist2(pos(lstVoxNums(2),:),pos)';
        
        [d1_sort, d1_i] = sort(d1);
        [d2_sort, d2_i]=  sort(d2);
        
        % exclude current Src pair
        d1_sort(d1_i == lstVoxNums(1) |  d1_i == lstVoxNums(2)) = [];
        d2_sort(d2_i == lstVoxNums(1) |  d2_i == lstVoxNums(2)) = [];
        
        d1_i(d1_i == lstVoxNums(1) |  d1_i == lstVoxNums(2)) = [];
        d2_i(d2_i == lstVoxNums(1) |  d2_i == lstVoxNums(2)) = [];
        
        
        neigh_distance    = 30; % max distance neighbors
        neigh_max         = 2;  % N neighbors per pair
        
        % Sort neighbors
        neig_1 = d1_i(d1_sort<neigh_distance);
        neig_2 = d2_i(d2_sort<neigh_distance);
        
 
        clear nei_hub1 nei_hub2
        
        % get strongest neighbors
        nei_hub1(:,1)       = sum(Pmat(neig_1,:),2);
        
        for n = 1:numel(neig_1)
            sIn.arrH        = arrH([lstVoxNums neig_1(n)],:,:);
            sOut            = constructMCMVWeights22sMER(sIn);
            nei_hub1(n,2)   = sOut.beamSNR;
        end
        
        if ~any(nei_hub1(:,1) ) && ~isempty(nei_hub1)% if neigh have no sig edges
            if ~any((nei_hub1(:,2)-3)/(Pz_ori-2)>1.45) % but contributes +45%
                str_nei = (nei_hub1(:,2)-3)/(Pz_ori-2)>1.45;
                nei_hub1 = nei_hub1(str_nei,:);
                neig_1 = neig_1(str_nei);
            else
                nei_hub1=[];
            end
        end
        
        % sort them based on connections and distance
        if ~isempty(nei_hub1) 
            [~, i1] = sortrows(nei_hub1,[-1 -2]  );
            Nnei1   = min(numel(i1), neigh_max);
        else
            i1      =[];
            Nnei1	=[];
        end
        
        neig_1      = neig_1(i1(1:Nnei1));
        
        [~,~,ix2]   = intersect(neig_1,neig_2); % take out common ones
        neig_2(ix2) = [];
        
        % get strongest neighbors
        nei_hub2(:,1)       = sum(Pmat(neig_2,:),2);
        for n = 1:numel(neig_2)
                sIn.arrH        = arrH([lstVoxNums neig_2(n)],:,:);
                sOut            = constructMCMVWeights22sMER(sIn);
                nei_hub2(n,2)   = sOut.beamSNR;
        end
        
        if  ~any(nei_hub2(:,1) ) && ~isempty(nei_hub2) % if neigh have no sig edges
            if ~any((nei_hub2(:,2)-3)/(Pz_ori-2)>1.45) % but contributes +45%
                str_nei = (nei_hub2(:,2)-3)/(Pz_ori-2)>1.45;
                nei_hub2 = nei_hub2(str_nei,:);
                neig_2 = neig_2(str_nei);
            else
                nei_hub2=[];
            end
        end
        
        % sort them based on connections and distance
        if ~isempty(nei_hub2)
            [~, i2] = sortrows(nei_hub2,[-1 -2]  );
            nnei2   = min(numel(i2), neigh_max);
        else
            i2      =[];
            nnei2   =[];
        end
        
        neig_2 = neig_2(i2(1:nnei2));
        
        % get neighbours pZ and LF
        src         = [lstVoxNums,neig_1',neig_2'];
        sIn.arrH    = arrH(src,:,:);
        sOut        = constructMCMVWeights22sMER(sIn);
        
        W(ii,:,:)   = sOut.arrW(:,1:2);
       
        
    else % no multisource
        sIn.arrH    =   arrH(lstVoxNums,:,:);
        sOut        = constructMCMVWeights22sMER(sIn);
        W(ii,:,:)   = sOut.arrW(:,1:2);
    end
    
    
end

end
