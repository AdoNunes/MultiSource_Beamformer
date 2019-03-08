# MultiSource_Beamformer
Multiple Constrained Minimum Variance (MCMV) beamformer for connectivity analysis


here it is provided:

- MCMV_BF.m for estimating strongest orientation and computing beamformer weights.

- MCMV_test.m a script that simulates a few signals and estimates connectivity using the single-source scalar LCMV and the multi-source MCMV beamformers. 
    - LF_4Src.mat are the leadfields necessary for the the simulation signals. They are quite close, some within 2 cm, good for seeing the effects of signal leakage. The lower the SNR the wider will be the nulls imposed by the MCMV, and less optimal will be for close by sources.  
    
- Analysis_scripts, folder with the main analyses done in the preprint:
    - Multiple constrained minimum variance beamformer (MCMV) performance in connectivity analyses
      - https://www.biorxiv.org/content/10.1101/567768v1
