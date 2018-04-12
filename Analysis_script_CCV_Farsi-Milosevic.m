% The following is an example script for detecting, sampling, and calibrating fluorescence from images of purified synaptic vesicles or clathrin coated vesicles as used in Farsi et al 2016 and Farsi et al 2018. 
% This script provides code to:
% - load raw data (.tif stacks)
% - perform single SV or CCV detection using wavelet based, multi-scale spot detection (see Olivio-Marin, J.C. Pattern Recognition 35(9):1989-1996. Sept. 2002)
% - determine SV or CCV centroid and local background
% - sample background corrected fluorescence intensity
% - calibration of fluorescence to pH or membrane potential
%
% Copyright (c) 2018, Zohreh Farsi, Andrew Woehler, Berlin Institute for Medical Systems Biology
% All rights reserved.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%% loading .tif image stacks
clear
load('MyColormaps_jetpbII.mat');
[FileName,PathName] = uigetfile('*.tif','Select tiff file'); 
fname = [PathName, FileName]; 
info = imfinfo(fname);
num_images = numel(info);

h = waitbar(0,'Loading image stack. Please wait...');
for k = 1:num_images
    F(:,:,k) = imread(fname, k); 
    waitbar(k/num_images)
end
F = double(F);
close(h)
%% Spot detection and spot mask creation
clear Fi
Fi = mean(F(:,:,1:10),3);

[detectionResults, detectionMask]=spotDetector(Fi); % (from Olivo-Marin, Pattern Recognition, 2002)
%subsize = 8; % leave a boarder without any spots so that the later sub-frame images do not extend beyond the frame
s_mask = im2bw(detectionMask);
s_mask(1:subsize, :) = 0; s_mask(:, 1:subsize) = 0;
s_mask(size(F,1)-subsize:size(F,1),:) = 0;
s_mask(:,size(F,2)-subsize:size(F,2)) = 0;

s_label = bwlabel(s_mask);

% Get stats on the detected spots
cmin = min(F(:)); cmax = max(F(:));
s_stats  = regionprops(s_label, Fi, 'WeightedCentroid', 'Area','Eccentricity');
centroids = cat(1, s_stats.WeightedCentroid);
F_area = cat(1, s_stats.Area);
F_ecc = cat(1, s_stats.Eccentricity);
% Selecting the detected spots based on their eccentricity and size
ecc_th = 0.8; (0 = CIRCLE, 1 = LINE)
area_hth = 20; 
selected = ones(size(F_area));
selected(F_ecc > ecc_th) = NaN;
selected(F_area > area_hth) = NaN;
cent_sel = centroids(isfinite(selected),:);
F_area_sel = F_area(isfinite(selected));
F_ecc_sel = F_ecc(isfinite(selected));
%%
clear F_m
icent=round(cent_sel); % round the centroid to the nearest int. coordinate
F_m = mean(F(:,:,1:size(F,3)),3);
%% take the bgp no. pixels with the lowest value as the local background
for i = 1:length(cent_sel);
    clear subsize
    subsize = 2*ceil(sqrt(F_area_sel(i)/pi)+0); %need to get even number, can add # to get buffer
    temp_m = F_m(icent(i,2)-subsize/2:icent(i,2)+subsize/2,icent(i,1)-subsize/2:icent(i,1)+subsize/2);
    
    bgp(i,:) = ceil(0.20*length(temp_m(:))); 
    temp_bgpv = sort(temp_m(:));
    bgpv = temp_bgpv(bgp(i));
    size_sub(i,:) = length(temp_m(:));
    
    sub_bg_mask = zeros(size(temp_m));
    sub_bg_mask(temp_m <= bgpv) = 1;

    sub_s_mask = zeros(size(temp_m));
    sub_s_mask(temp_m > bgpv) = 1;  
    
    for j = 1:num_images
        temp = (F(icent(i,2)-subsize/2:icent(i,2)+subsize/2,icent(i,1)-subsize/2:icent(i,1)+subsize/2,j));
        temp_bgs = temp - mean(temp(sub_bg_mask==1));
        Fmean(i,j)= mean(temp_bgs(sub_s_mask==1));
        Fsum(i,j)= sum(temp_bgs(sub_s_mask==1));
        Fmax(i,j)= max(temp_bgs(sub_s_mask==1));
           end
end
%% Initial fluorescence before chemical perturbation
clear Fs_pre
Fs_pre = mean(Fsum(:,1:10),2); 
%% normalized fluorescence
clear Fsum_norm
for i=1:size(Fsum,1);
    for j = 1:size(Fsum,2);
    Fsum_norm(i,j)=Fsum(i,j)./Fs_pre(i);
    end
end
%% fluorescence to pH conversion for spH-CCVs and spH-SVs
clear pH_Fsum
pka = 7.2;
Fmax_pH = (1+10^(pka-7.4));
for i=1:size(Fsum_norm,1);
    for j = 1:size(Fsum_norm,2);
    pH_Fsum(i,j) = pka-log10((Fmax_pH-Fsum_norm(i,j))/Fsum_norm(i,j));
    end
end
%% Fluorescence to membrane potential conversion for VF2.1.Cl-labeled CCVs and SVs
deltaF=mean(Fsum_norm(:,40:50),2)-mean(Fsum_norm(:,1:10),2); % fuorescence-chang inducedby a chemical perturbation
mp=deltaF*370.37; % membrane potential (for VF2.1.Cl (Miller, et al. PNAS, 2011)