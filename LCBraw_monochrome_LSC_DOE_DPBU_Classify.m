% LCBraw performs grayspot/blemish
% defects detection on flat field images
%
% INPUT:
%   IDraw: raw image data
%   bayerFormat: bayer format:
%       'bggr'
%       'rggb'
%       'grbg'
%       'gbrg'
%   ROISize:            size of roi in [m n] format
%   filterWidth:        the width of the filter.
%   threshold:          threshold limit for severe defects
%
% OUTPUT: See ERS
%
%  Author: Aaron Cong
%          March 2016
%  Algorithm: Cathy Zhang
%             May 2012
%
%******************************************************************
% ver C6.1-1: lens shading correction, removing cropping, apply corner
% masking based on field size, update faster binning (matlab)
% ver C6.2-1 2016-06-06 : add lens shading correction input, update
% preprocess method
% ver c3.0-1 2017-03-12: LCB DOE
%******************************************************************
function LCBraw_monochrome_LSC_DOE_DPBU_Classify(fid, index, filename, IDraw, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold,field_effective, option_lsc, ls_mat)
%clear all;
%close all;
%clc;

% parameters for getFile
%inputFormat = 'lg';
%outputFormat = 'raw10';
%bayerFormat = 'rggb';
%bitDepth = 10;
%pedestal = 0;

%roiSize = [7 7];
%HF:5 LF:9
%filterWidth = 9;you
%threshold = 16;
%field_effective = 0.93;
feature_roiSize = 50; % square size
%option_lsc = 3;
test_ver=3.01;
debug = 1; % from 1 to 0, Sam, 220411
%% select image files
% [filename, pathname] = uigetfile('*.raw','Please select the Shading images', 'MultiSelect', 'on');
% if iscell(filename)
%     FN = filename;
% else
%     FN{1} = filename;
% end
% numFiles = length(FN);
% %% calculation
% for j = 1:numFiles
%     Results{j}.fullFileName = [pathname,FN{j}];
%     Results{j}.fileName = FN{j};
%     % read in rawq data
%     IDraw = getFile([pathname,FN{j}], inputFormat, outputFormat);
%     IDraw = double(IDraw);
% end

[ls_mat] = lensShading_monochrome_noload_9x9(IDraw, bayerFormat, pedestal, bitDepth, []); % Sam, 220413
% [ls_mat] = lensShading_monochrome_noload_9x11(IDraw, bayerFormat, 0, 10, []);
% preprocess raw image
IDraw = preprocess(IDraw, bayerFormat, 'raw', bitDepth, pedestal, 83.5, false, true, option_lsc, ls_mat); % (READ POINT: IDraw)
[height width] = size(IDraw);

roiX_size = roiSize(2);
roiY_size = roiSize(1);

fw = filterWidth;

CornerROISize= [filterWidth+1,filterWidth+1];

cornerX_size = CornerROISize(2);
cornerY_size = CornerROISize(1);

% define the filter
h = [1/2, zeros(1, (fw-3)/2), -1, zeros(1, (fw-3)/2), 1/2];
%scale down input image
IDbin_all = binImage(IDraw, [roiY_size roiX_size]);
[rows cols cc] = size(IDbin_all);
output = []; %output array
for i=1:cc
    %select 1 channel
    IDbin = IDbin_all(:,:,i);
    IDbin_min = min(IDbin(:));
    IDbin_max = max(IDbin(:));
    IDbin_nor = (IDbin-IDbin_min) / (IDbin_max-IDbin_min);
    imgname = sprintf('%s01_%s_IDbin_normalize.png', index, filename);
    imwrite(IDbin_nor, imgname);
    
    %% ABC's flow, generate filtered image
    [I_filtered_c, I_filtered_c_cen, I_filtered_c_edg, I_filtered_c_cor, ...
     LCBraw_center, LCBraw_edge, LCBraw_corner] = ...
        ABC_Filter(IDbin, roiSize, fw, CornerROISize, h, threshold, rows, cols, field_effective);
    figure(9), subplot(2,5,1), imshow(IDbin, []), title('IDbin');
    figure(9), subplot(2,5,4), imshow(I_filtered_c_cen, []), title({'filtered_cen,';'IDbin'},'interpreter','none');
    figure(9), subplot(2,5,5), imshow(I_filtered_c_edg, []), title({'filtered_edg,';'IDbin'},'interpreter','none');    
    I_filtered_c_min = min(I_filtered_c(:));
    I_filtered_c_max = max(I_filtered_c(:));
    I_filtered_c_nor = (I_filtered_c-I_filtered_c_min) ./ (I_filtered_c_max-I_filtered_c_min);
    imgname = sprintf('%s02_%s_I_filtered_c_normalize.png', index, filename);
    imwrite(I_filtered_c_nor, imgname);
    
    %% ABC's flow, locate coordinate of maximum value, center
    [LCBraw_center_x, LCBraw_center_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        ABC_Locator(I_filtered_c_cen, roiSize, feature_roiSize, fw, rows, cols);

    %% ABC's flow, calculate particle statistics within ROI, center
    [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
        xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
    output_particle_cell = num2cell(output_particle.stats);
    [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
    figure(9), subplot(2,5,6), imshow(IDbin_roi), title({'roiBin_cen,';'enhanced IDbin'},'interpreter','none');
    imgname = sprintf('%s04_%s_IDbin_roi_center_%.4f.png', index, filename, particle_attenuation_raw);
    imwrite(IDbin_roi, imgname);
    
    %% ABC's flow, locate coordinate of maximum value, corner
    [LCBraw_corner_x, LCBraw_corner_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        ABC_Locator(I_filtered_c_cor, roiSize, feature_roiSize, fw, rows, cols);
    
    %% ABC's flow, locate coordinate of maximum value, edge
    [LCBraw_edge_x, LCBraw_edge_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        ABC_Locator(I_filtered_c_edg, roiSize, feature_roiSize, fw, rows, cols);

    %% ABC's flow, calculate particle statistics within ROI, edge
    if particle_intensity == 0
        [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
            xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
        output_particle_cell = num2cell(output_particle.stats);
        [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
        figure(9), subplot(2,5,7), imshow(IDbin_roi), title({'roiBin_edg,';'enhanced IDbin'},'interpreter','none');
        imgname = sprintf('%s05_%s_IDbin_roi_edge_%.4f.png', index, filename, particle_attenuation_raw);
        imwrite(IDbin_roi, imgname);
    end
    if debug
        figure(19);imshow(IDraw, []);
        truesize([uint16(height/2) uint16(width/2)]);
        hold on;
        if particle_area_x > 0 && particle_area_y > 0 && particle_h_raw > 0 && particle_w_raw > 0
        rectangle('position', [particle_area_x, particle_area_y, particle_h_raw, particle_w_raw], 'EdgeColor', 'r');
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 - 15), ['w: ' num2str(particle_w_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 + 15), ['h: ' num2str(particle_h_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        end
        scatter([LCBraw_center_x LCBraw_edge_x LCBraw_corner_x],[LCBraw_center_y LCBraw_edge_y LCBraw_corner_y],'yo');
        hold off;
        imgname = sprintf('%s08_%s_Result.png', index, filename);
%         saveas(figure(21), imgname);
    end

    
    
    %% DPBU LCBraw V2.0
    %% enhance contrast, Sam, 220428
    IDbin_min_ratio = 6.5;
    IDbin_max_ratio = 6.5;
    IDbin_avg = mean(IDbin(:));
    IDbin_std = std(IDbin(:));
    IDbin_min = IDbin_avg - IDbin_min_ratio*IDbin_std;
    IDbin_max = IDbin_avg + IDbin_max_ratio*IDbin_std;
    IDbin_hist = (IDbin-IDbin_min) ./ (IDbin_max-IDbin_min);
    figure(10), subplot(2,5,1), imshow(IDbin_hist, []), title({'enhanced IDbin,';'extend min/max'});
    imgname = sprintf('%s11_%s_IDbin_hist_normalize.png', index, filename);
    imwrite(IDbin_hist, imgname);
%     figure(10), subplot(1,2,1),
%     s = surf(IDbin_hist);
%     s.EdgeColor = 'none';



    % transfer to frequency domain, Sam, 220606
    H1 = fft2(double(IDbin));
    H2 = fftshift(H1);
    T1 = log(abs(H2));
    figure(100), subplot(2, 2, 3), imshow(T1, []); title('frequency domain');
        
    S1 = H2;
    [M,N] = size(S1);
    n1 = floor(M/2);
    n2 = floor(N/2);
    d0 = 5; % cutoff is the cutoff frequency of the filter 0 - 0.5
    n = 10; % n is the order of the filter, the higher n is the sharper the transition is. (n must be an integer >= 1).
    for i = 1:M     
        for j = 1:N   
            d = sqrt((i-n1)^2+(j-n2)^2);      
%             h = exp(-1/2*(d^2/d0^2)); % gaussian low-pass filter
            h = 1-exp(-1/2*(d^2/d0^2)); % gaussian high-pass filter
%             h = 1/(1 + (d/d0)^(2*n)); % butterworth low-pass filter
%             h = 1/(1 + (d0/d)^(2*n)); % butterworth high-pass filter
            S1(i,j) = h*S1(i,j);
        end
    end
    T2 = log(abs(S1));
    S2 = ifftshift(S1); 
    S3 = real(ifft2(S2));
%     S3_avg = mean(S3(:));
%     S3_std = std(S3(:));
%     S3_bias = 3;
%     S3_ulm = S3_avg + S3_bias*S3_std;
%     S3_dlm = S3_avg - S3_bias*S3_std;
%     S3 = (S3 - S3_dlm) .* (S3_ulm-S3_dlm);
%     S3(S3<0) = 0;
%     S3(S3>1) = 1;
    S3_min = min(S3(:));
    S3_max = max(S3(:));
    S3_nor = (S3-S3_min) ./ (S3_max-S3_min);
    figure(100), subplot(2, 2, 4), imshow(T2, []), title({'high-pass','frequency domain'});
    figure(100), subplot(2, 2, 1), imshow(IDbin_hist,[]), title('origin image');
    figure(100), subplot(2, 2, 2), imshow(S3_nor,[]), title('filterd image');
%     figure(101), imshow(IDbin_hist, []), title('origin image');
%     figure(102), imshow(S3, []), title('filterd image');

    IDbin_hist = S3_nor;
    %% sobel edge filter, Sam, 221107
    mx = [-1, 0, 1;
          -2, 0, 2;
          -1, 0, 1]; % Sobel Gx kernel
    my = mx'; % gradient Gy
    IDbin_hist_smooth = imgaussfilt(IDbin_hist, 'FilterSize', 7, 'Padding', 'symmetric');
    sobel_x = imfilter(double(IDbin_hist_smooth), mx, 'replicate');
    sobel_y = imfilter(double(IDbin_hist_smooth), my, 'replicate');
    sobel = sqrt(sobel_x.^2 + sobel_y.^2); % Find magnitude
    sobel_min = min(sobel(:));
    sobel_max = max(sobel(:));
    sobel_nor = (sobel-sobel_min) ./ (sobel_max-sobel_min);
%     sobel_16bit = uint16(sobel ./ max(sobel(:)) .* 2^16);
%     figure(10), subplot(2,5,2), title('sobel filter'), imshow(sobel_nor,[]);
    imgname = sprintf('%s12_%s_sobel_normalize.png', index, filename);
%     imwrite(sobel_8bit, imgname);

    %% create center, edge, and corner masks
    corner_size = (cornerX_size+cornerY_size) / 2;
    ULm = false(corner_size, corner_size);
    for y=1:corner_size
        for x=1:corner_size-y+1
            ULm(y,x) = 1;
        end
    end
    URm = fliplr(ULm);
    LLm = flipud(ULm);
    LRm = fliplr(LLm);
    UL = [1:corner_size;           1:corner_size];
    UR = [1:corner_size;           cols-corner_size+1:cols];
    LL = [rows-corner_size+1:rows; 1:corner_size];
    LR = [rows-corner_size+1:rows; cols-corner_size+1:cols];
    mask = false(size(IDbin));
    mask(UL(1,:), UL(2,:)) = ULm;
    mask(UR(1,:), UR(2,:)) = URm;
    mask(LL(1,:), LL(2,:)) = LLm;
    mask(LR(1,:), LR(2,:)) = LRm;
    maskULi = false(size(IDbin)); maskUL = false(size(IDbin));
    maskURi = false(size(IDbin)); maskUR = false(size(IDbin));
    maskLLi = false(size(IDbin)); maskLL = false(size(IDbin));
    maskLRi = false(size(IDbin)); maskLR = false(size(IDbin));
    maskULi(UL(1,:), UL(2,:)) = ~ULm; maskUL(UL(1,:), UL(2,:)) = ULm;
    maskURi(UR(1,:), UR(2,:)) = ~URm; maskUR(UR(1,:), UR(2,:)) = URm;
    maskLLi(LL(1,:), LL(2,:)) = ~LLm; maskLL(LL(1,:), LL(2,:)) = LLm;
    maskLRi(LR(1,:), LR(2,:)) = ~LRm; maskLR(LR(1,:), LR(2,:)) = LRm;
    
    %% create enhanced corner mask, Sam, 221123
    corner_size = (cornerX_size+cornerY_size) / 2;
    corner_size = corner_size * 1.5;
    ULm = false(corner_size, corner_size);
    for y=1:corner_size
        for x=1:corner_size-y+1
            ULm(y,x) = 1;
        end
    end
    URm = fliplr(ULm);
    LLm = flipud(ULm);
    LRm = fliplr(LLm);
    UL = [1:corner_size;           1:corner_size];
    UR = [1:corner_size;           cols-corner_size+1:cols];
    LL = [rows-corner_size+1:rows; 1:corner_size];
    LR = [rows-corner_size+1:rows; cols-corner_size+1:cols];
    mask_cor = false(size(IDbin));
    mask_cor(UL(1,:), UL(2,:)) = ULm;
    mask_cor(UR(1,:), UR(2,:)) = URm;
    mask_cor(LL(1,:), LL(2,:)) = LLm;
    mask_cor(LR(1,:), LR(2,:)) = LRm;
    
%     % multiply gain for corner, Sam, 221122
%     sobel_cor_gain = 1.25;
%     sobel_enh_gain(mask_cor) = sobel_enh_gain(mask_cor) .* sobel_cor_gain;
    
        %% enhance edge by sobel filter, Sam, 221107
%     sobel_smooth = imgaussfilt(sobel, 'FilterSize', 7, 'Padding', 'symmetric');
    sobel_max = max(sobel_nor(:));
%     sobel_smooth_noise_th = 0.0825; % precent of maximum pixel value
%     sobel_smooth_noise_th = 0.165; % precent of maximum pixel value
%     sobel_smooth_noise_th = 0; % precent of maximum pixel value
%     sobel_smooth(sobel_smooth < sobel_smooth_noise_th*sobel_smooth_max) = 0;
    sobel_smooth_noise_th   = [0 0.10 0.20 0.35 0.65 0.85 1.00];
    sobel_smooth_noise_gain = [0 0.15 0.15 0.35 0.70 1.00 1.00];
    sobel_smooth = zeros(size(sobel_nor));
    segment = length(sobel_smooth_noise_th) - 1;
    for idx=2:length(sobel_smooth_noise_th)
        th0 = sobel_smooth_noise_th(idx-1)*sobel_max;
        th1 = sobel_smooth_noise_th(idx)*sobel_max;
        ga = sobel_smooth_noise_gain(idx);
        mask = sobel_nor > th0 & sobel_nor<= th1;
        mask_cnt = sum(mask(:));
%         sobel_smooth(mask) = sobel_nor(mask) .* ga;
        % for debug used, Sam, 221123
        sobel_lv = sobel_nor;
        sobel_lv(~mask) = 0;
        sobel_lv_smooth = imgaussfilt(sobel_lv, 'FilterSize', 5, 'Padding', 'symmetric');
%         sobel_lv_smooth_mask = sobel_lv_smooth(sobel_lv_smooth>0);
%         sobel_smooth = sobel_smooth + sobel_lv_smooth;
        sobel_smooth(mask) = sobel_lv_smooth(mask) .* ga;
        title_1 = sprintf('range:%.2f-%.2f', th0, th1);
        title_2 = sprintf('gain:%.2f', ga);
        figure(31), subplot(2,segment,idx-1), imshow(sobel_lv_smooth), title({title_1,title_2});
        figure(31), subplot(2,segment,segment+idx-1), imshow(sobel_lv_smooth .* ga,[]), title({title_1,title_2});
    end
%     sobel_bias = 0.5;
    sobel_bias = 0.0;
%     sobel_gain = 1.25;
%     sobel_gain = 1.5;
    sobel_gain = 0.35;
    sobel_cor_gain = 0.55;
%     sobel_enh_gain = (sobel_smooth + sobel_bias) .* sobel_gain;
    sobel_enh_gain = sobel_smooth + sobel_bias;
    sobel_enh_gain(~mask_cor) = sobel_enh_gain(~mask_cor) .* sobel_gain;
    sobel_enh_gain(mask_cor) = sobel_enh_gain(mask_cor) .* sobel_cor_gain;
    sobel_enh_gain = imgaussfilt(sobel_enh_gain, 'FilterSize', 5, 'Padding', 'symmetric');
%     sobel_enh_gain = ((sobel_smooth ./ sobel_smooth_max) + sobel_bias) .* sobel_gain;
%     sobel_enh_gain = imgaussfilt(sobel_enh_gain, 'FilterSize', 7, 'Padding', 'symmetric');
    sobel_enh_gain_min = min(sobel_enh_gain(:));
    sobel_enh_gain_max = max(sobel_enh_gain(:));
    sobel_enh_gain_nor = (sobel_enh_gain-sobel_enh_gain_min) ./ (sobel_enh_gain_max-sobel_enh_gain_min);
%     sobel_enh_gain_nor = sobel_enh_gain;
%     figure(10), subplot(2,5,2), imshow(sobel_enh_gain_nor,[]), title('sobel filter');
    imgname = sprintf('%s13_%s_sobel_enhance_gain_normalize.png', index, filename);
    imwrite(sobel_enh_gain_nor, imgname);
    
    %% remove 4 corner, only for display, Sam, 221122
    sobel_enh_gain_rm4cor = sobel_enh_gain_nor .* ~mask_cor;
%     sobel_enh_gain_rm4cor(UL(1,:), UL(2,:)) = sobel_enh_gain_rm4cor(UL(1,:), UL(2,:)) .* UL_mask;
%     sobel_enh_gain_rm4cor(UR(1,:), UR(2,:)) = sobel_enh_gain_rm4cor(UR(1,:), UR(2,:)) .* UR_mask;
%     sobel_enh_gain_rm4cor(LL(1,:), LL(2,:)) = sobel_enh_gain_rm4cor(LL(1,:), LL(2,:)) .* LL_mask;
%     sobel_enh_gain_rm4cor(LR(1,:), LR(2,:)) = sobel_enh_gain_rm4cor(LR(1,:), LR(2,:)) .* LR_mask;
    figure(10), subplot(2,5,2), imshow(sobel_enh_gain_rm4cor,[]), title({'sobel filter as', 'gain map'});
    imgname = sprintf('%s13_%s_sobel_enhance_gain_rm4cor_normalize.png', index, filename);
%     imwrite(sobel_enh_gain_rm4cor_8bit, imgname);

    %% enhanced IDbin through edge filter gain, Sam, 221122
%     IDbin_hist_sobel_enh = zeros(size(IDbin));
%     IDbin_hist_sobel_enh(~mask) = IDbin_hist(~mask) - sobel_enh_gain(~mask);
%     IDbin_hist_sobel_enh(mask) = IDbin_hist(mask) + sobel_enh_gain(mask);
    IDbin_hist_sobel_enh = IDbin_hist - sobel_enh_gain;
    IDbin_hist_sobel_enh(IDbin_hist_sobel_enh<0) = 0; % Sam, 221123
%     IDbin_hist_sobel_enh(maskUL) = mean(IDbin_hist_sobel_enh(maskULi));
%     IDbin_hist_sobel_enh(maskUR) = mean(IDbin_hist_sobel_enh(maskURi));
%     IDbin_hist_sobel_enh(maskLL) = mean(IDbin_hist_sobel_enh(maskLLi));
%     IDbin_hist_sobel_enh(maskLR) = mean(IDbin_hist_sobel_enh(maskLRi));
    
    IDbin_hist_sobel_enh_min = min(IDbin_hist_sobel_enh(:));
    IDbin_hist_sobel_enh_max = max(IDbin_hist_sobel_enh(:));
    IDbin_hist_sobel_enh_nor = (IDbin_hist_sobel_enh-IDbin_hist_sobel_enh_min) ./ (IDbin_hist_sobel_enh_max-IDbin_hist_sobel_enh_min);
    figure(10), subplot(2,5,3), imshow(IDbin_hist_sobel_enh,[]), title({'enhanced IDbin,';'edge filter gain'});
    imgname = sprintf('%s14_%s_IDbin_hist_sobel_enh_normalize.png', index, filename);
    imwrite(IDbin_hist_sobel_enh_nor, imgname);

    %% ABC's flow, generate filtered image
    [I_filtered_c, I_filtered_c_cen, I_filtered_c_edg, I_filtered_c_cor, ...
     LCBraw_center, LCBraw_edge, LCBraw_corner] = ...
        ABC_Filter(IDbin_hist_sobel_enh, roiSize, fw, CornerROISize, h, threshold, rows, cols, field_effective);
    figure(10), subplot(2,5,4), imshow(I_filtered_c_cen, []), title({'filtered_cen,';'IDbin'},'interpreter','none');
%     figure(10), subplot(2,5,5), imshow(I_filtered_c_edg, []), title({'filtered_edg,';'IDbin'},'interpreter','none');    
    figure(10), subplot(2,5,5), imshow(I_filtered_c, []), title({'filtered_edg,';'IDbin'},'interpreter','none');    
    I_filtered_c_min = min(I_filtered_c(:));
    I_filtered_c_max = max(I_filtered_c(:));
    I_filtered_c_nor = (I_filtered_c-I_filtered_c_min) ./ (I_filtered_c_max-I_filtered_c_min);
    imgname = sprintf('%s15_%s_I_filtered_c_enhanced_normalize.png', index, filename);
    imwrite(I_filtered_c_nor, imgname);

    %% ABC's flow, locate coordinate of maximum value, center
    [LCBraw_center_x, LCBraw_center_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        DP_Locator(11, filename, 0, I_filtered_c_cen, roiSize, feature_roiSize, fw, rows, cols);

    %% DP's flow, calculate particle statistics within ROI, center
    [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
        xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
    output_particle_cell = num2cell(output_particle.stats);
    [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
    figure(10), subplot(2,5,6), imshow(IDbin_roi), title({'roiBin_cen,';'enhanced IDbin'},'interpreter','none');
    imgname = sprintf('%s16_%s_IDbin_enhanced_roi_center_%.4f.png', index, filename, particle_attenuation_raw);
%     imwrite(IDbin_roi, imgname);
    
    %% ABC's flow, locate coordinate of maximum value, corner
%     [LCBraw_corner_x, LCBraw_corner_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
%         DP_Locator(0, filename, 2, I_filtered_c_cor, roiSize, feature_roiSize, fw, rows, cols);
    [LCBraw_corner_x, LCBraw_corner_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        ABC_Locator(I_filtered_c_cor, roiSize, feature_roiSize, fw, rows, cols);
    
    %% DP's flow, locate coordinate of maximum value, edge
    [LCBraw_edge_x, LCBraw_edge_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        DP_Locator(12, filename, 1, I_filtered_c_edg, roiSize, feature_roiSize, fw, rows, cols);

    %% ABC's flow, calculate particle statistics within ROI, edge
    if particle_intensity == 0
        [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
            xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
        output_particle_cell = num2cell(output_particle.stats);
        [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
        figure(10), subplot(2,5,7), imshow(IDbin_roi), title({'roiBin_edg,';'enhanced IDbin'},'interpreter','none');
        imgname = sprintf('%s17_%s_IDbin_enhanced_roi_edge_%.4f.png', index, filename, particle_attenuation_raw);
%         imwrite(IDbin_roi, imgname);
    end
    if debug
        figure(15);imshow(IDraw, []);
        truesize([uint16(height/2) uint16(width/2)]);
        hold on;
        if particle_area_x > 0 && particle_area_y > 0 && particle_h_raw > 0 && particle_w_raw > 0
        rectangle('position', [particle_area_x, particle_area_y, particle_h_raw, particle_w_raw], 'EdgeColor', 'r');
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 - 15), ['w: ' num2str(particle_w_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 + 15), ['h: ' num2str(particle_h_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        end
        scatter([LCBraw_center_x LCBraw_edge_x LCBraw_corner_x],[LCBraw_center_y LCBraw_edge_y LCBraw_corner_y],'yo');
        hold off;
        imgname = sprintf('%s18_%s_Result.png', index, filename);
%         saveas(figure(21), imgname);
    end

    1==1;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% binning averages pixels in an m x n block into a reduced image size
% if image size is not a multiple of mxn, additional rows and
% columns are added to edge regions
%
% INPUT:
%   im: input image
%   blksz: [size_y size_x] size for each block
%
% OUTPUT:
%   RGB: RGB image
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imbin = binImage ( im, blksz)

% check if blksz has 2 elements
if length(blksz)~=2
    msg = 'Block size should be a 2 element vector.';
    error('MATLAB:myCode:input', msg);
end

% Image size
[H,W,c] = size(im);

% bin image with ROIsz blocks
[Hr,Wr] = deal(floor(H/blksz(1))*blksz(1), floor(W/blksz(2))*blksz(2));
borders = [ round(0.5*(H-Hr)), round(0.5*(W-Wr)) ];
imbin = zeros( Hr/blksz(1),Wr/blksz(2),c);
for i = 1:c
    
    % Average filter (faster than using For loop)
    meanIm = imfilter( im(borders(1)+1:borders(1)+Hr,borders(2)+1:borders(2)+Wr,i),...
        fspecial('average', blksz), 'symmetric','same' );
    imbin(:,:,i) = meanIm( ceil(blksz(1)/2):blksz(1):Hr, ceil(blksz(2)/2):blksz(2):Wr );
    
    % Special case for each borders
    xe = [1, borders(2)+blksz(2)+1:blksz(2):W-borders(2)-1, W+1];
    ye = [1, borders(1)+blksz(1)+1:blksz(1):H-borders(1)-1, H+1];
    
    % First and last column
    for j = 1:length(imbin(:,1,i))
        imbin(j,1,i) = mean2( im(ye(j):ye(j+1)-1, xe(1):xe(2)-1,i) );
        imbin(j,end,i) = mean2( im(ye(j):ye(j+1)-1, xe(end-1):xe(end)-1,i) );
    end
    
    % First and last row
    for j = 1:length(imbin(1,:,i))
        imbin(1,j,i) = mean2( im(ye(1):ye(2)-1, xe(j):xe(j+1)-1,i) );
        imbin(end,j,i) = mean2( im(ye(end-1):ye(end)-1, xe(j):xe(j+1)-1,i) );
    end
end
end

function IDbin_hist_rm4cor = Rm4Cor(IDbin_hist, cols, rows)
    num = 25;
    ext = 13;
    
    % top-left corner, Sam, 220511
    IDbin_hist_LT = IDbin_hist(1:num, 1:num); % top-left
    for k = 1:num
        non_zero = 0;
        for q = 1:k
            i = k-q+1;  j = q;
            if IDbin_hist_LT(j, i) < 0
                non_zero = 1;
            end
            IDbin_hist_LT(j, i) = 0;
        end
        % remove one more slash line, Sam, 220511
        if non_zero == 0
            for p = 1:ext
                for q = 1:(k+p)
                    i = (k+p)-q+1; j = q;
                    IDbin_hist_LT(j, i) = 0;
                end
            end
            break;
        end
    end
%     IDbin_hist_LT(IDbin_hist_LT>0) = 1; % Sam, 220525
%     IDbin_hist_LT = ~IDbin_hist_LT;
    
    % top-right corner, Sam, 220511
    IDbin_hist_RT = IDbin_hist(1:num, (cols-num+1):cols); % top-left
    for k = 1:num
        non_zero = 0;
        for q = 1:k
            i = (num-k+1)+q-1; j = q;
            if IDbin_hist_RT(j, i) < 0
                non_zero = 1;
            end
            IDbin_hist_RT(j, i) = 0;
        end
        % remove one more slash line, Sam, 220511
        if non_zero == 0
            for p = 1:ext
                for q = 1:(k+p)
                    i = (num-(k+p)+1)+q-1; j = q;
%                     fprintf('%d, %d\r', j, i);
                    IDbin_hist_RT(j, i) = 0;
                end
            end
            break;
        end
    end
%     IDbin_hist_RT(IDbin_hist_RT>0) = 1; % Sam, 220525
%     IDbin_hist_RT = ~IDbin_hist_RT;
    
    % bottom-left
    IDbin_hist_LB = IDbin_hist((rows-num+1):rows, 1:num); % top-left
    for k = 1:num
        non_zero = 0;
        for q = 1:k
            i = q; j = (num-k+1)+q-1;
%            fprintf('%d, %d\r', j, i);
            if IDbin_hist_LB(j, i) < 0
                non_zero = 1;
            end
            IDbin_hist_LB(j, i) = 0;
        end
        % remove one more slash line, Sam, 220511
        if non_zero == 0
            for p = 1:ext
                for q = 1:(k+p)
                    i = q; j = (num-(k+p)+1)+q-1;
%                     fprintf('%d, %d\r', j, i);
                    IDbin_hist_LB(j, i) = 0;
                end
            end
            break;
        end
    end
%     IDbin_hist_LB(IDbin_hist_LB>0) = 1; % Sam, 220525
%     IDbin_hist_LB = ~IDbin_hist_LB;
    
    % bottom-right
    IDbin_hist_RB = IDbin_hist((rows-num+1):rows, (cols-num+1):cols); % top-left
    for k = 1:num
        non_zero = 0;
        for q = 1:k
            i = (num-k+1)+q-1; j = num-q+1;
%             fprintf('%d, %d\r', j, i);
            if IDbin_hist_RB(j, i) < 0
                non_zero = 1;
            end
            IDbin_hist_RB(j, i) = 0;
        end
        % remove one more slash line, Sam, 220511
        if non_zero == 0
            for p = 1:ext
                for q = 1:(k+p)
                    i = (num-(k+p)+1)+q-1; j = num-q+1;
%                     fprintf('%d, %d\r', j, i);
                    IDbin_hist_RB(j, i) = 0;
                end
            end
            break;
        end
    end
%     IDbin_hist_RB(IDbin_hist_RB>0) = 1; % Sam, 220525
%     IDbin_hist_RB = ~IDbin_hist_RB;

%     cor_size_ext = 4;
%     corX_size_ext = cornerX_size + cor_size_ext;
%     corY_size_ext = cornerY_size + cor_size_ext;    
%     IDbin_hist_rm4cor(1:corY_size_ext, 1:corX_size_ext) = 0; % top-left
%     IDbin_hist_rm4cor(1:corY_size_ext, (cols-corX_size_ext+1):cols) = 0; % top-right
%     IDbin_hist_rm4cor((rows-corY_size_ext+1):rows, 1:corX_size_ext) = 0; % bottom-left
%     IDbin_hist_rm4cor((rows-corY_size_ext+1):rows, (cols-corX_size_ext+1):cols) = 0; % bottom-right

    IDbin_hist_rm4cor = IDbin_hist;
%     IDbin_hist_rm4cor = zeros(rows, cols); % Sam, 220525
    IDbin_hist_rm4cor(1:num, 1:num) = IDbin_hist_LT;
    IDbin_hist_rm4cor(1:num, (cols-num+1):cols) = IDbin_hist_RT;
    IDbin_hist_rm4cor((rows-num+1):rows, 1:num) = IDbin_hist_LB;
    IDbin_hist_rm4cor((rows-num+1):rows, (cols-num+1):cols) = IDbin_hist_RB;

%     figure(15);
%     subplot(2,2,1), imshow(IDbin_hist_LT), title('Left-Top');
%     subplot(2,2,2), imshow(IDbin_hist_RT), title('Right-Top');
%     subplot(2,2,3), imshow(IDbin_hist_LB), title('Left-Bottom');
%     subplot(2,2,4), imshow(IDbin_hist_RB), title('Right-Bottom');
%     imshow(IDbin_hist_rm4cor,[]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% thresholdROI
%   generates the ideal dot locations based on chart parameters
%
% INPUT:
%   IDy - ROI to be thresholded
%   bitDepth - BitDepth of the image
%   BorW - 1 for finding white dot, 0 for finding black dot
%
% OUTPUT:
%   OTSU_roi_bin - binary thresholded ROI
%
%   Author: Christos Hristovski
%           September 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [OTSU_roi_bin] = thresholdROI(IDy, bitDepth, BorW)
debug = 0;
%% B/W image based on Otsu's thresholding method
% Normalize image to 256 bins using given bitDepth
IDy_grey = uint8(IDy./(2^(bitDepth-8)));
counts = imhist(IDy_grey,256);
p = counts / sum(counts);
omega1 = 0;
omega2 = 1;
mu1 = 0;
mu2 = mean(IDy_grey(:));
for t = 1:256
    omega1(t) = sum(p(1:t));
    omega2(t) = sum(p(t+1:end));
    mu1(t) = sum(p(1:t).*(1:t)');
    mu2(t) = sum(p(t+1:end).*(t+1:256)');
end
sigma_b_squared_otsu = max(0,(mu1(end) .* omega1-mu1) .^2 ./(omega1 .* (1-omega1)));
[~,thres_level_otsu] = max(sigma_b_squared_otsu);

% Threshold using a buffered threshold that preserves more detail that
% normal by biasing towards the variance in the greys
if (BorW)
    OTSU_roi_bin = zeros(size(IDy_grey));
    thres_level_otsu = round(thres_level_otsu +  sqrt((var(sigma_b_squared_otsu)))/40);
    OTSU_roi_bin(IDy_grey <= thres_level_otsu) = 1;
else
    OTSU_roi_bin = ones(size(IDy_grey));
    thres_level_otsu = round(thres_level_otsu -  sqrt((var(sigma_b_squared_otsu)))/40);
    OTSU_roi_bin(IDy_grey >= thres_level_otsu) = 0;
end
if debug
    figure; imshow(IDy,[]);
    disp('thres_level_otsu, sqrt((var(sigma_b_squared_otsu)))/64');
    [thres_level_otsu, sqrt((var(sigma_b_squared_otsu)))/64]
    figure; imshow(OTSU_roi_bin,[]);
end
end

% LCBraw_center = max(I_filtered_c_cen(:));
% LCBraw_corner = max(I_filtered_c_cor(:));
% LCBraw_edge   = max(I_filtered_c_edg(:));
function [I_filtered_c, I_filtered_c_cen, I_filtered_c_edg, I_filtered_c_cor, ...
          LCBraw_center, LCBraw_edge, LCBraw_corner] = ...
    ABC_Filter(IDbin, roiSize, fw, CornerROISize, h, threshold, rows, cols, field_effective)
    roiX_size = roiSize(2);
    roiY_size = roiSize(1);

    cornerX_size = CornerROISize(2);
    cornerY_size = CornerROISize(1);

    %% horizontal direction
    % padding borders via extra interpolation before filtering
    IDbin_pad_h = zeros(size(IDbin,1), size(IDbin,2)+fw-1);
    for j = 1:size(IDbin,1)
        IDbin_pad_h(j,1:(fw+1)/2)       = interp1(0:(fw+1)/2-1, IDbin(j, 1:(fw+1)/2), -(fw-1)/2:0, 'linear', 'extrap');
        IDbin_pad_h(j,end-(fw-1)/2:end) = interp1(-(fw-1)/2:0, IDbin(j, end-(fw-1)/2:end), 0:(fw+1)/2-1, 'linear', 'extrap');
    end
    IDbin_pad_h(:,(fw+1)/2:end-(fw-1)/2) = IDbin;
    
    % 3x1 median filter (vertical)
    I_medfilt_h = medfilt2(IDbin_pad_h,[3 1]);
    
    % horizontal filtering
    I_filtered_pad_h = imfilter(I_medfilt_h, h, 'replicate');
    % crop back to original size
    I_filtered_h = I_filtered_pad_h(:,(fw+1)/2:end-(fw-1)/2);
    
    I_filtered_h(I_filtered_h<0) = 0;
    
    % save a backup before median filtering
    I_filtered_h_bk = I_filtered_h;
    
    % 3x3 median fiter on the filtered image
    I_filtered_h = medfilt2(I_filtered_h, [3 3]); % (READ POINT: I_filtered_h)
    
    %% vertical direction
    % padding borders via extra interpolation before filtering
    IDbin_pad_v = zeros(size(IDbin,1)+fw-1, size(IDbin,2));
    for j = 1:size(IDbin,2)
        IDbin_pad_v(1:(fw+1)/2,j)       = interp1(0:(fw+1)/2-1, IDbin(1:(fw+1)/2,j), -(fw-1)/2:0, 'linear', 'extrap');
        IDbin_pad_v(end-(fw-1)/2:end,j) = interp1(-(fw-1)/2:0, IDbin(end-(fw-1)/2:end,j), 0:(fw+1)/2-1, 'linear', 'extrap');
    end
    IDbin_pad_v((fw+1)/2:end-(fw-1)/2,:) = IDbin;
    
    % 1x3 median filter (horizontal)
    I_medfilt_v = medfilt2(IDbin_pad_v,[1 3]);
    
    % vertical filtering
    I_filtered_pad_v = imfilter(I_medfilt_v, h', 'replicate');
    % crop back to original size
    I_filtered_v = I_filtered_pad_v((fw+1)/2:end-(fw-1)/2,:);
    
    I_filtered_v(I_filtered_v<0) = 0;
    
    % save a backup before median filtering
    I_filtered_v_bk = I_filtered_v;
    
    % 3x3 median filter on the filtered image
    I_filtered_v = medfilt2(I_filtered_v, [3 3]); % (READ POINT: I_filtered_v)
    
    %% combine veritcal and horizontal results (edge, center, and corners are treated differently)
    
    I_filtered_c = (I_filtered_h+I_filtered_v)/2; % average
    % 4 edges
    I_filtered_c((fw+1)/2:end-(fw-1)/2,1:(fw+1)/2)        = I_filtered_v((fw+1)/2:end-(fw-1)/2,1:(fw+1)/2);
    I_filtered_c((fw+1)/2:end-(fw-1)/2,end-(fw-1)/2:end)  = I_filtered_v((fw+1)/2:end-(fw-1)/2,end-(fw-1)/2:end);
    I_filtered_c(1:(fw+1)/2, (fw+1)/2:end-(fw-1)/2)       = I_filtered_h(1:(fw+1)/2, (fw+1)/2:end-(fw-1)/2);
    I_filtered_c(end-(fw-1)/2:end, (fw+1)/2:end-(fw-1)/2) = I_filtered_h(end-(fw-1)/2:end, (fw+1)/2:end-(fw-1)/2);
    % 4 corners
    I_filtered_c(1:(fw+1)/2, 1:(fw+1)/2)             = min(I_filtered_h(1:(fw+1)/2, 1:(fw+1)/2),I_filtered_v(1:(fw+1)/2, 1:(fw+1)/2));
    I_filtered_c(1:(fw+1)/2, end-(fw-1)/2:end)       = min(I_filtered_h(1:(fw+1)/2, end-(fw-1)/2:end), I_filtered_v(1:(fw+1)/2, end-(fw-1)/2:end));
    I_filtered_c(end-(fw-1)/2:end, 1:(fw+1)/2)       = min(I_filtered_h(end-(fw-1)/2:end, 1:(fw+1)/2), I_filtered_v(end-(fw-1)/2:end, 1:(fw+1)/2));
    I_filtered_c(end-(fw-1)/2:end, end-(fw-1)/2:end) = min(I_filtered_h(end-(fw-1)/2:end, end-(fw-1)/2:end), I_filtered_v(end-(fw-1)/2:end, end-(fw-1)/2:end));
    
    %% keep values of severe defects
    % location of severe defects
    indx = I_filtered_h_bk >= threshold | I_filtered_v_bk >= threshold;
    
    % keep the values of those severe defects
    I_filtered_c(indx)   = max(I_filtered_h_bk(indx),I_filtered_v_bk(indx)); % (READ POINT: I_filtered_c)
    
    %% corner mask by radius
    % find centre region and calculate distances
    I_filter_med=median(I_filtered_c(:));
    centreX_hm = rows/2 + 0.5;
    centreY_hm = cols/2 + 0.5;
    halfDiag_hm = sqrt((centreX_hm-1)^2 + (centreY_hm-1)^2);
    for ic=1:rows
        for jc=1:cols
            Fpixel(ic,jc) = sqrt((centreX_hm-ic)^2 + (centreY_hm-jc)^2) / halfDiag_hm;
            if Fpixel(ic,jc)>=field_effective
                I_filtered_c(ic,jc)=I_filter_med;
                % debug plot corner
                %                 plot(jc,ic,'or')
            end
        end
    end
    % (READ POINT: I_filtered_c)
    
    %% Predefined corner size for reporting results of edge, center, and corners
    % define center
    CEN_ROW = cornerY_size+1:rows-cornerY_size;
    CEN_COL = cornerX_size+1:cols-cornerX_size;
    
    % define 4 corners
    UL = [1:cornerY_size;           1:cornerX_size];
    UR = [1:cornerY_size;           cols-cornerX_size+1:cols];
    LL = [rows-cornerY_size+1:rows; 1:cornerX_size];
    LR = [rows-cornerY_size+1:rows; cols-cornerX_size+1:cols];
    
    % define 4 edges
    EDG_L_ROW = cornerY_size+1:rows-cornerY_size;
    EDG_L_COL = 1:cornerX_size;
    EDG_R_ROW = cornerY_size+1:rows-cornerY_size;
    EDG_R_COL = cols-cornerX_size+1:cols;
    EDG_T_ROW = 1:cornerY_size;
    EDG_T_COL = cornerX_size+1:cols-cornerX_size;
    EDG_B_ROW = rows-cornerY_size+1:rows;
    EDG_B_COL = cornerX_size+1:cols-cornerX_size;
    
    %% create center, edge, and corner masks
    
    % center region
    % (make one image with the same size as the whole filtered
    % image, but 4 edges and corners are zeros, this is easier to find the location of centerMax)
    I_filtered_c_cen = zeros(size(I_filtered_c));
    I_filtered_c_cen(CEN_ROW,CEN_COL)   = I_filtered_c(CEN_ROW,CEN_COL);
    
    % edge region
    % (make one image with the same size as the whole filtered
    % image, but center and 4 corners are zeros, this is easier to find the location of edgeMax)
    I_filtered_c_edg = zeros(size(I_filtered_c));
    I_filtered_c_edg(EDG_L_ROW,EDG_L_COL) = I_filtered_c(EDG_L_ROW,EDG_L_COL);
    I_filtered_c_edg(EDG_R_ROW,EDG_R_COL) = I_filtered_c(EDG_R_ROW,EDG_R_COL);
    I_filtered_c_edg(EDG_T_ROW,EDG_T_COL) = I_filtered_c(EDG_T_ROW,EDG_T_COL);
    I_filtered_c_edg(EDG_B_ROW,EDG_B_COL) = I_filtered_c(EDG_B_ROW,EDG_B_COL);
    
    % 4 corners (will have different spec from cen/edge)
    % (make one image with the same size as the whole filtered
    % image, but center and 4 edges are zeros, this is easier to find the location of cornerMax)
    I_filtered_c_cor = zeros(size(I_filtered_c));
    I_filtered_c_cor(UL(1,:), UL(2,:)) = I_filtered_c(UL(1,:), UL(2,:));
    I_filtered_c_cor(UR(1,:), UR(2,:)) = I_filtered_c(UR(1,:), UR(2,:));
    I_filtered_c_cor(LL(1,:), LL(2,:)) = I_filtered_c(LL(1,:), LL(2,:));
    I_filtered_c_cor(LR(1,:), LR(2,:)) = I_filtered_c(LR(1,:), LR(2,:));

    % max intensity directly from difference image
    % these are the primary metrics for blemish detection
    LCBraw_center = max(I_filtered_c_cen(:));
    LCBraw_corner = max(I_filtered_c_cor(:));
    LCBraw_edge   = max(I_filtered_c_edg(:)); 
end

% Sam, 220510
function [col, row] = FindLoctionOfMaxVal(fig, filename, filterWidth, I_filtered_c, col, row, feature_roiSize, center_first)
% original location, Sam, 220512
% [col, row] = find((I_filtered_c_cen') == max(I_filtered_c_cen(:)),1);
if fig>0
    clf(figure(fig));
end

[rows, cols] = size(I_filtered_c);

filterWidth = 5;

if filterWidth == 5 % HF
    % select top N maximum pixel value, Sam, 220505
    [pix, ~] = sort(I_filtered_c(:), 'descend');
    idx = 1; search_num = 7;
    dist_th = 20;
    consider_sft = 2;
    %col_t = []; row_t = []; area_value_t = [];
    col_t = zeros(1, search_num);
    row_t = zeros(1, search_num);
    area_value_t = zeros(1, search_num);
    search_idx = 1;
%     while length(col_t) < search_num
    while search_idx <= search_num % Sam, 220517
        max_value_t = pix(idx);
        [col_, row_] = find((I_filtered_c') == max_value_t, 1);
%         [col_, row_] = find((I_filtered_c_cen_blur') == max_value_t, 1);
        
        if idx==1
            col_t(search_idx) = col_; row_t(search_idx) = row_; % Sam, 220517
%             col_t = [col_t, col_]; row_t = [row_t, row_];
            col_min = max(1,col_-consider_sft); col_max = min(col_+consider_sft,cols);
            row_min = max(1,row_-consider_sft); row_max = min(row_+consider_sft,rows);
            area_roi = I_filtered_c(row_min:row_max, col_min:col_max);
%             area_value_t(search_idx) = mean(area_roi(:)); % Sam, 220517
            area_value_t(search_idx) = sum(area_roi(:)); % Sam, 220517
%             area_roi_we = we .* area_roi;
%             area_value_t(search_idx) = sum(area_roi_we(:)); % Sam, 220519
            
            search_idx = search_idx + 1; % Sam, 220517
        else
            similar = 0;
            for h = 1:length(col_t)
                if area_value_t(h) == 0 % avoid area value with zero, Sam, 20519
                    continue;
                end
                dist = (col_-col_t(h))^2 + (row_-row_t(h))^2;
                if dist <= dist_th
                    similar = 1;
                    break;
                end
            end
            if similar == 0
                col_t(search_idx) = col_; row_t(search_idx) = row_; % Sam, 220517
%               col_t = [col_t, col_]; row_t = [row_t, row_];
                col_min = max(1,col_-consider_sft); col_max = min(col_+consider_sft,cols);
                row_min = max(1,row_-consider_sft); row_max = min(row_+consider_sft,rows);
                area_roi = I_filtered_c(row_min:row_max, col_min:col_max);
%                area_value_t(search_idx) = mean(area_roi(:)); % Sam, 220517
                area_value_t(search_idx) = sum(area_roi(:)); % Sam, 220517
%                 area_roi_we = we .* area_roi;
%                 area_value_t(search_idx) = sum(area_roi_we(:)); % Sam, 220519

                search_idx = search_idx + 1; % Sam, 220517
            end
        end
        
        idx = idx + 1;
    end
    
%     [area_value_sort, area_value_idx] = sort(area_value_t(:), 'descend');
    [~, area_value_idx] = sort(area_value_t(:), 'descend');
    
    % use full roi first if possible, Sam, 220517
    search_idx = 1;
    if center_first == 1
        for idx = 1:length(area_value_idx)
%             col = col_t(area_value_idx(idx));
%             row = row_t(area_value_idx(idx));

            xmin_bin_roi = max(col-round(feature_roiSize/2),1);
            xmax_bin_roi = min(col+round(feature_roiSize/2),cols);
            ymin_bin_roi = max(row-round(feature_roiSize/2),1);
            ymax_bin_roi = min(row+round(feature_roiSize/2),rows);

            w = round(ymax_bin_roi-ymin_bin_roi);
            h = round(xmax_bin_roi-xmin_bin_roi);
            if w >= feature_roiSize && h >= feature_roiSize
                search_idx = idx;
                break;
            end
        end
    end
    
%     fprintf(fid, '%s, %d, %d, %d, %d,', filename, col, row, col_t(area_value_idx(1)), row_t(area_value_idx(1)));

%     figure(8), subplot(1,3,1), % Sam, 220510
    if fig>0
        figure(fig),
        annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');

        subplot(1,2,1), imshow(I_filtered_c, []); title({'orignial possible', 'location'}, 'FontSize', 8); % Sam, 220512
        rectangle('position', [col-2, row-2, 5, 5], 'EdgeColor', 'r');
    %     set(gcf,'position',[800,500,400,400]),
        subplot(1,2,2), imshow(I_filtered_c, []); title({'optimized possible', 'location'}, 'FontSize', 8); % Sam, 220510
        rectangle('position', [col_t(1)-2, row_t(1)-2, 5, 5], 'EdgeColor', 'r');
        rectangle('position', [col_t(2)-2, row_t(2)-2, 5, 5], 'EdgeColor', 'g');
        rectangle('position', [col_t(3)-2, row_t(3)-2, 5, 5], 'EdgeColor', 'c');
        rectangle('position', [col_t(4)-2, row_t(4)-2, 5, 5], 'EdgeColor', 'y');
        rectangle('position', [col_t(5)-2, row_t(5)-2, 5, 5], 'EdgeColor', 'm');
        rectangle('position', [col_t(6)-2, row_t(6)-2, 5, 5], 'EdgeColor', 'b');
        rectangle('position', [col_t(7)-2, row_t(7)-2, 5, 5], 'EdgeColor', [0.4660 0.6740 0.1880]);
        %rectangle('position', [col_t(area_value_idx(1))-4, row_t(area_value_idx(1))-4, 9, 9], 'EdgeColor', [1 0.6823 0.7882], 'LineWidth', 1, 'LineStyle', ':');
        rectangle('position', [col_t(area_value_idx(search_idx))-4, row_t(area_value_idx(search_idx))-4, 9, 9], 'EdgeColor', [1 0.6823 0.7882], 'LineWidth', 1);
    end
    
    % replace to optimize location, Sam, 220505
%     fprintf('possible location:\r');
%     fprintf('  ABC : %d, %d\r', col, row);
%     fprintf('  DPBU: %d, %d\r', col_t(area_value_idx(1)), row_t(area_value_idx(1)));
    % can't find full roi, use first one, Sam, 220517
    col = col_t(area_value_idx(search_idx));
    row = row_t(area_value_idx(search_idx));
else
    I_filtered_c_blur = medfilt2(I_filtered_c, [3, 3], 'symmetric');
%     I_filtered_c_cor_blur = medfilt2(I_filtered_c_cor, [3, 3], 'symmetric');
%     I_filtered_c_edg_blur = medfilt2(I_filtered_c_edg, [3, 3], 'symmetric');

    % search maximum value from an roi instead of a pixel, Sam, 220427
    consider_sft = 2;
    roi_sum_max = -999;
    col_ = 0;
    row_ = 0;
    for j = 1+consider_sft:rows-consider_sft
        roi_y_min = j-consider_sft;
        roi_y_max = j+consider_sft;
        
        for i = 1+consider_sft:cols-consider_sft
            roi_x_min = i-consider_sft;
            roi_x_max = i+consider_sft;

%             roi_i = I_filtered_c_cen(roi_y_min:roi_y_max, roi_x_min:roi_x_max);
            roi_i = I_filtered_c_blur(roi_y_min:roi_y_max, roi_x_min:roi_x_max); % Sam, 220427
            roi_sum = sum(roi_i(:));
            
            if roi_sum > roi_sum_max
                roi_sum_max = roi_sum;
                col_ = i;
                row_ = j;
            end
        end
    end   
%     fprintf(fid, '%s,%d,%d,%d,%d,', filename, col, row, col_, row_);
    
    if fig>0
    %     figure(8), subplot(1,3,1), % Sam, 220510
        figure(fig),
        annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');

        subplot(1,2,1), imshow(I_filtered_c, []), title({'orignial possible', 'location'}, 'FontSize', 8), % Sam, 220512
        rectangle('position', [col-2, row-2, 5, 5], 'EdgeColor', 'r');
    %     set(gcf,'position',[800,500,400,400]),
        subplot(1,2,2), imshow(I_filtered_c, []), title({'optimized possible location'}, 'FontSize', 8), % Sam, 220510
    %     rectangle('position', [col-2, row-2, 5, 5], 'EdgeColor', 'r');
    %    rectangle('position', [col_-4, row_-4, 9, 9], 'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', ':');
        rectangle('position', [col_-4, row_-4, 9, 9], 'EdgeColor', [1 0.6823 0.7882], 'LineWidth', 1);
    end
    
    % replace col and row to col_ and row_, Sam, 220427
%     fprintf('possible location:\r');
%     fprintf('  ABC : %d, %d\r', col, row);
%     fprintf('  DPBU: %d, %d\r', col_, row_);
    col = col_;
    row = row_;
end

% saveas(figure(2), strcat(index, filename, '0_PossLoc.png')); % Sam, 220506
    
end

%% DP locator
function [LCBraw_x, LCBraw_y, col, row, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
    DP_Locator(fig, filename, region, I_filtered_c, roiSize, feature_roiSize, fw, rows, cols)

    %% find the location with max center score in the scalded image
    [col row] = find((I_filtered_c') == max(I_filtered_c(:)), 1);
    
    roiX_size = roiSize(2);
    roiY_size = roiSize(1);
    
    %% locate position of maximum value, Sam, 220510
    imgname = sprintf('3_%s', filename);
    switch region
        case 0 % center
            imgname = sprintf('3_0_cen_%s', filename);
        case 1 % edge
            imgname = sprintf('3_1_edg_%s', filename);
        case 2 % corner
            imgname = sprintf('3_2_cor_%s', filename);
    end   
%     [col_opt, row_opt] = FindLoctionOfMaxVal(index, imgname, fw, I_filtered_c, col, row, feature_roiSize, 0);
    [col_opt, row_opt] = FindLoctionOfMaxVal(fig, imgname, fw, I_filtered_c, col, row, feature_roiSize, 0);
    col = col_opt;
    row = row_opt;

    % trace back to the location in the full size image of current color plane
    defectLoc_edge_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_x = round(defectLoc_edge_q(1));
    LCBraw_y = round(defectLoc_edge_q(2));

    % get roi for feature detection
    % define roi coordinates
    xmin_bin_roi = max(col-round(feature_roiSize/2),1);
    xmax_bin_roi = min(col+round(feature_roiSize/2),cols);
    ymin_bin_roi = max(row-round(feature_roiSize/2),1);
    ymax_bin_roi = min(row+round(feature_roiSize/2),rows);
    
    % extract roi
    roi = I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
end

%% ABC locator
function [LCBraw_x, LCBraw_y, col, row, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
    ABC_Locator(I_filtered_c, roiSize, feature_roiSize, fw, rows, cols)
    %------------------------------------------
    % add defect location info for center
    %------------------------------------------
    roiX_size = roiSize(2);
    roiY_size = roiSize(1);
    
    % find the location with max center score in the scalded image
    [col row] = find((I_filtered_c') == max(I_filtered_c(:)),1);
    
    % trace back to the location in the full size image of current color plane
    defectLoc_edge_cen_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_x = round(defectLoc_edge_cen_q(1)) ;
    LCBraw_y = round(defectLoc_edge_cen_q(2)) ;
    
    % get roi for feature detection
    % define roi coordinates 
    xmin_bin_roi = max(col-round(feature_roiSize/2),1);
    xmax_bin_roi = min(col+round(feature_roiSize/2),cols);
    ymin_bin_roi = max(row-round(feature_roiSize/2),1);
    ymax_bin_roi = min(row+round(feature_roiSize/2),rows);
    
    % extract roi
    roi = I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
end

%% ABC Detector
% region: 0: center, 1: edge, 2: cornor
function [LCBraw_x, LCBraw_y, row, col, IDbin_roi, output_particle] = ...
    ABC_Detector(index, region, filename, IDraw, IDbin, I_filtered_c, roiSize, fw, rows, cols, feature_roiSize)
    
    [LCBraw_x, LCBraw_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        ABC_Locator(I_filtered_c, roiSize, feature_roiSize, fw, rows, cols);

    % calculate particle statistics within ROI
    [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
        xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
%     output_particle_cell = num2cell(output_particle.stats);
%     [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
end

%% DP Detector
% region: 0: center, 1: edge, 2: cornor
function [LCBraw_x, LCBraw_y, row, col, IDbin_roi, output_particle] = ...
    DP_Detector(index, region, filename, IDraw, IDbin, I_filtered_c, roiSize, fw, rows, cols, feature_roiSize)
    
    [LCBraw_x, LCBraw_y, row, col, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roi] = ...
        DP_Locator(I_filtered_c, roiSize, feature_roiSize, fw, rows, cols);

    % calculate particle statistics within ROI
    [IDbin_roi, output_particle] = getParticleStatsFromLCBmap( filename, roi, IDraw, IDbin, ...
        xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
%     output_particle_cell = num2cell(output_particle.stats);
%     [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% getParticleStatsFromLCBmap
%   calculate particle statics from LCB heat map to raw images
%
% INPUT:
%   roi - ROI from LCB heat map (binned size) for detecting feature
%   IDraw - Raw image (after preprocess)
%   xmin_bin_roi - xmin of ROI in binned image coordinates
%   xmax_bin_roi - xmax of ROI in binned image coordinates
%   ymin_bin_roi - ymin of ROI in binned image coordinates
%   ymax_bin_roi - ymax of ROI in binned image coordinates
%   roiX_size - roiX_size of Binning, LCBraw test input
%   roiY_size - roiY_size of Binning, LCBraw test input
%
% OUTPUT:
%   particle_intensity
%   image_level
%   particle_attenuation_raw
%   particle_h_raw
%   particle_w_raw
%   particle_area_raw
%
%   Author: Aaron Cong
%           March 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [IDbin, output] = getParticleStatsFromLCBmap(filename, roi, IDraw, IDbin, ...
    xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size)
debug = 0;
BorW=1; % always find brighter feature in LCBraw heat map
% clear & initialize output parameters
particle_intensity =0;
image_level=0;
particle_attenuation_raw=0;
particle_h_raw=0;
particle_w_raw=0;
particle_area_raw=0;
% get roi size
[hROI wROI] = size(roi);
% fileID = fopen('D:/M_roi.csv', 'w+b');
% for index = 1:size(roi,1)
%     fprintf(fileID, '%f,', roi(index, :));
%     fprintf(fileID, '\r\n');
% end
% fclose(fileID);
% Median filter to clean up noise or stray pixels
roi = medfilt2(roi, [3, 3]);
% Threshold ROI, invert if detecting white

%% convert to 8 bit
bitDepth=8;
roi=roi./max(roi(:))*2^bitDepth;
roiBin = thresholdROI(roi, bitDepth, BorW);
% fileID = fopen('D:/M_roiBin.csv', 'w+b')
% for index = 1:size(roiBin,1)
%     fprintf(fileID, '%d,', roiBin(index, :));
%     fprintf(fileID, '\r\n');
% end
% fclose(fileID);
if ~(BorW)
    roiBin = ~roiBin;
end
% Sam, 221102, material for request of Alex
% roiBin_full = zeros(187,157);
IDbin(ymin_bin_roi:ymax_bin_roi, xmin_bin_roi:xmax_bin_roi) = ~roiBin;
%imwrite(roiBin_full, sprintf('%s_roiBin_xmin=%d_ymin=%d_xmax=%d_ymax=%d.bmp', filename, xmin_bin_roi, ymin_bin_roi, xmax_bin_roi, ymax_bin_roi));
% roiBin = cleanROI(roiBin);
% find coordinates of particle
[y , x]=find(roiBin==0);
if isempty(y)
    ymax_bin_particle=0;
    ymin_bin_particle=0;
    xmax_bin_particle=0;
    xmin_bin_particle=0;
    particle_h=0;
    particle_w=0;
else
    ymax_bin_particle=max(y);
    ymin_bin_particle=min(y);
    xmax_bin_particle=max(x);
    xmin_bin_particle=min(x);
    particle_h=ymax_bin_particle-ymin_bin_particle;
    particle_w=xmax_bin_particle-xmin_bin_particle;
    % parcile size is limited to half size of the ROI detection area
    if particle_h>hROI/2 && particle_w>wROI/2
        ymax_bin_particle=0;
        ymin_bin_particle=0;
        xmax_bin_particle=0;
        xmin_bin_particle=0;
        particle_h=0;
        particle_w=0;
    end
end
xmin_raw_particle = 0;
ymin_raw_particle = 0;
% get particle size in raw dimension
if (xmin_bin_particle* xmax_bin_particle* ymin_bin_particle* ymax_bin_particle* particle_h* particle_w)~=0 % check if any output is 0, which is default flag for no particle found
    % convert x address from binned image coordinates to raw image coordinates
    xmin_raw_particle=(xmin_bin_roi+xmin_bin_particle-1)*roiX_size-roiX_size/2;
    xmax_raw_particle=(xmax_bin_roi+xmax_bin_particle-1)*roiX_size-roiX_size/2;
    particle_w_raw=particle_w*roiX_size;
    
    % convert y address from binned image coordinates to raw image coordinates
    ymin_raw_particle=(ymin_bin_roi+ymin_bin_particle-1)*roiY_size-roiY_size/2;
    ymax_raw_particle=(ymax_bin_roi+ymax_bin_particle-1)*roiY_size-roiY_size/2;
    particle_h_raw=particle_h*roiY_size;
    
    % calculate area of particle defined in unit of pixels
    particle_area_raw = particle_w_raw * particle_h_raw;
    
    % calculate particle intensity by taking the median of half height, half width aread in
    % the center of ROI
    
    IDparticle=IDraw(round(ymin_raw_particle):round(ymin_raw_particle+particle_h_raw),round(xmin_raw_particle):round(xmin_raw_particle+particle_w_raw));
    
    % calculate percentiles
    [pix, idx] = sort(IDparticle(:));
    numPix = length(pix);
    idx05 = round(0.05 * (numPix+1));
    particle_intensity=pix(max(idx05,1));
    
    % calcualte median of entire image (after LSC)
    image_level=median(IDraw(:));
    
    % calculate contrast of particle as the delta/image median (unit is
    % percentage) Negative values showing intensity of particle area is lower than image
    % level
    particle_attenuation_raw=particle_intensity/image_level-1;
    if (debug)
        fprintf('intensity ratio = %s \n',num2str(particle_attenuation_raw));
        figureDebug=figure('position', [0, 0, 1200, 350]) ;
        figure(figureDebug);
        subplot(1,3,1)
        % plot detected particle in binned image
        imagesc(roi);colorbar
        rectangle('position',[xmin_bin_particle,ymin_bin_particle,particle_w,particle_h],'EdgeColor','r');
        text(xmin_bin_particle,ymin_bin_particle-1,['w=',num2str(particle_w),' , ','h=',num2str(particle_h)],'Color','r')
        title('Deteced particle in filtered binned image')
        xlabel('ROI column index 1 -> 50 in binned image coordinates')
        ylabel('ROI row index 1 -> 50 in binned image coordinates')
        subplot(1,3,2)
        % plot detected particle in raw image & show contrast & size of
        % (height, width, area) in unit of pixels
        imagesc(IDraw);colorbar
        rectangle('position',[xmin_raw_particle,ymin_raw_particle,particle_w_raw,particle_h_raw],'EdgeColor','r');
        text(xmin_raw_particle,ymin_raw_particle-10,[num2str(particle_h_raw),',',num2str(particle_w_raw),',',num2str(particle_area_raw)],'Color','r')
        text(xmin_raw_particle+particle_w_raw/2,ymin_raw_particle+particle_h_raw/2-1,num2str(particle_attenuation_raw),'Color','r')
        title('Raw image heatmap with particle size & attenuation indication')
        xlabel('Column index 1 -> 1104 in raw image coordinates')
        ylabel('Row index 1 -> 1312 in raw image coordinates')
        subplot(1,3,3)
        % compare histogram of particle area & of entire raw image
        h1=histogram(IDraw(:));hold on
        h2=histogram(IDparticle(:));
        h1.BinLimits=[650 850];
        h1.BinWidth=5;
        h2.BinWidth=5;
        h1.Normalization= 'probability';
        h2.Normalization= 'probability';
        stem(particle_intensity,max([h1.Values,h2.Values]))
        text(particle_intensity-50,max([h1.Values,h2.Values]),['5% Hist =',num2str(particle_intensity)])
        grid on;
        legend('Image Histogram','Particle Histogram')
        title('particle histogram v.s. image histogram')
        xlabel('image pixel level [LSB]')
        ylabel('Normalized count by distribution')
    end
else
    if debug % add debug condition, Sam, 220411
        fprintf('no particle detected! \n')
        figureDebug =[];
    end
end
output.stats=[particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, xmin_raw_particle,ymin_raw_particle] ;
%output.debugPlot=figureDebug;
end
