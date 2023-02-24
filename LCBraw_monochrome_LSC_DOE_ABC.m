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

function [output_final] = LCBraw_monochrome_LSC_DOE_ABC(filename, IDraw, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold,field_effective, option_lsc, ls_mat)
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
%filterWidth = 9;
%threshold = 16;
%field_effective = 0.93;
feature_roiSize= 50; % square size
%option_lsc = 3;
test_ver=3.01;
debug = 0; % from 1 to 0, Sam, 220411
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
IDraw = preprocess(IDraw, bayerFormat, 'raw', bitDepth, pedestal, 83.5, false, true, option_lsc,ls_mat); % (READ POINT: IDraw)
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
[rows cols c] = size(IDbin_all);
output = []; %output array
for i=1:c
    %select 1 channel
    IDbin = IDbin_all(:,:,i);
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
    
    %% output
    
    % max intensity directly from difference image
    % these are the primary metrics for blemish detection
    LCBraw_center = max(I_filtered_c_cen(:));
    LCBraw_corner = max(I_filtered_c_cor(:));
    LCBraw_edge   = max(I_filtered_c_edg(:));
    
    %% Add location reporting
    % Sam, 221102, material of request from Alex
    imwrite(I_filtered_c, sprintf('%s_I_filtered_c.bmp', filename));
    
    %------------------------------------------
    % add defect location info for center
    %------------------------------------------
    % find the location with max center score in the scalded image
    [col row] = find((I_filtered_c_cen') == max(I_filtered_c_cen(:)),1);
    
    % trace back to the location in the full size image of current color plane
    defectLoc_edge_cen_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_center_x = round(defectLoc_edge_cen_q(1)) ;
    LCBraw_center_y = round(defectLoc_edge_cen_q(2)) ;
    
    % get roi for feature detection
    % define roi coordinates 
    xmin_bin_roi= max(col-round(feature_roiSize/2),1);
    xmax_bin_roi= min(col+round(feature_roiSize/2),cols);
    ymin_bin_roi=max(row-round(feature_roiSize/2),1);
    ymax_bin_roi=min(row+round(feature_roiSize/2),rows);
    
    % extract roi
    roi= I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
    
%     fileID = fopen('D:/M_I_filtered_c.csv', 'w+b');
%     for index = 1:size(I_filtered_c,1)
%         fprintf(fileID, '%f,', I_filtered_c(index, :));
%         fprintf(fileID, '\r\n');
%     end
%     fclose(fileID);
    % calculate particle statistics within ROI
    output_particle = getParticleStatsFromLCBmap( filename, roi, IDraw,xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
    output_particle_cell=num2cell(output_particle.stats);
    [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
    
    clear row col roi;
    
    %----------------------------------------
    % add defect location info for corners
    %----------------------------------------
    % determine the location with max corner score in the scaled image
    [col row] = find((I_filtered_c_cor') == max(I_filtered_c_cor(:)),1);
    
    % trace back to the location in the full size image of current color plane
    defectLoc_cor_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_corner_x = round(defectLoc_cor_q(1)) ;
    LCBraw_corner_y = round(defectLoc_cor_q(2)) ;
    clear row col;
    
    %----------------------------------------
    % add defect location info for edges
    %----------------------------------------
    % determine the location with max corner score in the scaled image
    [col row] = find((I_filtered_c_edg') == max(I_filtered_c_edg(:)),1);
    
    % trace back to the location in the full size image of current color plane
    defectLoc_edg_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_edge_x = round(defectLoc_edg_q(1)) ;
    LCBraw_edge_y = round(defectLoc_edg_q(2)) ;
    
    % if no center particle is found (conditon particle_intensity ==0) from previous center heat map detection, then start
    % edge particle detection
    if particle_intensity==0
        % get roi for feature detection
        % define roi coordinates
        xmin_bin_roi= max(col-round(feature_roiSize/2),1);
        xmax_bin_roi= min(col+round(feature_roiSize/2),cols);
        ymin_bin_roi=max(row-round(feature_roiSize/2),1);
        ymax_bin_roi=min(row+round(feature_roiSize/2),rows);
        % extract roi
        roi= I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
        % calculate particle statistics within ROI
        output_particle = getParticleStatsFromLCBmap( filename, roi, IDraw,xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
        output_particle_cell=num2cell(output_particle.stats);
        [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
    end
    if debug
        figure;imshow(IDraw, []);
        hold on;
        if particle_area_x > 0 && particle_area_y > 0 && particle_h_raw > 0 && particle_w_raw > 0
        rectangle('position', [particle_area_x, particle_area_y, particle_h_raw, particle_w_raw], 'EdgeColor', 'r');
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 - 15), ['w: ' num2str(particle_w_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        text(round(particle_area_x), round(particle_area_y + particle_h_raw / 2 + 15), ['h: ' num2str(particle_h_raw, '%d')], 'Color', 'g', 'FontSize', 10);
        end
        scatter([LCBraw_center_x LCBraw_edge_x LCBraw_corner_x],[LCBraw_center_y LCBraw_edge_y LCBraw_corner_y],'yo');
        hold off;
    end
%     clear row col roi;
    %----------------------------------------
    % output final results
    %----------------------------------------
    output_temp = [...
        LCBraw_center LCBraw_center_x LCBraw_center_y,...         % original LCB scores
        LCBraw_edge   LCBraw_edge_x   LCBraw_edge_y,...
        LCBraw_corner LCBraw_corner_x LCBraw_corner_y,...
        particle_intensity image_level particle_attenuation_raw,...         %particle statistics
        particle_h_raw particle_w_raw particle_area_raw];
    
    % output format = [R;Gr;Gb;B], not applicable to RAW output processing
    output = [output;output_temp]; % (READ POINT: ouput)
end

output_final.testLog=[output test_ver];
output_final.outputTrends=getLCBrawLSCOutputTrends(filterWidth);
% writeResultData(output_final.testLog, pathname, filename);
debug = 0;
if debug
    output_final.IDrawLSC=IDraw;
    % plot regular LCB heatmap
    output_final.debugPlot = figure('position', [0, 0, 1800, 900]) ; hold on; colour=[];
    subplot(1,2,1)
    set(gca,'YDir','reverse');
    %scale=std2(I_filtered_c)*20;
    scale=5;
    imagesc(I_filtered_c,[0 scale]), colorbar; hold on;
    title({'Heat map of filtered binned image',colour});
    xlabel('Column index 1 -> 157 in binned image coordinates')
    ylabel('Row index 1 -> 187 in binned image coordinates')
    text(round(defectLoc_edge_cen_q(1)/roiX_size+1/2), round(defectLoc_edge_cen_q(2)/roiY_size+1/2),[num2str(max(I_filtered_c_cen(:)),'%6.4f')],'Color','r', 'FontSize', 15);
    text(round(defectLoc_edg_q(1)/roiX_size+1/2), round(defectLoc_edg_q(2)/roiY_size+1/2),[num2str(max(I_filtered_c_edg(:)),'%6.4f')],'Color','g', 'FontSize', 15);
    text(round(defectLoc_cor_q(1)/roiX_size+1/2), round(defectLoc_cor_q(2)/roiY_size+1/2),[num2str(max(I_filtered_c_cor(:)),'%6.4f')],'Color','c', 'FontSize', 15);
    subplot(1,2,2)
%     imshow(IDraw, []);
     surf(1:size(IDraw,2), 1:size(IDraw,1), IDraw);colorbar;shading interp
%     imagesc(IDraw,[700 850]), colorbar;
    xlabel('Column index 1 -> 1104 in raw image coordinates')
    ylabel('Row index 1 -> 1312 in raw image coordinates')
    title('Heat map of raw image after shading correction')
     hgexport(output_final.debugPlot, [pathname '/' filename '.png'], hgexport('factorystyle'), 'Format', 'png');
%     saveas(output_final.debugPlot, [pathname '/' filename '.png']);
%     print(output_final.debugPlot, '-dpng','-r72', [pathname '/' filename '.png']);
    % plot particle detection information
%     if isempty(output_particle.debugPlot)
%         output_final.particlePlot= [];
%     else
%         output_final.particlePlot= output_particle.debugPlot;
%     end
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
    for j = 1:length(imbin(:,1,i));
        imbin(j,1,i) = mean2( im(ye(j):ye(j+1)-1, xe(1):xe(2)-1,i) );
        imbin(j,end,i) = mean2( im(ye(j):ye(j+1)-1, xe(end-1):xe(end)-1,i) );
    end
    
    % First and last row
    for j = 1:length(imbin(1,:,i));
        imbin(1,j,i) = mean2( im(ye(1):ye(2)-1, xe(j):xe(j+1)-1,i) );
        imbin(end,j,i) = mean2( im(ye(end-1):ye(end)-1, xe(j):xe(j+1)-1,i) );
    end
end
end

function [output_trends] = getLCBrawLSCOutputTrends(filterWidth)
output_trends= {...
    ['LCBraw_FW',num2str(filterWidth),'_cen_max_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_cen_max_x_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_cen_max_y_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_edge_max_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_edge_max_x_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_edge_max_y_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_corner_max_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_corner_max_x_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_corner_max_y_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_particle_intensity_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_image_level_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_particle_attenuation_raw_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_particle_h_raw_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_particle_w_raw_[LS]'],...
    ['LCBraw_FW',num2str(filterWidth),'_particle_area_raw_[LS]'],...
    'test_ver'...
    };
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
function output = getParticleStatsFromLCBmap(filename,roi,IDraw,xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi,roiX_size, roiY_size)
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
roiBin_full = zeros(187,157);
roiBin_full(ymin_bin_roi:ymax_bin_roi, xmin_bin_roi:xmax_bin_roi) = roiBin;
imwrite(roiBin_full, sprintf('%s_roiBin_xmin=%d_ymin=%d_xmax=%d_ymax=%d.bmp', filename, xmin_bin_roi, ymin_bin_roi, xmax_bin_roi, ymax_bin_roi));
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

function [bw] = cleanROI(roiBIN)
roiBIN = ~roiBIN;
[h w] = size(roiBIN);
Connected_ROI_Bin = int16(zeros(h,w));
% Find all non-zero pixels
[y_i x_i] = find(roiBIN >= 1);
[linearIndexes] = find(roiBIN >= 1);
currentComponent = 1;
componentConnections = [];
% First pass component labelling
for counter = 1:length(linearIndexes)
    current_x = x_i(counter);
    current_y = y_i(counter);
    % check only rows / columns that are not on the border
    if ((current_y > 1) && (current_y < h) && (current_x > 1) && (current_x < w))
        % Check surrounding components for any labels
        if (Connected_ROI_Bin(current_y+1,current_x-1)==0) && ...
                (Connected_ROI_Bin(current_y,current_x-1)==0) && ...
                (Connected_ROI_Bin(current_y-1,current_x-1)==0) && ...
                (Connected_ROI_Bin(current_y-1,current_x)==0)
            % Component is not labelled, so label it
            Connected_ROI_Bin(current_y, current_x) = currentComponent;
            % Make a component connected to itself (base case for tree search)
            componentConnections(currentComponent,currentComponent)=currentComponent;
            % Increment component counter
            currentComponent = currentComponent+1;
        else
            % Inherit the highest component value in the neighboring region
            componentList = [Connected_ROI_Bin(current_y+1,current_x-1), ...
                Connected_ROI_Bin(current_y,current_x-1), ...
                Connected_ROI_Bin(current_y-1,current_x-1), ...
                Connected_ROI_Bin(current_y-1,current_x)];
            % Inherit the lowest component identifier
            indexes = find(componentList);
            Connected_ROI_Bin(current_y, current_x) = min(componentList(indexes));
            % Record the equivalences in an upper triangular matrix
            for i = 1:length(indexes)
                componentConnections(Connected_ROI_Bin(current_y, current_x),componentList(indexes(i))) = componentList(indexes(i));
            end
        end
    else
        % Borders are labelled as 1's
        %Connected_ROI_Bin(current_y, current_x) = 1;
        componentList = [];
        if(((current_y+1)<h)&&((current_x-1)>0))
            if (Connected_ROI_Bin(current_y+1,current_x-1)>0)
                componentList = [componentList, Connected_ROI_Bin(current_y+1,current_x-1)];
            end
        end
        if (current_x-1)>0
            if (Connected_ROI_Bin(current_y,current_x-1)>0)
                componentList = [componentList, Connected_ROI_Bin(current_y,current_x-1)];
            end
        end
        if(((current_y-1)>0)&&((current_x-1)>0))
            if (Connected_ROI_Bin(current_y-1,current_x-1)>0)
                componentList = [componentList, Connected_ROI_Bin(current_y-1,current_x-1)];
            end
        end
        if (current_y-1)>0
            if (Connected_ROI_Bin(current_y-1,current_x)>0)
                componentList = [componentList, Connected_ROI_Bin(current_y-1,current_x)];
            end
        end
        if isempty(componentList)
            % Component is not labelled, so label it
            Connected_ROI_Bin(current_y, current_x) = currentComponent;
            % Make a component connected to itself (base case for tree search)
            componentConnections(currentComponent,currentComponent)=currentComponent;
            % Increment component counter
            currentComponent = currentComponent+1;
        else
            indexes = find(componentList);
            Connected_ROI_Bin(current_y, current_x) = min(componentList(indexes));
            % Record the equivalences in an upper triangular matrix
            for i = 1:length(indexes)
                componentConnections(Connected_ROI_Bin(current_y, current_x),componentList(indexes(i))) = componentList(indexes(i));
            end
        end
    end
end
% Collapse equivalences into unique trees only
for i = (currentComponent-1):-1:1
    % Loop through each labelled component backwards and set the value
    % to the minimum detected connected component
    lowestJointComponent = min(find(componentConnections(:,i)));
    Connected_ROI_Bin(Connected_ROI_Bin == i) = lowestJointComponent;
end
% Anything connected to the border is masked as zero
%Connected_ROI_Bin(Connected_ROI_Bin == 1) = 0;

% Anything component than 5 pixels gets set to zero
remainderComponents = unique(Connected_ROI_Bin(find(Connected_ROI_Bin)));
components = length(remainderComponents);
if (components > 1)
    for i = 1:components
        indexesOfComponent = find(Connected_ROI_Bin==remainderComponents(i));
        ComponentLength(i) = length(indexesOfComponent);
    end
    [~,index2keep] = max(ComponentLength);
    for i = 1:components
        if i == index2keep
            continue
        end
        indexesOfComponent = find(Connected_ROI_Bin==remainderComponents(i));
        Connected_ROI_Bin(indexesOfComponent) = 0;
    end
end

% Anything left is set to one
Connected_ROI_Bin(Connected_ROI_Bin > 1) = 1;
bw = ~Connected_ROI_Bin;
end
