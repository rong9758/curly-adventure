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

function [output_final] = LCBraw_monochrome_LSC_DOE_DPBU(fid, index, filename, IDraw, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold,field_effective, option_lsc, ls_mat)

%roiSize = [7 7];
%HF:5 LF:9
%filterWidth = 9;you
%threshold = 16;
%field_effective = 0.93;
feature_roiSize = 50; % square size
%option_lsc = 3;
test_ver=3.01;
debug = 0; % from 1 to 0, Sam, 220411

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

    IDbin_o = IDbin; % preserve original IDbin, Sam, 220512

    % enhance contrast, Sam, 220428
    IDbin_avg = mean(IDbin(:));
    IDbin_std = std(IDbin(:));
    IDbin_min = IDbin_avg - 6*IDbin_std;
    IDbin_max = IDbin_avg + 6*IDbin_std;
    IDbin_hist = (IDbin-IDbin_min) ./ (IDbin_max-IDbin_min);

    % get mask of corner, Sam, 220511
    IDbin_hist_rm4cor = Rm4Cor(IDbin_hist, cols, rows);
    
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
       
    %% Add location reporting
    
    %------------------------------------------
    % add defect location info for center
    %------------------------------------------
    % find the location with max center score in the scalded image
    [col row] = find((I_filtered_c_cen') == max(I_filtered_c_cen(:)),1);

    % Sam, 220510
    [col_, row_] = FindLoctionOfMaxVal(index, strcat(filename, '_1_cen'), fw, I_filtered_c_cen, col, row, feature_roiSize, 0);
    fprintf(fid, '%s,%d,%d,%d,%d,', filename, col, row, col_, row_);
    fprintf('possible location of center:\r');
    fprintf('  ABC : %d, %d\r', col, row);
    fprintf('  DPBU: %d, %d\r', col_, row_);
    col = col_;
    row = row_;    
    LCBraw_center = I_filtered_c_cen(row, col); % update pixel value of center, Sam, 220523
    
    
    % trace back to the location in the full size image of current color plane
    defectLoc_edge_cen_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_center_x = round(defectLoc_edge_cen_q(1)) ;
    LCBraw_center_y = round(defectLoc_edge_cen_q(2)) ;
    
    % get roi for feature detection
    % define roi coordinates
    xmin_bin_roi = max(col-round(feature_roiSize/2),1);
    xmax_bin_roi = min(col+round(feature_roiSize/2),cols);
    ymin_bin_roi = max(row-round(feature_roiSize/2),1);
    ymax_bin_roi = min(row+round(feature_roiSize/2),rows);
    
    % extract roi
    roi = I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
    
    % calculate particle statistics within ROI
    output_particle = getParticleStatsFromLCBmap( fid, index, strcat(filename, '_1_cen'), fw, IDbin_o, IDbin_hist, IDbin_hist_rm4cor, roi, IDraw, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size, feature_roiSize);
    output_particle_cell=num2cell(output_particle.stats);
    [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
    
    clear row col roi;
   
    %----------------------------------------
    % add defect location info for corners
    %------------------?----------------------
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
    [col, row] = find((I_filtered_c_edg') == max(I_filtered_c_edg(:)),1);
%     [col row] = find((I_filtered_c_edg_blur') == max(I_filtered_c_edg_blur(:)),1); % Sam, 220425
    
    % Sam, 220510
    [col_, row_] = FindLoctionOfMaxVal(index, strcat(filename, '_2_edg'), fw, I_filtered_c_edg, col, row, feature_roiSize, 0);
%     fprintf(fid, '%s,%d,%d,%d,%d,', filename, col, row, col_, row_);
    fprintf('possible location of edge:\r');
    fprintf('  ABC : %d, %d\r', col, row);
    fprintf('  DPBU: %d, %d\r', col_, row_);
    col = col_;
    row = row_;
    LCBraw_edge = I_filtered_c_edg(row, col); % update pixel value of center, Sam, 220523

    % trace back to the location in the full size image of current color plane
    defectLoc_edg_q = [col*roiX_size-roiX_size/2, row*roiY_size-roiY_size/2];
    % trace back to the location in the full size image of bayer raw (not
    % consider the bayer pattern here, so there will be 1 pixel difference for
    % certain color plane)
    LCBraw_edge_x = round(defectLoc_edg_q(1)) ;
    LCBraw_edge_y = round(defectLoc_edg_q(2)) ;
    
%     disp( strcat('center, particle_intensity:', num2str(particle_intensity))); % Sam, 220424
%     disp(particle_intensity); % Sam, 220424
    % if no center particle is found (conditon particle_intensity ==0) from previous center heat map detection, then start
    % edge particle detection
    if particle_intensity==0
        % get roi for feature detection
        % define roi coordinates
        xmin_bin_roi= max(col-round(feature_roiSize/2),1);
        xmax_bin_roi= min(col+round(feature_roiSize/2),cols);
        ymin_bin_roi=max(row-round(feature_roiSize/2),1);
        ymax_bin_roi=min(row+round(feature_roiSize/2),rows);
        
        % Sam, 220425
%         disp( strcat('edge, col/row: ', num2str(col), ',', num2str(row)) );
%         disp( strcat('edge, roiW/H: ', num2str(xmax_bin_roi-xmin_bin_roi+1), ',', num2str(xmax_bin_roi-xmin_bin_roi+1)) );
%         
        % extract roi
        roi= I_filtered_c(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
        
%         % Sam, 220516
%         [ret, xmin_bin_particle, xmax_bin_particle, ymin_bin_particle, ymax_bin_particle, particle_w, particle_h] = ParticleDetector(IDraw, IDbin, IDbin_hist_rm4cor, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);

        % calculate particle statistics within ROI 
        %output_particle = getParticleStatsFromLCBmap( strcat(filename, '_2'), roi,IDraw,xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size);
        output_particle = getParticleStatsFromLCBmap( '', index, strcat(filename, '_2_edg'), fw, IDbin_o, IDbin_hist, IDbin_hist_rm4cor, roi, IDraw, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size, feature_roiSize);
        output_particle_cell=num2cell(output_particle.stats);
        [particle_intensity, image_level, particle_attenuation_raw, particle_h_raw, particle_w_raw, particle_area_raw, particle_area_x, particle_area_y] = deal(output_particle_cell{:});
        
%         disp( strcat('edge, particle_intensity:', num2str(particle_intensity))); % Sam, 220424
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
adapthisteq
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
function output = getParticleStatsFromLCBmap(fid, index, filename, fw, IDbin_o, IDbin_hist, IDbin_hist_rm4cor, roi, IDraw, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi,roiX_size, roiY_size, feature_roiSize, w_h_ratio_th, std_dist_to_mc_th, r2_edge_filtered_th, roi_top_inc_num_th)
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
roi_o = roi; % Sam, 220610
roi = medfilt2(roi, [3, 3]);
% roi_med = roi; % preserve roi, Sam, 220425
% Threshold ROI, invert if detecting white

%% convert to 8 bit
bitDepth=8;
% roi=roi./max(roi(:))*2^bitDepth;
roi = roi./max(roi(:)) * (2^bitDepth - 1); % Sam, 220425
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
        % Sam, 220516    
        [ret, xmin_bin_particle, xmax_bin_particle, ymin_bin_particle, ymax_bin_particle, particle_w, particle_h, ...
         IDbin_o_roi_particle_avg, IDbin_o_roi_border_avg, IDbin_hist_roi_particle_avg, IDbin_hist_roi_border_avg] = ...
            ParticleDetector(index, filename, IDraw, roi_o, IDbin_o, IDbin_hist, IDbin_hist_rm4cor, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size, feature_roiSize);
        
        fprintf(fid, '%.3f,%.3f,%.3f,%.3f,', ...
                IDbin_o_roi_particle_avg, IDbin_o_roi_border_avg, ...
                IDbin_hist_roi_particle_avg, IDbin_hist_roi_border_avg ); 

        if ret == 0 % it's nosie, Sam, 220505
            ymax_bin_particle=0;
            ymin_bin_particle=0;
            xmax_bin_particle=0;
            xmin_bin_particle=0;
            particle_h=0;
            particle_w=0;
        end
    else
        % Sam, 220623
        fprintf(fid, '%.3f,%.3f,%.3f,%.3f,', '', '', '', '' );
    end
end
xmin_raw_particle = 0;
ymin_raw_particle = 0;
% get particle size in raw dimension
if (xmin_bin_particle* xmax_bin_particle* ymin_bin_particle* ymax_bin_particle* particle_h* particle_w)~=0 % check if any output is 0, which is default flag for no particle found
    % Sam, 220510
    clf(figure(8));
    figure(8),
    annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');
    subplot(1,2,1), imshow(IDbin_o,[]), title({'roi of original IDbin'}, 'FontSize', 8);
    rectangle('position', [xmin_bin_roi, ymin_bin_roi, xmax_bin_roi-xmin_bin_roi+1, ymax_bin_roi-ymin_bin_roi+1], 'EdgeColor', 'b');
    rectangle('position', [xmin_bin_roi+xmin_bin_particle-1, ymin_bin_roi+ymin_bin_particle-1, particle_w, particle_h], 'EdgeColor', 'r');

    subplot(1,2,2), imshow(IDbin_hist), title({'roi of enhanced IDbin'}, 'FontSize', 8);
    rectangle('position', [xmin_bin_roi, ymin_bin_roi, xmax_bin_roi-xmin_bin_roi+1, ymax_bin_roi-ymin_bin_roi+1], 'EdgeColor', 'b');
    rectangle('position', [xmin_bin_roi+xmin_bin_particle-1, ymin_bin_roi+ymin_bin_particle-1, particle_w, particle_h], 'EdgeColor', 'r');
    
%     saveas(figure(8), strcat(filename, '_2_IDbin.png')); % Sam, 220510   
    saveas(figure(8), strcat(index, filename, '_2_IDbin.png')); % Sam, 220510   
    
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
        
%     my1 = [ 0  0  1  1;
%             1  1  2  2;
%             2  2  4  4;
%             4  4  0  0;
%             0  0 -4 -4;
%            -4 -4 -2 -2;
%            -2 -2 -1 -1;
%            -1 -1 -0 -0; ];
%     mx1 = my1';
%     mx2 = fliplr(mx1);
%     my2 = fliplr(my1);
%     edge_x1 = imfilter(double(IDparticle), mx1, 'symmetric');
%     edge_y1 = imfilter(double(IDparticle), my1, 'symmetric');  
%     edge_x2 = imfilter(double(IDparticle), mx2, 'symmetric');
%     edge_y2 = imfilter(double(IDparticle), my2, 'symmetric');
%     edges = sqrt(edge_x1.^2 + edge_y1.^2 + edge_x2.^2 + edge_y2.^2);
%     figure(99), imshow(edges,[]);
    
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

% calculate attenuation, Sam, 220429
function particle_attenuation_raw = CalAttenuation(IDraw, IDparticle, without_zero)
    % calculate percentiles
    if without_zero == 0
        [pix, ~] = sort(IDparticle(:));
    else
        [pix, ~] = sort(IDparticle(IDparticle>0));
    end
    numPix = length(pix);
    idx05 = round(0.05 * (numPix+1));
    particle_intensity=pix(max(idx05,1));
    
    % calcualte median of entire image (after LSC)
    image_level=median(IDraw(:));
    
    % calculate contrast of particle as the delta/image median (unit is
    % percentage) Negative values showing intensity of particle area is lower than image
    % level
    particle_attenuation_raw=particle_intensity/image_level-1;
end

function [mx, my] = CalCenterOfGravity(A)
x = 1 : size(A, 2); % Columns.
y = 1 : size(A, 1); % Rows.
[X, Y] = meshgrid(x, y);
meanA = mean(A(:));
mx = mean(A(:) .* X(:)) / meanA;
my = mean(A(:) .* Y(:)) / meanA;
end

function std2 = CalStdOfDistWithCenterOfGravity(A)
    [mx, my] = CalCenterOfGravity(A);
    [col, row] = find(A > 0);

    dist2 = zeros(length(col), 1);
    for k = 1:length(col)
        x = col(k);
        y = row(k);
        
        dist2(k, 1) = (mx-x)^2 + (my-y)^2;
    end
    
    std2 = std(dist2);
end

% Sam, 220510
function [col, row] = FindLoctionOfMaxVal(index, filename, filterWidth, I_filtered_c, col, row, feature_roiSize, center_first)
% original location, Sam, 220512
% [col, row] = find((I_filtered_c_cen') == max(I_filtered_c_cen(:)),1);
clf(figure(2));

[rows, cols] = size(I_filtered_c);

filterWidth = 5;

if filterWidth==5 % HF
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
            if similar==0
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
    figure(2),
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

    % replace to optimize location, Sam, 220505
%     fprintf('possible location:\r');
%     fprintf('  ABC : %d, %d\r', col, row);
%     fprintf('  DPBU: %d, %d\r', col_t(area_value_idx(1)), row_t(area_value_idx(1)));
    % can't find full roi, use first one, Sam, 220517
    col = col_t(area_value_idx(search_idx));
    row = row_t(area_value_idx(search_idx));
else
    I_filtered_c_cen_blur = medfilt2(I_filtered_c, [3, 3], 'symmetric');
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
            roi_i = I_filtered_c_cen_blur(roi_y_min:roi_y_max, roi_x_min:roi_x_max); % Sam, 220427
            roi_sum = sum(roi_i(:));
            
            if roi_sum > roi_sum_max
                roi_sum_max = roi_sum;
                col_ = i;
                row_ = j;
            end
        end
    end   
%     fprintf(fid, '%s,%d,%d,%d,%d,', filename, col, row, col_, row_);
    
%     figure(8), subplot(1,3,1), % Sam, 220510
    figure(2),
    annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');
    
    subplot(1,2,1), imshow(I_filtered_c, []), title({'orignial possible', 'location'}, 'FontSize', 8), % Sam, 220512
    rectangle('position', [col-2, row-2, 5, 5], 'EdgeColor', 'r');
%     set(gcf,'position',[800,500,400,400]),
    subplot(1,2,2), imshow(I_filtered_c, []), title({'optimized possible location'}, 'FontSize', 8), % Sam, 220510
%     rectangle('position', [col-2, row-2, 5, 5], 'EdgeColor', 'r');
%    rectangle('position', [col_-4, row_-4, 9, 9], 'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', ':');
    rectangle('position', [col_-4, row_-4, 9, 9], 'EdgeColor', [1 0.6823 0.7882], 'LineWidth', 1);
    
    % replace col and row to col_ and row_, Sam, 220427
%     fprintf('possible location:\r');
%     fprintf('  ABC : %d, %d\r', col, row);
%     fprintf('  DPBU: %d, %d\r', col_, row_);
    col = col_;
    row = row_;
end

saveas(figure(2), strcat(index, filename, '_0_PossLoc.png')); % Sam, 220506
    
end

function IDbin_hist_rm4cor = Rm4Cor(IDbin_hist, cols, rows)
    num = 25;
    ext = 10;
    
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

% new edge filter, Sam, 220516
function [ret, xmin_bin_particle, xmax_bin_particle, ymin_bin_particle, ymax_bin_particle, particle_w, particle_h, ...
          IDbin_o_roi_particle_avg, IDbin_o_roi_border_avg, IDbin_hist_roi_particle_avg, IDbin_hist_roi_border_avg] = ...
    ParticleDetector(index, filename, IDraw, roi_o, IDbin_o, IDbin_hist, IDbin_hist_rm4cor, xmin_bin_roi, xmax_bin_roi, ymin_bin_roi, ymax_bin_roi, roiX_size, roiY_size, feature_roiSize)
    clf(figure(8));
    clf(figure(9));
    figure(9), annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');

    title_font_size = 7;
    
    % initial output data
    xmin_bin_particle = 0;
    ymin_bin_particle = 0;
    xmax_bin_particle = 0;
    ymax_bin_particle = 0;
    particle_w = 0;
    particle_h = 0;
    
    IDbin_o_roi_particle_avg = 0;
    IDbin_o_roi_border_avg = 0;
    IDbin_hist_roi_particle_avg = 0;
    IDbin_hist_roi_border_avg = 0;
    
    ret = 0;
    % return of hough circle, 1: catch some circles, 0: didn't catch any circle, Sam, 220519
    ret_hough_circle = 0;
    
    % return of context of possible object, 1: has something, 0: nothing left
    ret_possible_object = 0;

    % crop mask of remove 4 corner
    IDbin_hist_rm4cor_roi = IDbin_hist_rm4cor(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi);
    IDbin_hist_rm4cor_roi(IDbin_hist_rm4cor_roi>0) = 1;
    
    % crop roi of IDbin_hist
    IDbin_o_roi = IDbin_o(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi); % original IDbin
    IDbin_hist_roi = IDbin_hist(ymin_bin_roi:ymax_bin_roi,xmin_bin_roi:xmax_bin_roi); % enhanced IDbin
    % imshow
    subplot(3, 5, 1), imshow(IDbin_hist_roi, []); title({'roi of enhanced', 'IDbin'}, 'FontSize', title_font_size);

    % get size of roi
    [rows, cols, ~] = size(IDbin_hist_roi);
    
    % transfer to frequency domain, Sam, 220606
% %     H1=fft2(im2double(imfilter(IDbin_o_roi, ones(3, 3)./9, 'replicate'))); %FFT
%     H1=fft2(double(IDbin_o_roi)); %FFT
%     H2=fftshift(H1); %FFT???? 
%     T=log(abs(H2)); %??????
%     subplot(3, 5, 9), imshow(T, []); title('frequency domain');
%     
%     S1=H2; %G??????????
%     [M,N]=size(S1);
%     n1=floor(M/2);
%     n2=floor(N/2);
%     d0=3; %GLPF???d0=5?15?30
%     n=10;
%     for i=1:M     
%         for j=1:N   
%             d=sqrt((i-n1)^2+(j-n2)^2);      
%             
% %             h=exp(-1/2*(d^2/d0^2)); % gaussian low-pass filter
% %             h=1-exp(-1/2*(d^2/d0^2)); % gaussian high-pass filter
% %             h=1/(1 + (d/d0)^(2*n)); % butterworth low-pass filter
%             h=1/(1 + (d0/d)^(2*n)); % butterworth high-pass filter
% 
%             S1(i,j)=h*S1(i,j);
%         end
%     end
%     S2=ifftshift(S1); 
%     S3=real(ifft2(S2));
% %     S3_avg = mean(S3(:));
% %     S3_std = std(S3(:));
% %     S3_bias = 3;
% %     S3_ulm = S3_avg + S3_bias*S3_std;
% %     S3_dlm = S3_avg - S3_bias*S3_std;
% %     S3 = (S3 - S3_dlm) .* (S3_ulm-S3_dlm);
% %     S3(S3<0) = 0;
% %     S3(S3>1) = 1;
%     subplot(3, 5, 10), imshow(S3, []), title('filterd image');

    % edge filter
%     my = [ 1  2  1;
%            2  4  2;
%            0  0  0;
%           -2 -4 -2;
%           -1 -2 -1; ];
%     my = [ 1  2  1;
%            2  3  2;
%            3  4  3;
%            0  0  0;
%           -3 -4 -3;
%           -2 -3 -2;
%           -1 -2 -1; ];
%     mx = my'

% Sam, 220601
%         my1 = [ 0  0  0  1  1;
%                 1  1  0  2  2;
%                 2  2  0  4  4;
%                 4  4  0  8  8;
%                 8  8  0  0  0;
%                 0  0  0 -8 -8;
%                -8 -8  0 -4 -4;
%                -4 -4  0 -2 -2;
%                -2 -2  0 -1 -1;
%                -1 -1  0  0  0;];
%         my1 = [ 0  0  1  1;
%                 1  1  2  2;
%                 2  2  4  4;
%                 4  4  8  8;
%                 8  8  0  0;
%                 0  0 -8 -8;
%                -8 -8 -4 -4;
%                -4 -4 -2 -2;
%                -2 -2 -1 -1;
%                -1 -1  0  0;];
        my1 = [ 0  0  1  1;
                1  1  2  2;
                2  2  4  4;
                4  4  0  0;
                0  0 -4 -4;
               -4 -4 -2 -2;
               -2 -2 -1 -1;
               -1 -1 -0 -0; ];
%         my11 = [ 0  0 -1 -1;
%                 -1 -1 -2 -2;
%                 -2 -2 -4 -4;
%                 -4 -4  0  0;
%                  0  0  4  4;
%                  4  4  2  2;
%                  2  2  1  1;
%                  1  1  0  0; ];
%         my1 = [ 0  0  1  1;
%                 1  1  2  2;
%                 2  2  4  4;
%                 4  4  0  0;
%                 0  0 -4 -4;
%                -4 -4 -2 -2;
%                -2 -2 -1 -1;
%                -1 -1 -0 -0; ];
%         my1 = [ 0  0  0  1  1;
%                 0  0  1  2  2;
%                 1  1  2  8  8;
%                 2  2  8  0  0;
%                 8  8  0 -8 -8;
%                 0  0 -8 -2 -2;
%                -8 -8 -2 -1 -1;
%                -2 -2 -1 -0 -0;
%                -1 -1 -0 -0 -0; ];
%         my1 = [ 0  1;
%                1  2;
%                2  4;
%                4  0;
%                0 -4;
%               -4 -2;
%               -2 -1;
%               -1 -0];
%        my = [  0  0  1;
%                0  1  2;
%                1  2  4;
%                2  4  0;
%                4  0 -4;
%                0 -4 -2;
%               -4 -2 -1;
%               -2 -1 -0;
%               -1 -0 -0;];
%        my = [  0  0;
%                0 -1;
%                1 -2;
%                2 -4;
%                4  0;
%                0  4;
%               -4  2;
%               -2  1;
%               -1  0;];
        mx1 = my1';
%         mx11 = my11';
        mx2 = fliplr(mx1);
%         mx22 = fliplr(mx11);
        my2 = fliplr(my1);
%         my22 = fliplr(my11);
%     mx = [2 0 -2]; my = mx';
%     mx = [1 1 2 0 -2 -1 -1]; my = mx';
%     mx = [1 2 2 0 -2 -2 -1]; my = mx'; % Sam 220524
%     mx = [2 2 4 0 -4 -2 -2]; my = mx';
%     mx = [1 2 8 0 -8 -2 -1]; my = mx';
%     mx = [1 1 8 0 -8 -1 -1]; my = mx';
%    mx = [1 8 0 -8 -1]; my = mx';
%     edge_x1 = imfilter(double(IDbin_hist_roi), mx1, 'replicate');
%     edge_y1 = imfilter(double(IDbin_hist_roi), my1, 'replicate');
%     edge_x2 = imfilter(double(IDbin_hist_roi), mx2, 'replicate');
%     edge_y2 = imfilter(double(IDbin_hist_roi), my2, 'replicate');
    edge_x1 = imfilter(double(IDbin_hist_roi), mx1, 'symmetric');
%     edge_x11 = imfilter(double(IDbin_hist_roi), mx11, 'symmetric');
    edge_y1 = imfilter(double(IDbin_hist_roi), my1, 'symmetric');
%     edge_y11 = imfilter(double(IDbin_hist_roi), my11, 'symmetric');    
    edge_x2 = imfilter(double(IDbin_hist_roi), mx2, 'symmetric');
%     edge_x22 = imfilter(double(IDbin_hist_roi), mx22, 'symmetric');
    edge_y2 = imfilter(double(IDbin_hist_roi), my2, 'symmetric');
%     edge_y22 = imfilter(double(IDbin_hist_roi), my22, 'symmetric');
%     edge_x1 = imfilter(double(S3), mx1, 'replicate');
%     edge_y1 = imfilter(double(S3), my1, 'replicate');
%     edge_x2 = imfilter(double(S3), mx2, 'replicate');
%     edge_y2 = imfilter(double(S3), my2, 'replicate');
    edges = sqrt(edge_x1.^2 + edge_y1.^2 + edge_x2.^2 + edge_y2.^2);
    edges = edges.*IDbin_hist_rm4cor_roi;
%     edges = sqrt(max(edge_y1.^2, edge_y2.^2) + max(edge_x1.^2, edge_x2.^2));
%     edges_minx = min(max(edge_x1, edge_x11), min(edge_x2, edge_x22));
%     edges_miny = min(max(edge_y1, edge_y11), min(edge_y2, edge_y22));
%     edges_min = min(edges_minx, edges_miny);
%     edges_maxx = max(max(edge_x1, edge_x11), max(edge_x2, edge_x22));
%     edges_maxy = max(max(edge_y1, edge_y11), max(edge_y2, edge_y22));
%     edges_max = max(edges_maxx, edges_maxy);
%     edges = edges_max.*IDbin_hist_rm4cor_roi + edges_min.*~IDbin_hist_rm4cor_roi;
%     edges(edges<0) = 0;
    
    edge_avg = mean(edges(:));
    edge_std = std(edges(:));
    fprintf('avg/std of edge(IDbin_hist): %f / %f\n', edge_avg, edge_std);
    
%     mask = ~IDbin_hist_rm4cor_roi;
%     edge_cor = min(edge_x .* mask, edge_x .* mask);
    % imshow
%     subplot(3, 5, 2), imshow(edge, []); title({'IDbin_hist edge filter'}, 'FontSize', title_font_size); 
    subplot(3, 5, 2), imshow(edges, []); title({'edge filtered', 'of enhanced IDbin'}, 'FontSize', title_font_size); 
%     subplot(3, 5, 5), imshow(edge+edge2, []); title({'IDbin_hist edge filter'}, 'FontSize', title_font_size); 
%     
%     bitDepth=8;
%     edges_nor = edges./max(edges(:)) * (2^bitDepth - 1); % Sam, 220425
%     edges_bin = ~thresholdROI(edges_nor, bitDepth, 1);
%     subplot(3, 5, 3), imshow(edges_bin, []);

    % preserve top 5% data
%     [edge_pix, ~] = sort(edge(:), 'descend');    
%     pre_ratio = 0.05; % Sam, 220510
%     pre_val = edge_pix( round(pre_ratio*length(edge_pix)) ); % preserve top 10% pixel value
    pre_ratio = 1.70;
    pre_val = edge_avg + pre_ratio*edge_std;
    fprintf('pre_val: %.2f\n', pre_val);
    pre = edges;
    pre(pre<pre_val) = 0;
    pre(pre>0) = 1;
    % imshow
%     subplot(3, 5, 3), imshow(pre, []); title({'preserve top', sprintf('%.2f%% pixel', pre_ratio*100)}, 'FontSize', title_font_size);
    subplot(3, 5, 3), imshow(pre, []); title({'preserve', sprintf('avg + %.2f% x std pixels', pre_ratio)}, 'FontSize', title_font_size);
    
    % remove small object
    object_small_th = 7;
    pre_rm = bwareaopen(pre, object_small_th);
    % imshow
    subplot(3, 5, 4), imshow(pre_rm, []); title({'remove object', sprintf('smaller than %d',object_small_th)}, 'FontSize', title_font_size);
    
    % morphology operator
    nhood3 = [0 1 0;
              1 1 1;
              0 1 0 ];
%     pre_rm_close = imclose(pre_rm, nhood3);
    pre_rm_close = pre_rm;
%     pre_rm_close = pre;
    
    % get rough coordinates of particle, Sam, 220519
    [y, x] = find(pre_rm_close==1);
    xmin_bin_particle = min(x);
    xmax_bin_particle = max(x);
    ymin_bin_particle = min(y);
    ymax_bin_particle = max(y);
    particle_w = xmax_bin_particle - xmin_bin_particle;
    particle_h = ymax_bin_particle - ymin_bin_particle;
        
    pre_rm_close_h = Thinning(pre_rm_close); % thinning, Sam, 220606
    pre_rm_close_h = pre_rm_close_h(3:rows+2, 3:cols+2);
    subplot(3, 5, 5), imshow(pre_rm_close_h, []); title({'skeletonize of', 'roi'}, 'FontSize', title_font_size);
    
    % extract object with maximum ares
    pre_rm_lb = bwlabel(pre_rm_close_h);
    pre_rm_sts = regionprops(logical(pre_rm_lb), 'Area', 'Image', 'BoundingBox');
    lb_cnt = length(pre_rm_sts); % Sam, 220604
    [~, pre_rm_idx] = sort([pre_rm_sts.Area], 'descend');
    if isempty(pre_rm_idx)
        return;
%     elseif length(pre_rm_idx) > 3 % Sam, 220605
%         return;
    end
%     lb_size = 1;
%     for i=1:lb_cnt
%         if pre_rm_sts(i).BoundingBox(3)<(cols/2.0) || pre_rm_sts(i).BoundingBox(4)<(rows/2.0)
%             lb_size = 0;
%             break;
%         end
%     end
    
    % use hough circle to locate possible particle, Sam, 220518
%     HoughCircle2(pre_rm_close_ext, [6, 15]);
%     [c, r] = imfindcircles(pre_rm_close,[5, 20]);
%     [c, r] = imfindcircles(pre_rm_close_ext, [6, 18], 'Method', 'TwoStage', 'Sensitivity', 0.93); % Sam, 220528
%     c = c - ext;
    ext = ceil(max(rows, cols) / 3.0); % extend border, Sam, 220603
    pre_rm_close_ext = zeros(rows+ext*2, cols+ext*2);
    pre_rm_close_ext(ext+1:ext+rows, ext+1:ext+cols) = pre_rm_close_h;
%     mask = zeros(rows+ext*2, cols+ext*2);
%     mask(ext+1:ext+rows, ext+1:ext+cols) = 1;

    [cx, cy, rr, cnt, er, sr, ar] = HoughCircle2(pre_rm_close_ext, 5, 13, 1, 3, 0.85); % Sam, 220608
%     [cx, cy, rr, cnt, cm] = HoughCircle3(pre_rm_close_ext, 5, 15, 1, 2, 0.85);
%     fprintf('(%d): cm: %f\n', 1, cm);
%     [cx, cy, rr, cnt, prob] = HoughCircle(pre_rm_close_ext, mask, 5, 15, 1, 2, 0.85);
%     [~, prob_idx] = sort(prob, 'descend');
    [~, er_idx] = sort(er); % Sam, 220608

    % draw circle as mask, Sam, 220610
    mask_particle = zeros(size(pre_rm_close_ext));
    mask_particle = drawellipse(cy(er_idx(1)), cx(er_idx(1)), rr(er_idx(1))-2, rr(er_idx(1))-2, 0, mask_particle, 1);
    mask_particle = mask_particle(ext+1:ext+rows, ext+1:ext+cols);
    mask_background = zeros(size(pre_rm_close_ext));
    mask_background = ~drawellipse(cy(er_idx(1)), cx(er_idx(1)), rr(er_idx(1)), rr(er_idx(1)), 0, mask_background, 1);
    mask_background = mask_background(ext+1:ext+rows, ext+1:ext+cols);
    cx = cx-ext;
    cy = cy-ext;
    fprintf('(%d): rr, cnt, prob: %d, %d, %f\n', er_idx(1), rr(er_idx(1)), cnt(er_idx(1)), er(er_idx(1)) );

    IDbin_hist_roi_particle   = IDbin_hist_roi .* mask_particle .* IDbin_hist_rm4cor_roi;
    IDbin_hist_roi_background = IDbin_hist_roi .* mask_background .* IDbin_hist_rm4cor_roi;
    IDbin_hist_roi_particle_avg = mean(IDbin_hist_roi_particle(IDbin_hist_roi_particle>0));
    IDbin_hist_roi_background_avg = mean(IDbin_hist_roi_background(IDbin_hist_roi_background>0));
    fprintf('ratio of particle and background: %.4f\n', IDbin_hist_roi_particle_avg/IDbin_hist_roi_background_avg);
    
    ret_hough_circle = 1;
    
    xmin_bin_particle = max(1,cx(er_idx(1))-rr(er_idx(1)));
    xmax_bin_particle = min(cx(er_idx(1))+rr(er_idx(1)),cols);
    ymin_bin_particle = max(1,cy(er_idx(1))-rr(er_idx(1)));
    ymax_bin_particle = min(cy(er_idx(1))+rr(er_idx(1)),rows);
    particle_w = xmax_bin_particle - xmin_bin_particle;
    particle_h = ymax_bin_particle - ymin_bin_particle;
%     if ~isempty(c)
%         ret_hough_circle = 1;
%     end
    % imshow
%     if ret_hough_circle
%         subplot(3, 5, 6), imshow(pre_rm_close, []); title({'close operation', 'with hough circle'}, 'FontSize', title_font_size);
        subplot(3, 5, 6), imshow(pre_rm_close_h, []); title({'hough circle', 'detection'}, 'FontSize', title_font_size);
        
        th = 0:pi/50:2*pi;
        hold on;
%         for p = 1:length(r)
%         for p = 1:length(rr)
% %             x = r(p) * cos(th) + c(p,1);
% %             y = r(p) * sin(th) + c(p,2);
%             x = rr(p) * cos(th) + cx(p);
%             y = rr(p) * sin(th) + cy(p);
%             plot(x, y, 'Color', 'm');
%         end
        % final circle
        x = rr(er_idx(1)) * cos(th) + cx(er_idx(1));
        y = rr(er_idx(1)) * sin(th) + cy(er_idx(1));
        plot(x, y, 'Color', 'm');
        
        rectangle('position', [cx(er_idx(1))-rr(er_idx(1))-10, cy(er_idx(1))-rr(er_idx(1))-10, 2*rr(er_idx(1))+20, 2*rr(er_idx(1))+20], 'EdgeColor', 'm');
        hold off;
%     end



    sft = 5;
    xmin_zoom = max(1,cx(er_idx(1))-rr(er_idx(1))-sft);
    ymin_zoom = max(1,cy(er_idx(1))-rr(er_idx(1))-sft);
    xmax_zoom = min(cx(er_idx(1))+rr(er_idx(1))+sft,cols);
    ymax_zoom = min(cy(er_idx(1))+rr(er_idx(1))+sft,rows);
    zoom_w = xmax_zoom-xmin_zoom+1;
    zoom_h = ymax_zoom-ymin_zoom+1;
    IDbin_o_roi_zoom = IDbin_o_roi(ymin_zoom:ymax_zoom, xmin_zoom:xmax_zoom);
    IDbin_hist_roi_zoom = IDbin_hist_roi(ymin_zoom:ymax_zoom, xmin_zoom:xmax_zoom);

    
%       mx = [ 3  10  3;
%              0   0  0;
%             -3 -10 -3]; my = mx';
%     mx = [ 3  10   0;
%           10   0  -10;
%            0 -10  -3];
%     my = [  0  10   3;
%           -10   0  10;
%            -3 -10   0];
    % Sam, 220614
    mx = [ 1  2  0;
           2  0 -2;
           0 -2 -1];
    my = [ 0  2  1;
          -2  0  2;
          -1 -2  0];
%     mx = [4 4 8 16 0 -16 -8 -4 -4]; my = mx';
%     mx = [1 2 0 -2 -1 -1]; my = mx';
%     mx = [1 0 -1;
%           2 0 -2;
%           1 0 -1]; my = mx';
%     my = [ 1  2  1;
%            2  4  2;
%            0  0  0;
%           -2 -4 -2;
%           -1 -2 -1; ];
%     mx = my';
%     my = [ 1  2  1;
%            2  3  2;
%            3  4  3;
%            0  0  0;
%           -3 -4 -3;
%           -2 -3 -2;
%           -1 -2 -1; ];
%     my = [ 1  1  1;
%            2  2  2;
%            4  4  4;
%            0  0  0;
%           -4 -4 -4;
%           -2 -2 -2;
%           -1 -1 -1; ];
%     mx = my';
    edge_x = imfilter(double(IDbin_hist_roi_zoom), mx, 'symmetric');
    edge_y = imfilter(double(IDbin_hist_roi_zoom), my, 'symmetric');
    edges_zoom = sqrt(edge_x.^2 + edge_y.^2);
    
    IDbin_hist_rm4cor_roi_zoom = IDbin_hist_rm4cor_roi(ymin_zoom:ymax_zoom, xmin_zoom:xmax_zoom);
    IDbin_hist_rm4cor_roi_zoom(IDbin_hist_rm4cor_roi_zoom>0) = 1;

    edges_zoom = edges_zoom .* IDbin_hist_rm4cor_roi_zoom;
    
%     subplot(3,5,11), imshow(0.5*edges_zoom.*IDbin_hist_roi_zoom,[]);
    
    edges_zoom = edges_zoom + edges(ymin_zoom:ymax_zoom, xmin_zoom:xmax_zoom); % Sam, 220614
    subplot(3,5,7), imshow(edges_zoom,[]); title({'extract roi of', 'result of hough circle'}, 'FontSize', title_font_size);
    
    % binearize, Sam, 220608
    bitDepth=8;
    edges_zoom_nor = edges_zoom ./ max(edges_zoom(:)) * (2^bitDepth - 1); % Sam, 220425
    edges_zoom_nor_bin = ~thresholdROI(edges_zoom_nor, bitDepth, 1); 
    subplot(3,5,8), imshow(edges_zoom_nor_bin,[]); title({'binarize of', 'extract roi'}, 'FontSize', title_font_size);

    % remove small object, Sam, 220615
    object_small_th = 7;
    edges_zoom_nor_bin = bwareaopen(edges_zoom_nor_bin, object_small_th);
    subplot(3,5,9), imshow(edges_zoom_nor_bin,[]); title({'remove object', sprintf('smaller than %d',object_small_th)}, 'FontSize', title_font_size);
    % imshow
    
    % skeletonize
    edges_zoom_nor_bin = Thinning(edges_zoom_nor_bin); % thinning, Sam, 220606
    edges_zoom_nor_bin = edges_zoom_nor_bin(3:zoom_h+2, 3:zoom_w+2);
    subplot(3,5,10), imshow(edges_zoom_nor_bin,[]); title({'skeletonize of', 'extract roi'}, 'FontSize', title_font_size);
    
    % hough circle
    ext = ceil(max(zoom_h, zoom_w) / 2.0); % extend border, Sam, 220603
    edges_zoom_ext = zeros(zoom_h+ext*2, zoom_w+ext*2);
    edges_zoom_ext(ext+1:ext+zoom_h, ext+1:ext+zoom_w) = edges_zoom_nor_bin;
    [cx2, cy2, rr2, cnt2, er2, sr2, ar2] = HoughCircle2(edges_zoom_ext, 5, 13, 1, 3, 0.85); % Sam, 220608
    [~, er_idx2] = sort(er2); % Sam, 220608
%     [edges_zoom_nor_bin_row, edges_zoom_nor_bin_col] = find(edges_zoom_nor_bin);
%     fprintf('cnt2: %d, er: %.3f\n', cnt2(er_idx(1)), er(er_idx(1)) );
    fprintf('(%d): rr2, cnt2, er2: %d, %d, %f\n', er_idx2(1), rr2(er_idx2(1)), cnt2(er_idx2(1)), er2(er_idx2(1)) );
    % replace result of hough circle, Sam, 220615
%     er_idx = er_idx2;
%     cx = xmin_zoom - 1 + cx2 - ext;
%     cy = ymin_zoom - 1 + cy2 - ext;
%     rr = rr2;
%     cnt = cnt2;
%     er = er2;
%     sr = sr2;
%     ar = ar2;
%     fprintf('cnt2: %d, er2: %.3f\n', cnt2(er2_idx(1)), er2(er2_idx(1)) );
    
    % genertae circle bounding box, Sam, 220518
    if ret_hough_circle
        mask_particle_border = false(rows, cols);
        mask_particle_inside = false(rows, cols);
        % Sam, 220615
        rr_ou = rr(er_idx(1))+4;
        rr_in = rr(er_idx(1))-4;
        for th = 0:pi/50:2*pi
%             x = round(r(1)*cos(th)+c(1,1));
%             y = round(r(1)*sin(th)+c(1,2));
%             r = rr(er_idx(1))+4;
            r = rr_ou;
            x = round(r*cos(th)+cx(er_idx(1)));
            y = round(r*sin(th)+cy(er_idx(1)));
%             x = min(max(1,x),cols);
%             y = min(max(1,y),rows);
            if ~(x<1 || x>cols || y<1 || y>rows) % Sam, 220604
                mask_particle_border(y,x) = 1;
            end
            
%             r = rr(er_idx(1))-4;
            r = rr_in;
            x = round(r*cos(th)+cx(er_idx(1)));
            y = round(r*sin(th)+cy(er_idx(1))); 
            if ~(x<1 || x>cols || y<1 || y>rows) % Sam, 220604
                mask_particle_inside(y,x) = 1;
            end
        end
        
        % morphologh operator, Sam, 220518
%     nhood5 = ones(5,5);
        mask_particle_border_dilate = logical(imdilate(mask_particle_border, nhood3));
        mask_particle_inside_dilate = logical(imdilate(mask_particle_inside, nhood3));
%         border_size = 1.4; % Sam, 220524
%         circle_bin_dilate = imresize(circle_bin, border_size, 'method', 'nearest');
%     circle_bbox_x = floor(c(1,1)-r(1));
%     circle_bbox_y = floor(c(1,2)-r(1));
%     circle_bbox_w = ceil(r(1)*2);
%     circle_bbox_h = ceil(r(1)*2);
%     circle_bbox = [circle_bbox_x, circle_bbox_y, circle_bbox_w, circle_bbox_h];
        % imshow
%         subplot(3, 5, 7), imshow(mask_particle_border_dilate, []); title({'use dilate of', 'Hough circle', 'as mask(border)'}, 'FontSize', title_font_size);
%         subplot(3, 5, 8), imshow(mask_particle_inside_dilate, []); title({'use dilate of', 'Hough circle', 'as mask(inside)'}, 'FontSize', title_font_size);
        
        th = 0:pi/50:2*pi;
        subplot(3, 5, 11), imshow(IDbin_hist_roi,[]);
        hold on;
        x = (rr_ou+1) * cos(th) + cx(er_idx(1));
        y = (rr_ou+1) * sin(th) + cy(er_idx(1));
        plot(x,y,'b');
        x = (rr_ou-1) * cos(th) + cx(er_idx(1));
        y = (rr_ou-1) * sin(th) + cy(er_idx(1));
        plot(x,y,'b');
        x = (rr_in+1) * cos(th) + cx(er_idx(1));
        y = (rr_in+1) * sin(th) + cy(er_idx(1));
        plot(x,y,'r');
        x = (rr_in-1) * cos(th) + cx(er_idx(1));
        y = (rr_in-1) * sin(th) + cy(er_idx(1));
        plot(x,y,'r');
        % draw hough circle, Sam, 220615
        x = rr(er_idx(1)) * cos(th) + cx(er_idx(1));
        y = rr(er_idx(1)) * sin(th) + cy(er_idx(1));
        plot(x, y, 'Color', 'c');
        hold off;
        title({'use dilate of', 'Hough circle', 'as mask'}, 'FontSize', title_font_size);
    
%     m2 = [ 1  2  1;
%            0  0  0;
%           -1 -2 -1];
%     m3 = [ 1  0 -1;
%            2  0 -2;
%            1  0 -1];
%     m4 = [ 1  2  0;
%            2  0 -2;
%            0 -2 -1];
%     m5 = [ 0  2  1;
%           -2  0  2;
%           -1 -2  0];
%     edge2 = imfilter(double(IDbin_hist_roi_zoom), m2, 'symmetric');
%     edge3 = imfilter(double(IDbin_hist_roi_zoom), m3, 'symmetric');
%     edge4 = imfilter(double(IDbin_hist_roi_zoom), m4, 'symmetric');
%     edge5 = imfilter(double(IDbin_hist_roi_zoom), m5, 'symmetric');
%     edgeN = sqrt(edge2.^2+edge3.^2+edge4.^2+edge5.^2);
%     subplot(3,5,11), imshow(edgeN,[]);
        

        IDbin_o_roi_border = IDbin_o_roi .* mask_particle_border_dilate;
        IDbin_o_roi_inside = IDbin_o_roi .* mask_particle_inside_dilate;
        IDbin_o_roi_border_avg = mean(IDbin_o_roi_border(IDbin_o_roi_border>0));
        IDbin_o_roi_inside_avg = mean(IDbin_o_roi_inside(IDbin_o_roi_inside>0));
%         IDbin_o_roi_sum = IDbin_o_roi_border + IDbin_o_roi_inside; % for debug
        IDbin_o_roi_border_inside = IDbin_o_roi_inside_avg / IDbin_o_roi_border_avg;
        fprintf('ratio of inside and border(original): %.3f\n', IDbin_o_roi_border_inside);
        
        IDbin_hist_roi_border = IDbin_hist_roi .* mask_particle_border_dilate;
        IDbin_hist_roi_inside = IDbin_hist_roi .* mask_particle_inside_dilate;
        IDbin_hist_roi_border_avg = mean(IDbin_hist_roi_border(IDbin_hist_roi_border>0));
        IDbin_hist_roi_inside_avg = mean(IDbin_hist_roi_inside(IDbin_hist_roi_inside>0));
%         IDbin_hist_roi_sum = IDbin_hist_roi_border + IDbin_hist_roi_inside; % for debug
        IDbin_hist_roi_border_inside = IDbin_hist_roi_inside_avg / IDbin_hist_roi_border_avg;
        fprintf('ratio of inside and border(enhanced): %.3f\n', IDbin_hist_roi_border_inside);
    end   
    
    % calculate attenuation, Sam, 220516
    % convert x address from binned image coordinates to raw image coordinates
    xmin_raw_particle=(xmin_bin_roi+xmin_bin_particle-1)*roiX_size-roiX_size/2;
%     xmax_raw_particle=(xmax_bin_roi+xmax_bin_particle-1)*roiX_size-roiX_size/2;
    particle_w_raw=particle_w*roiX_size;
    % convert y address from binned image coordinates to raw image coordinates
    ymin_raw_particle=(ymin_bin_roi+ymin_bin_particle-1)*roiY_size-roiY_size/2;
%     ymax_raw_particle=(ymax_bin_roi+ymax_bin_particle-1)*roiY_size-roiY_size/2;
    particle_h_raw=particle_h*roiY_size;
    % calculate particle intensity by taking the median of half height, half width aread in the center of ROI
    IDparticle = IDraw(round(ymin_raw_particle):round(ymin_raw_particle+particle_h_raw),round(xmin_raw_particle):round(xmin_raw_particle+particle_w_raw));
    attenuation = CalAttenuation(IDraw, IDparticle, 0);
    fprintf('attenuation of object with maximum area: %.4f\n', attenuation);
    
    % output
    if ret_hough_circle == 1 % Sam, 220604
%     if ret_hough_circle == 1 && ret_possible_object == 1
%         axis_th = 5;
%         dist_avg_th = 1.5;
        er_th = 1.5; % Sam, 220608
        sr_th = 1.5; % Sam, 220608
%         prob_th = 0.3;
        IDbin_o_roi_border_inside_th = -1;
        IDbin_hist_roi_border_inside_th = 0.95;
        attenuation_th = -1; %attenuation = 0.02;
        lb_cnt_th = 5;
        cnt_th = 45;
        ar_th = 0.6; % Sam, 220608
%         ratio_par_bor_th = 0.996; % Sam, 220524
        cond0 = 0; cond0_r = 'X'; % mean of error of hough circle, Sam, 220608
        cond1 = 0; cond1_r = 'X'; % std of error of hough circle
        cond2 = 0; cond2_r = 'X'; % ratio of mean between inside and border of IDbin_o
        cond3 = 0; cond3_r = 'X'; % ratio of mean between inside and border of IDbin_hist
        cond4 = 0; cond4_r = 'X'; % attenuation
        cond5 = 0; cond5_r = 'X'; % appear count in hough space
        cond6 = 0; cond6_r = 'X'; % number of count of hough circle
        cond7 = 0; cond7_r = 'X'; % cover rate of angle, Sam, 220608
        
        if er(er_idx(1)) <= er_th || er_th == -1
            cond0 = 1;
            cond0_r = 'O';
            if er_th == -1
                cond0_r = '..';
            end
        end
        
        if sr(er_idx(1)) <= sr_th || sr_th == -1
            cond1 = 1;
            cond1_r = 'O';
            if sr_th == -1
                cond1_r = '..';
            end
        end

        if IDbin_o_roi_border_inside < IDbin_o_roi_border_inside_th || IDbin_o_roi_border_inside_th == -1
            cond2 = 1;
            cond2_r = 'O';
            if IDbin_o_roi_border_inside_th == -1
                cond2_r = '..';
            end
        end
        
        if IDbin_hist_roi_border_inside < IDbin_hist_roi_border_inside_th || IDbin_hist_roi_border_inside_th == -1
            cond3 = 1;
            cond3_r = 'O';
            if IDbin_hist_roi_border_inside_th == -1
                cond3_r = '..';
            end
        end         

        if abs(attenuation) < attenuation_th || attenuation_th == -1
            cond4 = 1;
            cond4_r = 'O';
            if attenuation_th == -1
                cond4_r = '..';
            end
        end
        
        if lb_cnt <= lb_cnt_th || lb_cnt_th == -1
            cond5 = 1;
            cond5_r = 'O';
            if lb_cnt_th == -1
                cond5_r = '..';
            end
        end
        
        if cnt(er_idx(1)) > cnt_th || cnt_th == -1
            cond6 = 1;
            cond6_r = 'O';
            if cnt_th == -1
                cond6_r = '..';
            end
        end
        
        if ar(er_idx(1)) > ar_th || ar_th == -1
            cond7 = 1;
            cond7_r = 'O';
            if ar_th == -1
                cond7_r = '..';
            end
        end
        
%         % Sam, 220524
%         ratio_par_bor = IDbin_hist_roi_particle_avg / IDbin_hist_roi_border_avg;
%         if ratio_par_bor < ratio_par_bor_th || ratio_par_bor_th == -1
%             cond4 = 1;
%             cond4_r = 'O';
%             if ratio_par_bor_th == -1
%                 cond4_r = '..';
%             end
%         end

        if cond0 + cond1 + cond2 + cond3 + cond4 + cond5 + cond6 + cond7 == 8
%         if cond1 + cond2 + cond3 == 3
            ret = 1;
        end

        cond0_s = sprintf('(%s) mean of error of hough circle: %.3f (%.3f)', cond0_r, er(er_idx(1)), er_th);
        cond1_s = sprintf('(%s) std of error of hough circle : %.3f (%.3f)', cond1_r, sr(er_idx(1)), sr_th);
        cond2_s = sprintf('(%s) ratio of mean between inside and border of IDbin(o)   : %.5f (%.5f)', cond2_r, IDbin_o_roi_border_inside, IDbin_o_roi_border_inside_th);
        cond3_s = sprintf('(%s) ratio of mean between inside and border of IDbin(hist): %.5f (%.5f)', cond3_r, IDbin_hist_roi_border_inside, IDbin_hist_roi_border_inside_th);
        cond4_s = sprintf('(%s) attenuation: %.5f (%.5f)', cond4_r, attenuation, attenuation_th);
        cond5_s = sprintf('(%s) appear count in hough space: %d (%d)', cond5_r, lb_cnt, lb_cnt_th);
        cond6_s = sprintf('(%s) number of count of hough circle: %d (%d)', cond6_r, cnt(er_idx(1)), cnt_th);
        cond7_s = sprintf('(%s) cover rate of angle of hough circle: %.3f (%.3f)', cond7_r, ar(er_idx(1)), ar_th);
        
%         cond4_s = sprintf('(%s) ratio of mean between particle and its border: %.4f (%.4f)', cond4_r, ratio_par_bor, ratio_par_bor_th);
        if ret == 1
            ret_s = sprintf('=> it''s a particle' );
            annotation('textbox', [.29 0.15 0.7 0.2], 'String', {ret_s, cond0_s, cond1_s, cond2_s, cond3_s, cond4_s, cond5_s, cond6_s, cond7_s}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'r');
        else
            ret_s = sprintf('=> it''s noise' );
            annotation('textbox', [.29 0.15 0.7 0.2], 'String', {ret_s, cond0_s, cond1_s, cond2_s, cond3_s, cond4_s, cond5_s, cond6_s, cond7_s}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'b');
        end

        % Sam, 220510
        if ret == 1
            rec_c = 'r';
        else
            rec_c = 'c';
        end
    else
        rec_c = 'c';
    end
    
%     if ret == 1
    figure(8),
    annotation('textbox', [0 0.815 0.7 0.2], 'String', {strrep(filename,'_','\_')}, 'LineStyle', 'none', 'FontSize', 8, 'Color', 'k');
    
    subplot(1,2,1), imshow(IDbin_o,[]), title({'roi of original IDbin'}, 'FontSize', 8);
    rectangle('position', [xmin_bin_roi, ymin_bin_roi, xmax_bin_roi-xmin_bin_roi+1, ymax_bin_roi-ymin_bin_roi+1], 'EdgeColor', 'b');
    rectangle('position', [xmin_bin_roi+xmin_bin_particle-1, ymin_bin_roi+ymin_bin_particle-1, particle_w, particle_h], 'EdgeColor', rec_c);

    subplot(1,2,2), imshow(IDbin_hist), title({'roi of enhanced IDbin'}, 'FontSize', 8);
    rectangle('position', [xmin_bin_roi, ymin_bin_roi, xmax_bin_roi-xmin_bin_roi+1, ymax_bin_roi-ymin_bin_roi+1], 'EdgeColor', 'b');
    rectangle('position', [xmin_bin_roi+xmin_bin_particle-1, ymin_bin_roi+ymin_bin_particle-1, particle_w, particle_h], 'EdgeColor', rec_c);
    
    saveas(figure(9), strcat(index, filename, '_1_particle_detector.png')); % Sam, 220510
    saveas(figure(8), strcat(index, filename, '_2_IDbin.png')); % Sam, 220510
%     end
end

function [cx, cy, rr, cnt, prob] = HoughCircle(img, mask, minr, maxr, stepr, stepa, percent)
r = ceil((maxr-minr)/stepr); % number of radius need to be calculate
t = ceil(360.0/stepa); % number of theta need to be calculate
[rows, cols] = size(img);
hspace = zeros(rows, cols, r);

[row, col] = find(img);
num = length(row);

% calculate hough transform
% a = x - r * cos(a)
% b = y - r * sin(a)
for i=1:num
    for j=1:r
        for k=1:t
            aa = round(row(i) - (minr+(j-1)*stepr) * cosd(k*stepa));
            bb = round(col(i) - (minr+(j-1)*stepr) * sind(k*stepa));
            if aa>0 && aa<=rows && bb>0 && bb<=cols % must inside image
                hspace(aa,bb,j) = hspace(aa,bb,j)+1;
            end
        end
    end
end

% find maximum count in hough space
%par = max(max(max(hspace)));
par = max(hspace(:));
par_th = par*percent;
[m2, n2, r2] = size(hspace);

cx = []; cy = []; rr = [];
cnt = [];
% num = m2 * n2 * r2;

for i=1:m2
    for j=1:n2
        for k=1:r2
            if hspace(i,j,k) >= par_th
                cx = [cx j]; cy = [cy i];
                rr = [rr minr+k*stepr];
                cnt = [cnt hspace(i,j,k)];
            end
        end
    end
end

prob = zeros(1,length(rr));
% figure(4);
for i=1:length(rr)
    rc = ceil(rr(i));
    rc_ext = rc + 7;
    xs = max(1,cx(i) - rc_ext); xe = min(cx(i) + rc_ext,cols);
    ys = max(1,cy(i) - rc_ext); ye = min(cy(i) + rc_ext,rows);
    img_crop = img(ys:ye, xs:xe);
%     subplot(1,2,1); imshow(img_crop);
    
    img_fit = zeros(rows, cols);
    for a=0:pi/50:2*pi
        x = round(cx(i) - rr(i)*cos(a));
        y = round(cy(i) - rr(i)*sin(a));
        if x<1 || x>cols || y<1 || y>rows % Sam, 220604
            continue;
        end
        img_fit(y,x) = 1;
    end
    nhood = [0 1 0;
             1 1 1;
             0 1 0;];
%     img_fit_dilate = imdilate(img_fit, nhood);
%     img_fit_dilate_mask = img_fit_dilate .* mask;
%     img_fit_crop = img_fit_dilate_mask(ys:ye, xs:xe);
    img_fit_crop = img_fit(ys:ye, xs:xe);
%     subplot(1,2,2); imshow(img_crop + 0.5 * img_fit_crop);
    
    coeff = corr2(img_crop, img_fit_crop);
    prob(i) = power(coeff, 2);
%     fprintf('rsq: %f\n', rsq);
end

end

% use distance with circle to calculate similarity, Sam, 220608
% cx: x, center of circle
% cy: y, center of circle
% rr: radius of circle
% cnt: appear count in hough space
% er: error of points with hough circle
% sr: std of points with hough circle
% ar: cover rate of angle of points
function [cx, cy, rr, cnt, er, sr, ar] = HoughCircle2(img, minr, maxr, stepr, stepa, percent)
tic
r = ceil((maxr-minr)/stepr); % number of radius need to be calculate
inters_t = ceil(360.0/stepa); % number of theta need to be calculate
[rows, cols] = size(img);
hspace = zeros(rows, cols, r);

[row, col] = find(img);
num = length(row);

% calculate hough transform
% a = x - r * cos(a)
% b = y - r * sin(a)
for i=1:num
    for j=1:r
        for k=1:inters_t
            aa = round(row(i) - (minr+(j-1)*stepr) * cosd(k*stepa));
            bb = round(col(i) - (minr+(j-1)*stepr) * sind(k*stepa));
            if aa>0 && aa<=rows && bb>0 && bb<=cols % must inside image
                hspace(aa,bb,j) = hspace(aa,bb,j)+1;
            end
        end
    end
end

% find maximum count in hough space
par = max(hspace(:));
par_th = par*percent;
[m2, n2, r2] = size(hspace);
max_num = m2*n2*r2;

% cx = []; cy = []; rr = [];
% cnt = [];
cx_ = zeros(1,max_num); cy_ = zeros(1,max_num);
rr_ = zeros(1,max_num);
cnt_ = zeros(1,max_num);

idx = 0;
for i=1:m2
    for j=1:n2
        for k=1:r2
            if hspace(i,j,k) >= par_th
%                 cx = [cx, j]; cy = [cy, i];
%                 rr = [rr, minr+k*stepr];
%                 cnt = [cnt, hspace(i,j,k)];
                idx = idx + 1;
                cx_(idx) = j; cy_(idx) = i;
                rr_(idx) = minr+k*stepr;
                cnt_(idx) = hspace(i,j,k);
            end
        end
    end
end
cx = cx_(1:idx); cy = cy_(1:idx);
rr = rr_(1:idx);
cnt = cnt_(1:idx);
cnt_num = idx;

% generate mask
fit = zeros(rows, cols);
for i=1:cnt_num
    for inters_t = 0:pi/36:2*pi
%     for inters_t = 0:5:360
        x = round( rr(i) * cos(inters_t) + cx(i) );
        y = round( rr(i) * sin(inters_t) + cy(i) );
        if x<1 || x>cols || y<1 || y>rows
            continue;
        end
        fit(y,x) = 1;
    end
end

% dilate the mask
nhood = [0 1 0;
         1 1 1;
         0 1 0];
fit_dilate = imdilate(fit, nhood);

% do connected component, Sam, 220606
lb = bwlabel(img);
sts = regionprops(logical(lb), 'Area', 'Image', 'BoundingBox');
if isempty(sts)
    return;
end

% get all object with intersection
inters = false(rows, cols);
cc_idx = [];
for i=1:length(sts)
    cc = ismember(lb, i); % object with maximum area
    inter = cc .* fit_dilate;
    if sum(inter(:)) > 0
        inters = inters | cc;
        cc_idx = [cc_idx, i];
    end
end
% [inters_row, inters_col] = find(inters);
% inters_num = length(inters_row);

er = zeros(1,cnt_num);
sr = zeros(1,cnt_num);
ar = zeros(1,cnt_num);
for i=1:cnt_num
    e = [];
    a_int = 1.414213/rr(i) / pi*180.0; % calculate interval of angle by radius, Sam, 220610
    a_num = round(360/a_int);
    a = zeros(1,a_num); % 0~360, Sam, 220609
    for j=1:length(cc_idx)
        % calculate error
        cc = ismember(lb, cc_idx(j)); % object with maximum area
        [cc_row, cc_col] = find(cc);
        d = sqrt( (cc_col-cx(i)).^2 + (cc_row-cy(i)).^2 );
        e = [e, (abs(d-rr(i)))'];
        
        % calculate angle, Sam, 220610
        t = atan2(cc_row-cy(i), cc_col-cx(i))/pi*180.0 + 180.0;
        t_round = min(round(t/a_int)+1, a_num);
        a(t_round) = 1;
    end

    % Sam, 220610
    er(i) = mean(e);
    sr(i) = std(e);
    ar(i) = sum(a(:))/a_num;
end
toc
end

% use Frechet distance to decide similarity of circle, Sam, 220607
function [cx, cy, rr, cnt, cm] = HoughCircle3(img, minr, maxr, stepr, stepa, percent)
debug = 1;

r = ceil((maxr-minr)/stepr); % number of radius need to be calculate
inters_t = ceil(360.0/stepa); % number of theta need to be calculate
[rows, cols] = size(img);
hspace = zeros(rows, cols, r);

[row, col] = find(img);
num = length(row);

% calculate hough transform
% a = x - r * cos(a)
% b = y - r * sin(a)
for i=1:num
    for j=1:r
        for k=1:inters_t
            aa = round(row(i) - (minr+(j-1)*stepr) * cosd(k*stepa));
            bb = round(col(i) - (minr+(j-1)*stepr) * sind(k*stepa));
            if aa>0 && aa<=rows && bb>0 && bb<=cols % must inside image
                hspace(aa,bb,j) = hspace(aa,bb,j)+1;
            end
        end
    end
end

% find maximum count in hough space
%par = max(max(max(hspace)));
par = max(hspace(:));
par_th = par*percent;
[m2, n2, r2] = size(hspace);

cx = []; cy = []; rr = [];
cnt = [];
% prob = []; % Sam, 220606
% num = m2 * n2 * r2;

for i=1:m2
    for j=1:n2
        for k=1:r2
            if hspace(i,j,k) >= par_th
                cx = [cx j]; cy = [cy i];
                rr = [rr minr+k*stepr];
                cnt = [cnt hspace(i,j,k)];
            end
        end
    end
end

% do connected component, Sam, 220606
lb = bwlabel(img);
sts = regionprops(logical(lb), 'Area', 'Image', 'BoundingBox');
if isempty(sts)
    return;
end
    
% find maximum count in hough space, Sam, 220606
[~, cnt_idx] = sort(cnt, 'descend');
% draw corresponding circle
fit = false(rows, cols);
fit_col = round( rr(cnt_idx(1)) * cos(0) + cx(cnt_idx(1)) );
fit_row = round( rr(cnt_idx(1)) * sin(0) + cy(cnt_idx(1)) );

% fit_sort = zeros(rows, cols);
% fit_sort(fit_row(1), fit_col(1)) = 1; % Sam, 220607
% figure(4), subplot(1,2,1), imshow(fit_sort);
for inters_t = 0:pi/50:2*pi
    x = round( rr(cnt_idx(1)) * cos(inters_t) + cx(cnt_idx(1)) );
    y = round( rr(cnt_idx(1)) * sin(inters_t) + cy(cnt_idx(1)) );
    if x<1 || x>cols || y<1 || y>rows
        continue;
    end
    fit(y,x) = 1;
    % store point
    if x~=fit_row(end) || y~=fit_col(end)
        fit_col = [fit_col x];
        fit_row = [fit_row y];
        
%         fit_sort(y, x) = 1; % Sam, 220607
%         figure(4), subplot(1,2,1), imshow(fit_sort);
    end
end
fit_num = length(fit_row); % Sam, 220607

nhood = [0 1 0;
         1 1 1;
         0 1 0];
fit_dilate = imdilate(fit, nhood);
 
% get all object with intersection
inters = false(rows, cols);
for i=1:length(sts)
    cc = ismember(lb, i); % object with maximum area
    inter = cc .* fit_dilate;
    if sum(inter(:)) > 0
        inters = inters | cc;
    end
end

% % sorting objects that intersect with fitting result, Sam 220607
% [row_i, col_i] = find(inters);
% num_i = length(row_i);
% inters_row = zeros(1,num_i);
% inters_col = zeros(1,num_i);
% 
% inters_row(1) = row_i(1);
% inters_col(1) = col_i(1);
% row_i = row_i(2:end);
% col_i = col_i(2:end);
% 
% % inters_sort = zeros(rows, cols);
% % inters_sort(inters_row(1), inters_col(1)) = 1; % Sam, 220607
% % figure(4), subplot(1,2,1), imshow(inters_sort);
% for i=2:num_i
%     d2_min = 999;
%     d2_min_idx = 0;
%     % search closest point
%     num_j = length(row_i);
%     for j=1:num_j
%         d2 = (inters_row(i-1)-row_i(j))^2 + (inters_col(i-1)-col_i(j))^2;
%         if d2<d2_min
%             d2_min = d2;
%             d2_min_idx = j;
%         end
%     end
%     % store closest point
%     inters_row(i) = row_i(d2_min_idx);
%     inters_col(i) = col_i(d2_min_idx);
%     % remove closest point
%     if d2_min_idx==1 && num_j>1
%         row_i = row_i(2:end);
%         col_i = col_i(2:end);
%     elseif d2_min_idx==num_j
%         row_i = row_i(1:end-1);
%         col_i = col_i(1:end-1);
%     else
%         row_i = [row_i(1:d2_min_idx-1); row_i(d2_min_idx+1:end);];
%         col_i = [col_i(1:d2_min_idx-1); col_i(d2_min_idx+1:end);];
%     end
%     
% %     inters_sort(inters_row(i), inters_col(i)) = 1; % Sam, 220607
% %     figure(4), subplot(1,2,1), imshow(inters_sort);
% end
% inters_num = length(inters_row);

% % search closest point of inters with first fit curve, Sam, 220607
% d2_min = 999;
% d2_min_idx = 0;
% for i=1:inters_num
%     d2 = (inters_row(i)-fit_row(1))^2 + (inters_col(i)-fit_col(1))^2;
%     if d2<d2_min
%         d2_min = d2;
%         d2_min_idx = i;
%     end
% end
% 
% % inters_1 = [inters_col(d2_min_idx), inters_row(d2_min_idx)];
% if d2_min_idx==1
%     inters_1f = [inters_col(2), inters_row(2)]; % forward
%     inters_1b = [inters_col(end), inters_row(end)]; % backward
% elseif d2_min_idx==inters_num
%     inters_1f = [inters_col(1), inters_row(1)]; % forward
%     inters_1b = [inters_col(end-1), inters_row(end-1)]; % backward
% else
%     inters_1f = [inters_col(d2_min_idx+1), inters_row(d2_min_idx+1)]; % forward
%     inters_1b = [inters_col(d2_min_idx-1), inters_row(d2_min_idx-1)]; % backward
% end
% 
% fit_1f = [fit_col(min(5,fit_num)) fit_row(min(5,fit_num))];
% d2_f = (fit_1f(1)-inters_1f(1))^2 + (fit_1f(2)-inters_1f(2))^2;
% d2_b = (fit_1f(1)-inters_1b(1))^2 + (fit_1f(2)-inters_1b(2))^2;
% if d2_f<d2_b % forward, Sam, 220607
%     if d2_min_idx==1
%         ;
%     elseif d2_min_idx==inters_num
%         inters_col = [inters_col(d2_min_idx), inters_col(1:d2_min_idx-1)];
%         inters_row = [inters_row(d2_min_idx), inters_row(1:d2_min_idx-1)];
%     else
%         inters_col = [inters_col(d2_min_idx:end), inters_col(1:d2_min_idx-1)];
%         inters_row = [inters_row(d2_min_idx:end), inters_row(1:d2_min_idx-1)];
%     end
% else
%     if d2_min_idx==1
%         inters_col = [fliplr(inters_col(2:end)), inters_col(1)];
%         inters_row = [fliplr(inters_row(2:end)), inters_row(1)];
%     elseif d2_min_idx==inters_num
%         inters_col = fliplr(inters_col);
%         inters_row = fliplr(inters_row);
%     else
%         inters_col = [fliplr(inters_col(1:d2_min_idx)), fliplr(inters_col(d2_min_idx+1:end))];
%         inters_row = [fliplr(inters_row(1:d2_min_idx)), fliplr(inters_row(d2_min_idx+1:end))];
%     end
% end

% sort points of intersection by its angle between center of circle, Sam, 220607
[row_i, col_i] = find(inters);
inters_t = atan2(row_i - cy(cnt_idx(1)), col_i - cx(cnt_idx(1)));
[inters_t_val, inters_t_idx] = sort(inters_t);
inters_col = col_i(inters_t_idx)';
inters_row = row_i(inters_t_idx)';
inters_num = length(inters_row);

fit_t = atan2(fit_row - cy(cnt_idx(1)), fit_col - cx(cnt_idx(1)));
[~, fit_t_idx] = min(abs(fit_t-inters_t_val(1)));
if fit_t_idx==1
    ;
elseif fit_t_idx==fit_num
    fit_col = [fit_col(end), fit_col(1:end-1)];
    fit_row = [fit_row(end), fit_row(1:end-1)];
else
    fit_col = [fit_col(fit_t_idx:end), fit_col(1:fit_t_idx-1)];
    fit_row = [fit_row(fit_t_idx:end), fit_row(1:fit_t_idx-1)];
end

if debug % Sam, 220608
    % debug
    fit_sort = zeros(rows, cols);
    for i=1:fit_num
        fit_sort(fit_row(i), fit_col(i)) = 1; % Sam, 220607
        figure(4), subplot(1,2,1), imshow(fit_sort);
    end

    % debug
    inters_sort = zeros(rows, cols);
    for i=1:inters_num
        inters_sort(inters_row(i), inters_col(i)) = 1; % Sam, 220607
        figure(4), subplot(1,2,2), imshow(inters_sort);
    end
end

% combine as two dimaintion array, Sam, 220607
fit_sort = [fit_col' fit_row'];
inters_sort = [inters_col' inters_row'];

% calculate discrete Frechet distance, Sam, 220607
% t = 0:pi/8:2*pi;
% y = linspace(1,3,6);
% P = [(2:7)' y']+0.3.*randn(6,2);
% Q = [t' sin(t')]+2+0.3.*randn(length(t),2);
P = fit_sort;
Q = inters_sort;
[cm, cSq] = DiscreteFrechetDist(P, Q);

if debug % Sam, 220608
    % plot result
    clf(figure(5));
    plot(Q(:,1),Q(:,2),'o-r','linewidth',1,'markerfacecolor','r')
    hold on
    plot(P(:,1),P(:,2),'o-b','linewidth',1,'markerfacecolor','b')
    title(['Discrete Frechet Distance of curves P and Q: ' num2str(cm)])
    legend('Q','P','location','best')
    line([2 cm+2],[0.5 0.5],'color','m','linewidth',1)
    text(2,0.4,'dFD length')
    for i=1:length(cSq)
      line([P(cSq(i,1),1) Q(cSq(i,2),1)],...
           [P(cSq(i,1),2) Q(cSq(i,2),2)],...
           'color',[0 0 0]+(i/length(cSq)/1.35));
    end
    axis equal
    % display the coupling sequence along with eac7h distance between points
    disp([cSq sqrt(sum((P(cSq(:,1),:) - Q(cSq(:,2),:)).^2,2))])
end

% Sam, 220608
cx = cx(cnt_idx(1));
cy = cy(cnt_idx(1));
rr = rr(cnt_idx(1));
cnt = cnt(cnt_idx(1));

end

% Sam, 220602
function [cx, cy, aa, bb, theta] = HoughEllipse(img)
min_axis = 6;
max_axis = 18;
accuracy = 10;
hist_th = 0;
[rows, cols] = size(img);
[row, col] = find(img);
num = length(row); 

accumulator = [];
bin_size = accuracy*accuracy;
max_axis2 = min_axis^2;

cx = [];
cy = [];
aa = [];
bb = [];
theta = [];

for p1=1:num
    for p2=1:p1
        x1 = col(p1); y1 = row(p1);
        x2 = col(p2); y2 = row(p2);
        dx = x1-x2;   dy = y1-y2;
        a = 0.5 * sqrt(dx^2+dy^2);
        if a>=min_axis && a<=max_axis
            x0 = 0.5 * (x1+x2);
            y0 = 0.5 * (y1+y2);
            for p=1:num
                dx = col(p)-x0;
                dy = row(p)-y0;
                d = sqrt(dx^2+dy^2);
                if d>min_axis
                    dx = col(p)-x2;
                    dy = row(p)-y2;
                    f2 = dx^2+dy^2;
                    c = (a^2+d^2-f2) / (2*a*d);
                    c2 = c^2;
                    k = a^2-d^2*c2;
                    if k>0 && c2<1
                        b2 = a^2*d^2*(1-c2) / k;
                        if b2 > max_axis2
                            accumulator = [accumulator b2];
                        end
                    end
                end
            end
            
            if ~isempty(accumulator)
                bins = [];
                acc_max = max(accumulator);
                for d=0:bin_size:(acc_max+bin_size)
                    bins = [bins d];
                end
                
                hist = zeros(1,length(bins)-1);
                for p=1:length(accumulator)
                    a = accumulator(p);
                    for j=1:length(bins)-1
                        left = bins(j);
                        right = bins(j+1);
                        if a>=left && a<right
                            hist(j) = hist(j)+1;
                            break;
                        end
                    end
                    if a==bins(end)
                        hist(end-1) = hist(end-1)+1;
                    end
                end
                
                [hist_max, hist_max_idx] = max(hist);
                if hist_max > hist_th
                     cx = [cx x0]; cy = [cy y0];
                     aa = [aa a]; bb = [bb sqrt(bins(hist_max_idx))];
                     t = atan2(y2-y1, x2-x1);
                     theta = [theta t];
                end
                accumulator = [];
            end
        end
    end
end

end


% Sam, 220602.
function [cx, cy, a, b, theta] = HoughEllipse2(img, axis_min, axis_max, stepa, percent)
t = ceil(180.0/stepa); % number of theta need to be calculate
[rows, cols] = size(img);
[row, col] = find(img);
num = length(row);

dists_max = zeros(rows, cols);
dists = zeros(1,num);

% calculate distance between points of image and edge
for i=1:rows
    for j=1:cols
        for k=1:num
            dx = j-col(k);
            dy = i-row(k);
            dists(k) = sqrt(dx^2 + dy^2);
        end
        dist_max = max(dists);
        dists_max(i,j) = dist_max;
    end
end

% calculate long axis
dist_min = min(dists_max(:));
[cy, cx] = find(dists_max==dist_min);

% angle = zeros(1,num);
% for i=1:num
%     angle(i) = atan2(cy-row(i), cx-col(i));
% end
% angle_sort = sort(angle(:));
% figure(5);
% hold on;
% for i=1:num
%     plot(i, angle_sort(i)/pi*180.0, 'ro');
% end
% hold off;



aa = ceil(dist_min);
hspace = zeros(aa, t);
for i=1:num
    for j=1:t
        p = cx(1);
        q = cy(1);
        x = col(i);
        y = row(i);
        c = cosd(j*stepa);
        s = sind(j*stepa);
        p1 = ( (x-p)*c + (y-q)*s)^2 / aa^2;
        p2 = (-(x-p)*s + (y-q)*c)^2;
        bb = ceil( sqrt(p2/(1-p1)) );
        if bb>0 && bb<=aa
            hspace(bb,j) = hspace(bb,j)+1;
        end
    end
end

% find maximum count in hough space
par = max(max(hspace));
par_th = percent*par;
[m2, n2] = size(hspace);

a = aa;
b = [];
theta = [];

for i=1:m2
    for j=1:n2
        if hspace(i,j) >= par_th
            b = [b i];
            theta = [theta j*stepa];
        end
    end
end

end

function [cx, cy] = RandHoughEllipse(img)
[rows, cols] = size(img);
[row, col] = find(img);
num = length(row);

neighbor = 5;
iter = floor(0.75*num);
cx = zeros(1,iter);
cy = zeros(1,iter);

figure(5);
imshow(img);
hold on;

for i=1:iter
    rnd = randi([1 num], 1, 3);

    pts_x = zeros(1,3); 
    pts_y = zeros(1,3);
    bs = [];
    ks = [];
    for j=1:3
        pts_x(j) = col(rnd(j));
        pts_y(j) = row(rnd(j));
        plot(pts_x, pts_y, 'b.');
            
        xs = max(1, pts_x(j) - (neighbor-1)/2);
        xe = min(xs + neighbor, cols);
        ys = max(1, pts_y(j) - (neighbor-1)/2);
        ye = min(ys + neighbor, rows);
        
        [y, x] = find(img(ys:ye, xs:xe)>0);
        x = x + xs;
        y = y + ys;
%         plot(x, y, 'b.');
%         [k, b] = fitline(x, y);
        p = polyfit(x,y,1);
        k = p(1);
        b = p(2);
%         xx = 1:cols;
%         yy = xx*k+b;
%         plot(xx,yy,'g');
        
        
        ks = [ks k];
        bs = [bs pts_y(j) - k*pts_x(j)];
    end
    
    % find center
    line01 = [1 2];
    line12 = [2 3];
    lines = [line01; line12];
    intss_x = zeros(1,length(lines));
    intss_y = zeros(1,length(lines));
    slopes = zeros(1,length(lines));
    intcs =  zeros(1,length(lines));
    for j=1:length(lines)
        line = lines(j,:);
        k1 = ks(line(1)); k2 = ks(line(2));
        b1 = bs(line(1)); b2 = bs(line(2));
        
        ints_x = (b2-b1) / (k1-k2);
        ints_y = k1*ints_x + b1;
        intss_x(j) = ints_x;
        intss_y(j) = ints_y;
        
        mid_x = (pts_x(line(1)) + pts_x(line(2))) / 2.0;
        mid_y = (pts_y(line(1)) + pts_y(line(2))) / 2.0;
        
        slopes(j) = (mid_y - ints_y) / (mid_x - ints_x);
        intcs(j) = ints_y - slopes(j) * ints_x;
    end
    
    k1 = slopes(1); k2 = slopes(2);
    b1 = intcs(1); b2 = intcs(2);
    cx(i) = (b2-b1) / (k1-k2);
    cy(i) = k1 * cx(i) + b1;
    plot(cx(i), cy(i), 'rx');
end

x = 0;
y = 0;
sum=0;
for i=1:iter
    if isnan(cx(i)) || isnan(cy(i))
        continue;
    end
    y = y + cy(i);
    x = x + cx(i);
    sum = sum + 1;
end
xx = x/sum;
yy = y/sum;
plot(xx, yy, 'g+');

hold off;

end

function [k,b] = fitline(x,y)
n=length(x);
x=reshape(x,n,1);%?????
y=reshape(y,n,1);
A=[x,ones(n,1)];
bb=y;
B=A'*A;
bb=A'*bb;
yy=B\bb;
k=yy(1);
b=yy(2);
end

%% This is the implement of Hilditch thinning algorithm
% Author: Junjun Li
% Email: junlee_happy@163.com
% Source: Linear Skeletons from Square Cupboards (C.Judith Hilditch)
function PIC = Thinning(PIC)
% %% ????????
% PIC = imread('shou.bmp');
% 
% if islogical(PIC) == 0
%     PIC = im2bw(PIC);
% end

%% ???????PIC?????????????????
[imageRow, imageCol] = size(PIC);
midPIC = zeros([imageRow+4, imageCol+4]);
midPIC(3:end-2, 3:end-2) = PIC;
PIC = midPIC;

PIC2 = PIC;                   % pic2??MATLAB????????????????  

%% Hilditch??
C = 1;    % ????

while C
    C = 0;  
    [row, col] = find(PIC == 1);    % ??????1?????
    PIC_DEL = ones(size(PIC));
    count = 0; % ????????????????0???????
    
    for i = 1:length(row)
        
        % ????1????
        X = row(i);
        Y = col(i);
        
        % 3*3????????????????P??
        % X4 X3 X2
        % X5 P  X1
        % X6 X7 X8
        P = [PIC(X, Y+1), PIC(X-1, Y+1), PIC(X-1, Y), ...
             PIC(X-1, Y-1), PIC(X, Y-1), PIC(X+1, Y-1), ...
             PIC(X+1, Y), PIC(X+1, Y+1), PIC(X, Y+1)];
        
        % condition2?X1?X3?X5?X7???1
        C2 = P(1) + P(3) + P(5) + (7);
        if C2 == 4
            continue;
        end
        
        % condition3???????X1~X8??????1
        if sum(P(1:8)) < 2
            continue;
        end
        
        % condition4?????????????????????????1???remove
        w = 0;
        for m = X-1:X+1
            for n = Y-1:Y+1
                if m == X && n == Y
                    continue;
                end
                if PIC(m, n) == 1 && PIC_DEL(m, n) == 1
                    w = w + 1;
                end
            end
        end
        if w < 1
            continue;
        end
        
        % condition5???????8?????????
        
        if c_number(P) ~= 1
            continue;
        end
        
        % condition6??????????????????????????X3?X5
        if PIC_DEL(X-1, Y) == 0
            P3 = [PIC(X, Y+1), PIC(X-1, Y+1), 0, ...
                  PIC(X-1, Y-1), PIC(X, Y-1), PIC(X+1, Y-1), ...
                  PIC(X+1, Y), PIC(X+1, Y+1), PIC(X, Y+1)];
             if c_number(P3) ~= 1
                 continue;
             end
        end
        
        % condition6?X5
        if PIC_DEL(X, Y-1) == 0
            P5 = [PIC(X, Y+1), PIC(X-1, Y+1), PIC(X-1, Y), ...
                  PIC(X-1, Y-1), 0, PIC(X+1, Y-1), ...
                  PIC(X+1, Y), PIC(X+1, Y+1), PIC(X, Y+1)];
             if c_number(P5) ~= 1
                 continue;
             end
        end
        
        PIC_DEL(X, Y) = 0;
        count = count + 1;
    end
    if count ~= 0
        PIC = PIC .* PIC_DEL;    % ??PIC?????????PIC???????
        C = 1;
    end
end

end

function result=c_number(P)
N = 0;
for k = 1:4
    temp = ~P(2*k-1) - ~P(2*k-1)*~P(2*k)*~P(2*k+1);
    N = N + temp;
end
result = N;
end

% Sam, 220606
% by : Peijin Zhang
% Calculates the curvature using three points, written in Python and MATLAB
% https://github.com/peijin94/PJCurvature
function [kappa,norm_k] = PJcurvature(x,y)
    x = reshape(x,3,1);
    y = reshape(y,3,1);
    t_a = norm([x(2)-x(1),y(2)-y(1)]);
    t_b = norm([x(3)-x(2),y(3)-y(2)]);
    
    M =[[1, -t_a, t_a^2];
        [1, 0,    0    ];
        [1,  t_b, t_b^2]];

    a = M\x
    b = M\y

    kappa  = 2.*(a(3)*b(2)-b(3)*a(2)) / (a(2)^2.+b(2)^2.)^(1.5);
    norm_k =  [b(2),-a(2)]/sqrt(a(2)^2.+b(2)^2.);
end