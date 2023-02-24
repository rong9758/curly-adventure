% lensShading calculates the lens shading map
%
% INPUT:
%   IDraw: raw image data
%   bayerFormat: bayer format:
%       'bggr'
%       'rggb'
%       'grbg'
%       'gbrg'
%   rois: An array containing ROIs defined by ERS
%        Format:[wBayer hBayer] for each ROI
% OUTPUT:
%
%  Author: Eugene Lam
%          November 2012
%******************************************************************
% ver 6.2-1 2016-07-08: version number initialize
%******************************************************************

function output  = lensShading_monochrome_noload(IDraw, bayerFormat, pedestal, bitDepth, rois)
    %% change for monochrome: will measure lens shading rather than color
    %% shading 
% parameters for getFile
%bayerFormat = 'rggb'; % Sam, 220413
%bitDepth = 10; % Sam, 220413
%pedestal = 0; % have substracted 64 in getfile % Sam, 220413
debug =0;
roisW = [0.005 0.12875 0.2525 0.37625 0.5   0.62375 0.7475 0.87125 0.995];
roisH = [0.005 0.12875 0.2525 0.37625 0.5   0.62375 0.7475 0.87125 0.995];

h = length(roisH);
w = length(roisW);
for j = 1:h
    for i = 1:w
        rois(i + (j-1)*w,1) = roisW(i);
        rois(i + (j-1)*w,2) = roisH(j);  
        
        if (i == 1 && j == 1) 
            rois(i + (j-1)*w,2) = 0.0248;
            rois(i + (j-1)*w,1) = 0.0248;
        end
        if (i == 1 && j == h)
            rois(i + (j-1)*w,2) = 0.9752;
            rois(i + (j-1)*w,1) = 0.0248;
        end
        if (i == w && j == 1)
            rois(i + (j-1)*w,2) = 0.0248;
            rois(i + (j-1)*w,1) = 0.9752;
        end
        if (i == w && j == h)
            rois(i + (j-1)*w,2) = 0.9752;
            rois(i + (j-1)*w,1) = 0.9752;
        end
    end
end

%     IDbayer = preprocess(IDraw, bayerFormat, 'bayer', bitDepth, pedestal, 0, false, true);
    IDraw = preprocess(IDraw, bayerFormat, 'raw', bitDepth, pedestal, 83.5, false, true); % must input bitDepth and pedestal, Sam, 220413
    if debug
        [h w c]=size(IDraw);
        figure;imshow(IDraw,[])
        title('colourShading on the juliet flat-field image; blue label is color median values')
        rectangle('Position', [1,1,73,60], 'EdgeColor', 'y');
        rectangle('Position', [w-72,1,73,60], 'EdgeColor', 'y');
        rectangle('Position', [1,h-59,73,60], 'EdgeColor', 'y');
        rectangle('Position', [w-72,h-59,73,60], 'EdgeColor', 'y');
        text(1, 1+30,['73x62 pix'],'Color','g', 'FontSize', 12);
        text(w-72,1+30,['73x62 pix'],'Color','g', 'FontSize', 12);
        text(1,h-59+30,['73x62 pix'],'Color','g', 'FontSize', 12);
        text(w-72,h-59+30,['73x62 pix'],'Color','g', 'FontSize', 12);
    end
    % end: change for monochrome
    %%
    [h w c] = size(IDraw);
    roiSize = 10;
    
    % calculate the locations for the ROI
    %   [x y xLeft yTop xRight yBottom]
    coordROI = zeros(length(rois),6);
    for i = 1:length(rois)
        roiX = round(w * rois(i,1)); % centre X value
        roiY = round(h * rois(i,2)); % centre Y value
        coordROI(i,:) = [roiX ...
                         roiY ...
                         int16(roiX - roiSize/2 + 1.4999) ...
                         int16(roiY - roiSize/2 + 1.4999) ...
                         int16(roiX + roiSize/2 - 0.5) ...
                         int16(roiY + roiSize/2 - 0.5)];
    end
    % need to provide additional clipping in case ROI size exceeds bounds
    coordROI(:,[3 4]) = max(coordROI(:,[3 4]), 1);
    coordROI(:,5) = min(coordROI(:,5), w);
    coordROI(:,6) = min(coordROI(:,6), h);
    
    % calculate the medians for all these locations
    elem_x = length(unique(rois(:,1))) - 2;     % number of grid values on x axis
    elem_y = length(rois) / elem_x;         % number of grid values on y axis
    medianROI = zeros(elem_y, elem_x, c);
    for i = 1:elem_y
        for j = 1:elem_x
            % create the ROI
            indexLocation = j + (i-1) * elem_x;
            ROI = IDraw(coordROI(indexLocation,4):coordROI(indexLocation,6),...
                        coordROI(indexLocation,3):coordROI(indexLocation,5),:);
            % reshape for median calculation
            %% change for monochrome 
            ROIreshape = reshape(ROI, [size(ROI,1)*size(ROI,2) 1]);
            medianROI(i,j) = median(double(ROIreshape(:,:)));
            % can produce decimals since we have an even number of elements            
            % end change for monochrome
            
            %% debug plot 
            if debug
            rectangle('Position', [ coordROI(indexLocation,1), coordROI(indexLocation,2),abs( coordROI(indexLocation,3)-coordROI(indexLocation,1))+1,abs(coordROI(indexLocation,4)-coordROI(indexLocation,2))+1], 'EdgeColor', 'k');
            text(coordROI(indexLocation,1),coordROI(indexLocation,2) ,num2str(median(double(ROIreshape(:,:)))),'Color','b', 'FontSize', 12);
            end
        end
    end    
    scaledROI = ones(size(medianROI));
    %% change for monochrome
%     for i = 1:c
%         scaledROI(:,:,i) = uint8(250 * double(medianROI(:,:,i)) / double(max(max(medianROI(:,:,i)))));
%     end
    scaledROI = uint8(250 * double(medianROI) / double(max(max(medianROI))));
    % end: change for monochrome
    %%
    % additional code for confirming NVM values
    %% change for monochrome
    scaledROI = uint8(scaledROI);
    scaledROI_reshape = scaledROI;
    
%     for i = 1:size(scaledROI, 1)
%         %% change for monochrome
%         scaledROI_reshape = [scaledROI_reshape; scaledROI(i,:)];
%         %%
%         % end change for monochrome
%     end
    % end: change for monochrome
    %%

    scaledROI_reshape = scaledROI_reshape';
    scaledROI_hex = dec2hex(scaledROI_reshape(:));
    
    new_ROI = double(scaledROI_reshape);
    ls_mat = [rois(:,1), rois(:,2),reshape(new_ROI, [size(scaledROI,1)*size(scaledROI,2) 1])];
   
    %output = scaledROI;
    output = ls_mat;

end
