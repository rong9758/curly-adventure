% preprocess takes raw image data, applies any preprocessing, and reformats
% to desired output
%
% INPUT:
%   ID: image data
%   bayerFormat: bayer format:
%       'bggr'
%       'rggb'
%       'grbg'
%       'gbrg'
%   outputFormat: requested output
%       'raw'
%       'bayer'
%       'rgb'
%       'yuv'
%   pedestal: amount of pedestal to add
%   FOV: nominal design FOV for the program, set to 0 for no lens shading
%   AWB: determine whether to apply software AWB or not
%   signed: boolean for signed or not
%   option_lsc: method for shading correction
%   ls_mat: shading metrics in the format of
%       [roiX1, roiY1, roiMean1
%        roiX2, roiY2, roiMean2
%        ...     ...     ...]
%   3rd dimension is quad space i.e. TL, TR, BL, BR
%
% OUTPUT:
%   ID: image array as requested
%
%  Author: Aaron Cong
%          March 2016
%******************************************************************
% ver 6.1-1 Lens Shading Correction Upate
% ver 6.2-1 2016-06-06 : per module Lens Shading Correction Update
% ver 5.0-1 2016-09-07 : add embedded data line cropping
% ver 5.1-1 2016-12-01 : update lsc curve fitting from 6th to 12th
% ver 5.1-2 2016-12-06 : Correct fitting coefficient expansion
% ver 5.1-3 2017-01-06 : add condition to cropping: 1314
%******************************************************************

function [ID] = preprocess(raw, bayerFormat, outputFormat, outputBitDepth, pedestal, FOV, AWB, signed, option_lsc, ls_mat) 
% parameters
debug = false;
ID = double(raw);

% Cropping embedded data line, Before cropping 1314 x 1104, After cropping,
% 1312 x 1104
if size(ID,1)==1314
    ID=ID(3:end,:);
end

% pedestal formatting
ID = ID + pedestal;
if ~signed
    ID(ID<0) = 0; % don't clip
end

% white balance (colour plane agnostic)
if AWB
    ID = whiteBalance(ID);
    ID(ID>2^outputBitDepth-1) = 2^outputBitDepth-1; % don't oversaturate
end

% convert to desired output format
switch lower(outputFormat)
    case {'raw'}
        ID = ID;
    case {'bayer'}
        ch1 = ID(1:2:end,1:2:end);
        ch2 = ID(1:2:end,2:2:end);
        ch3 = ID(2:2:end,1:2:end);
        ch4 = ID(2:2:end,2:2:end);
        switch (bayerFormat)
            case 'rggb'
                R=ch1; Gr=ch2; Gb=ch3; B=ch4;
            case 'bggr'
                B=ch1; Gb=ch2; Gr=ch3; R=ch4;
            case 'gbrg'
                Gb=ch1; B=ch2; R=ch3; Gr=ch4;
            case 'grbg'
                Gr=ch1; R=ch2; B=ch3; Gb=ch4;
            case 'monochrome'
                R=ch1; Gr=ch2; Gb=ch3; B=ch4;
        end
        ID = [];
        ID(:,:,1) = R;
        ID(:,:,2) = Gr;
        ID(:,:,3) = Gb;
        ID(:,:,4) = B;
end

% lens shading
if (FOV~=0)
%     fprintf('nargin: %d', nargin);
    if nargin> 9
        switch option_lsc
            case (0) % original dummy shading correction
                ID = lensShading(ID,FOV);
            case (1) % formula based shading correction
                ID = lensShading_plus(ID);
            case (2) % bilinear shading correction
                ID = lensShading_bilinear(ID, bayerFormat, pedestal, outputBitDepth, ls_mat(:,1:2));
            case (3) % flat field per module shading correction
                ID = lensShading_nvm(ID, ls_mat);
        end
    else
        disp('Incomplete lens shading input!')
    end
    ID(ID>2^outputBitDepth-1) = 2^outputBitDepth-1; % don't oversaturate
end

ID = round(ID);
end

%%
% lensShading applies lens shading and returns the image back
%
% INPUT:
%   ID_preLSC: raw image input (this part is critical)
%   ls_mat: lens shading matrix in the formation of [x           ]
%
% OUTPUT:
%   ID_postLSC: corrected values with lens shading applied
%
%  Author: Aaron Cong
%          May 2016
%******************************************************************
function ID_postLSC=lensShading_nvm(ID_preLSC, ls_mat)
createREF=0;
alldots=0;% overwritting curve fitting with all dots fitting (temporary)
debug=0;
[h w c] = size(ID_preLSC);
colorSet= {'r','g','c','b'};
for cIdx =1:c
    % find centre region and calculate distances
    centreX = w/2 + 0.5;
    centreY = h/2 + 0.5;
    halfDiag = sqrt((centreX-1)^2 + (centreY-1)^2);
    % overwrite with calculated shading oc
    %     centreX=ls_oc(1);
    %     centreY=ls_oc(2);
    % calculate the locations for the ROI
    %   [x y, r]
    rois=ls_mat(:,1:2,cIdx);
    coordROI = zeros(length(rois),4);
    for i = 1:length(rois)
        roiX = round(w * rois(i,1)); % roi centre X value
        roiY = round(h * rois(i,2)); % roi centre Y value
        roiR = sqrt((centreX-roiX)^2 + (centreY-roiY)^2); % roi Radius to optical center (X, Y)
        roiR_norm = roiR/ halfDiag; % normalized distance by half diagnol field of view (0.0 ~ 1.0)
        coordROI(i,:) = [roiX, roiY, roiR, roiR_norm];
    end
    % non-normalized shading values (0 ~ 256) & roi distance to center
    ls_data=[coordROI(:,4) ls_mat(:,3,cIdx)];
    
    % symmetric fitting curve around center
    ls_data_mirror=[-coordROI(:,4) ls_mat(:,3,cIdx)];
    ls_data=[ls_data;ls_data_mirror];
    
    % 12th order polynomial fitting
    order_fit=12; % use 12th order fitting
    [ls_fit,S]=polyfit(ls_data(:,1),ls_data(:,2),order_fit); % (READ POINT: ls_fit, S)
    coeff=ls_fit;

    %% calculate shading correction using all points
    if alldots
    maxID=mean(mean(ID_preLSC(h/2-5:h/2+5,w/2-5:w/2+5)));
    for i=1:w
        for j=1:h
            Fpixel(j,i) = sqrt((centreX-i)^2 + (centreY-j)^2) / halfDiag;
            ls_data_all((i-1)*h+j,1)=Fpixel(j,i);
            ls_data_all((i-1)*h+j,2)=ID_preLSC(j,i)/maxID*250;
            i
        end
    end
    
    % symmetric fitting curve around center
    ls_data_mirror=[-ls_data_all(:,1) ls_data_all(:,2)];
    ls_data=[ls_data_all;ls_data_mirror];
    
    % 12th order polynomial fitting
    order_fit=12; % use 12th order fitting
    [ls_fit,S]=polyfit(ls_data(:,1),ls_data(:,2),order_fit); % (READ POINT: ls_fit, S)
    coeff=ls_fit;
    end
    
    %% debut plot
    if debug
        % debug plot
        figureRaw=figure('position', [0, 0, 1500, 1500]) ;
        x1=linspace(0,0.98);
        y1=polyval(ls_fit,x1);
        plot(ls_data(:,1),ls_data(:,2),'o','Color',colorSet{cIdx});hold on; plot(x1,y1,'--','Color',colorSet{cIdx});grid on
        ls_data_4channel(:,1)=ls_data(:,1);
        ls_data_4channel(:,cIdx+1)=ls_data(:,2);
        ls_fit_4channel(:,cIdx)=ls_fit';
        csvwrite('ls_data_all.csv',ls_data_4channel);
        csvwrite('ls_fit_all.csv',ls_fit_4channel);
        legend('lsData','lsFit','lensModel');
        title('Shading curve fitting comparison')
    end
    %% calculate correction coefficient
    for i=1:w
        for j=1:h
            Fpixel(j,i) = sqrt((centreX-i)^2 + (centreY-j)^2) / halfDiag;
            lsREF(j,i)=0;
            % fitting model coeffcient
            for idx_coeff=1:length(coeff)
                % create reference light field images
                lsREF(j,i) = lsREF(j,i) + coeff(idx_coeff)*Fpixel(j,i)^(order_fit-idx_coeff+1);
            end
            % create gain correction per pixels
            lscApply_nvm(j,i)= 1/lsREF(j,i);
        end
    end
    % normalized by fitting peak
    fit_peak=min(min(lscApply_nvm));
    % apply correction coefficient
    ID_postLSC(:,:,cIdx) = lscApply_nvm.* double(ID_preLSC(:,:,cIdx))/fit_peak;% (READ POINT: ID_postLSC)
    if createREF
        % create referece images
        ID_refLSC(:,:,cIdx) = lsREF*fit_peak*max(max(ID_preLSC(:,:,cIdx))); % (READ POINT: ID_refLSC)
    end
end
end

% lensShading applies lens shading and returns the image back
%
% INPUT:
%   ID: bayer image input (this part is critical)
%   relativeIllumination: estimated fallout towards the edges
%
% OUTPUT:
%   ID: corrected values with lens shading applied
%
%  Author: Aaron Cong
%          March 2016
%******************************************************************

function ID_postLSC = lensShading_plus(ID_preLSC)

[h w c] = size(ID_preLSC);

% find centre region and calculate distances
centreX = w/2 + 0.5;
centreY = h/2 + 0.5;
halfDiag = sqrt((centreX-1)^2 + (centreY-1)^2);

% calculate correction coefficient
for i=1:w
    for j=1:h
        Fpixel(j,i) = sqrt((centreX-i)^2 + (centreY-j)^2) / halfDiag;
        
        % define model coeffcient
        % model y = -44.566x4 + 209.89x3 - 230.61x2 - 1.9002x + 99.521
        c4= -44.566;
        c3= 209.89;
        c2= - 230.61;
        c1= - 1.9002;
        c0= 99.521;
        
        lscApply_plus(j,i) = 100/(c4*Fpixel(j,i)^4+c3* Fpixel(j,i)^3+c2* Fpixel(j,i)^2+c1*Fpixel(j,i)^1+c0);
    end
end
% apply correction coefficient
ID_postLSC = lscApply_plus.* double(ID_preLSC);

end

%%
% lensShading applies lens shading and returns the image back
%
% INPUT:
%   ID: bayer image input (this part is critical)
%   relativeIllumination: estimated fallout towards the edges
%
% OUTPUT:
%   ID: corrected values with lens shading applied
%
%  Author: Eugene Lam
%          September 2012
%******************************************************************

function ID_postLSC = lensShading(ID_preLSC, FOV)

[h w c] = size(ID_preLSC);

% find centre region and calculate distances
centreX = w/2 + 0.5;
centreY = h/2 + 0.5;
halfDiag = sqrt((centreX-1)^2 + (centreY-1)^2);

% calculate correction coefficient
for i=1:w
    for j=1:h
        FOVpixel(j,i) = 180/pi * atan(tan(pi/180*(FOV/2)) * sqrt((centreX-i)^2 + (centreY-j)^2) / halfDiag);
        lscApply(j,i) = 1/(cos(pi/180 * FOVpixel(j,i)))^4;
    end
end

% apply correction coefficient
ID_postLSC = lscApply .* double(ID_preLSC);

end
%%
% lensShading applies lens shading and returns the image back
%
% INPUT:
%   ID: bayer image input (this part is critical)
%   relativeIllumination: estimated fallout towards the edges
%
% OUTPUT:
%   ID: corrected values with lens shading applied
%
%  Author: Aaron Cong
%          March 2016
%******************************************************************

function ID_postLSC = lensShading_bilinear(ID_preLSC, bayerFormat, pedestal, bitDepth, roishading)

[h w c] = size(ID_preLSC);

% calculate monochrome lens shading map
roishading_values = lensShading_monochrome(ID_preLSC, bayerFormat, pedestal, bitDepth, roishading);
roishading_values = double(roishading_values)/250;

%!!! pay attention, if shifting corner, the bilinear interpolaiton need to update.pad to obtain border map
roishading_valuespd =  padarray(roishading_values,[1 1],'replicate','both');
roishading_x = padarray(roishading(1:size(roishading_values,2),1),1,1,'both')';
roishading_x(1) = 0;
roishading_y = padarray(unique(roishading(:,2)),1,1,'both')';
roishading_y(1) = 0;

% 2D bilinear interpolate color shading map to obtain per pixel shading
perPixelLS = bilinearinterp2D(roishading_x,roishading_y,roishading_valuespd,1:w,1:h);

% lens shading correction
ID_postLSC = ID_preLSC./perPixelLS;
ID_postLSC = min(ID_postLSC,2^bitDepth-1);

%%
% bilinearinterp2D performs 2D bilinear interpolation using grd values at points specified by
% matrices x and y to find grd_interped at the points in matrices xo and yo.
%
% INPUT:
%   x, y: sample locations of input grd, 1D matrix, 0<=x<=1, 0<=y<=1
%   grd: image with values at sample locations specified by x and y
%   xo, yo: cooridnate locations of interpolated output grd_interped,
%   1D matrix, x=1:wo, y=1:ho, where [ho wo] = size(grd_interped)
%
% OUTPUT:
%   grd_interped: 2D bilinear interpolated image
%
% COMMENTS:
%   Note that the algorithm uses simple bilinear interpolation. Some special considerations:
%       - Pixel values for a colour which is the origin colour have no
%       interpolation
%       - Borders use a truncated interpolation, averaging only the
%       directly adjacent pixels of that colour
%       - No intepolation is done outside of sampled input points, 0 will
%       be returned, e.g. grd_interped(1,1) = 0 if x(1)>1, y(1)>1
%
%  Author: Maria Zhao
%          April 2013
%******************************************************************

    function grd_interped = bilinearinterp2D(x,y,grd,xo,yo)
        
        % set up after interpolation output size
        wo = numel(xo); ho = numel(yo);
        
        % set up input coordinates and size
        [hi wi ci] = size(grd);
        xi = round(x*wo); yi = round(y*ho);
        if xi(1)==0, xi(1)=1; end
        if yi(1)==0, yi(1)=1; end
        
        % no interpolation done outside of sampled inputs
        xoinit = xi(1); xoend = xi(end);
        yoinit = yi(1); yoend = yi(end);
        
        % initialize output
        grd_interped = zeros(ho,wo,ci);
        
        for i=1:ci
            % input cooridinates counter
            xiCounter = 1;yiCounter = 1;
            
            % first calculate column interpolation
            colInterped = zeros(hi,wo);
            for j=xoinit:xoend
                xp = xo(j); xin = xi(xiCounter);
                if (xp==xin)   % assign input grid column to output column if output x coordinate same as input x coordinate
                    colInterped(:,j) = grd(:,xiCounter,i);
                    xiCounter = xiCounter+1;
                else           % otherwise bilinear interpolate using input grid column values to get output column
                    xinPrev = xi(xiCounter-1);
                    colInterped(:,j) = (xin-xp)/(xin-xinPrev)*grd(:,xiCounter-1,i) + (xp-xinPrev)/(xin-xinPrev)*grd(:,xiCounter,i);
                end
            end
            
            % row interpolation on column interpolated output to get final output
            for k=yoinit:yoend
                yp = yo(k); yin = yi(yiCounter);
                if (yp==yin)   % assign input grid row to output row if output y coordinate same as input y coordinate
                    grd_interped(k,:,i) = colInterped(yiCounter,:);
                    yiCounter = yiCounter+1;
                else           % otherwise bilinear interpolate using input row values to get output row
                    yinPrev = yi(yiCounter-1);
                    grd_interped(k,:,i) = (yin-yp)/(yin-yinPrev)*colInterped(yiCounter-1,:) + (yp-yinPrev)/(yin-yinPrev)*colInterped(yiCounter,:);
                end
            end
        end
        
    end


end

%%
% whiteBalance applies whiteBalance and returns the image back
%
% INPUT:
%   ID: bayer image input (this part is critical)
%
% OUTPUT:
%   ID: corrected values with white balance applied
%
%  Author: Eugene Lam
%          September 2012
%******************************************************************

function ID_postWB = whiteBalance(ID_preWB)

% separate out the different bayer planes
ID_plane(:,:,1) = ID_preWB(1:2:end,1:2:end);
ID_plane(:,:,2) = ID_preWB(1:2:end,2:2:end);
ID_plane(:,:,3) = ID_preWB(2:2:end,1:2:end);
ID_plane(:,:,4) = ID_preWB(2:2:end,2:2:end);
[h w c] = size(ID_plane);

blockSizeR = 100; blockSizeC = 100;
centre = round([h/2 w/2]);

centreMean = mean(mean(ID_plane(centre(1)-blockSizeR/2+1:centre(1)+blockSizeR/2, centre(2)-blockSizeC/2+1:centre(2)+blockSizeC/2, :)));
%AE_Level = 0.75; bitDepth = 10; AE_Target = AE_Level * 2^bitDepth;
%centreMean = AE_Target / (mean(mean(ID_plane(centre(1)-blockSizeR/2+1:centre(1)+blockSizeR/2,centre(2)-blockSizeC/2+1:centre(2)+blockSizeC/2,:))));
%balance = centreMean;

balance = max(centreMean) / centreMean;

% apply correction factor for each plane
for i=1:c
    ID_plane(:,:,i) = ID_plane(:,:,i) * balance(:,:,i);
end

% recombine into a single RAW plane
ID_postWB = ID_preWB;
ID_postWB(1:2:end,1:2:end) = ID_plane(:,:,1);
ID_postWB(1:2:end,2:2:end) = ID_plane(:,:,2);
ID_postWB(2:2:end,1:2:end) = ID_plane(:,:,3);
ID_postWB(2:2:end,2:2:end) = ID_plane(:,:,4);

%ID_postWB = round(ID_postWB);
end

%%
% bilinearDemosiac takes a raw image and outputs an RGB image
%
% INPUT:
%   ID: raw image data
%   bayerPattern:
%       rggb
%       bggr
%       gbrg
%       grbg
%
% OUTPUT:
%   RGB: RGB image
%
% COMMENTS:
%   Note that the algorithm uses simple bilinear interpolation. Some special considerations:
%       - Pixel values for a colour which is the origin colour have no
%       interpolation
%       - Borders use a truncated interpolation, averaging only the
%       directly adjacent pixels of that colour
%
%  Author: Eugene Lam
%          September 2012
%******************************************************************

function IDrgb = bilinearDemosaic (IDraw, bayerPattern)

[h w c] = size(IDraw);

switch bayerPattern
    case 'rggb'
        red_mask = repmat([1 0; 0 0], h/2+1, w/2+1);
        green_mask = repmat([0 1; 1 0], h/2+1, w/2+1);
        blue_mask = repmat([0 0; 0 1], h/2+1, w/2+1);
        
    case 'bggr'
        blue_mask = repmat([1 0; 0 0], h/2+1, w/2+1);
        green_mask = repmat([0 1; 1 0], h/2+1, w/2+1);
        red_mask = repmat([0 0; 0 1], h/2+1, w/2+1);
        
    case 'gbrg'
        green_mask = repmat([1 0; 0 1], h/2+1, w/2+1);
        blue_mask = repmat([0 1; 0 0], h/2+1, w/2+1);
        red_mask = repmat([0 0; 1 0], h/2+1, w/2+1);
        
    case 'grbg'
        green_mask = repmat([1 0; 0 1], h/2+1, w/2+1);
        red_mask = repmat([0 1; 0 0], h/2+1, w/2+1);
        blue_mask = repmat([0 0; 1 0], h/2+1, w/2+1);
        
    case 'monochrome'
        red_mask = repmat([1 0; 0 0], h/2+1, w/2+1);
        green_mask = repmat([0 1; 1 0], h/2+1, w/2+1);
        blue_mask = repmat([0 0; 0 1], h/2+1, w/2+1);
        
    otherwise
        error('Error: Unrecognized bayer pattern.');
end

red_mask = red_mask(1:h,1:w);
green_mask = green_mask(1:h,1:w);
blue_mask = blue_mask(1:h,1:w);

R=IDraw.*red_mask;
G=IDraw.*green_mask;
B=IDraw.*blue_mask;

% Interpolation for the green at the missing points
G= G + conv2(G, [0 1 0; 1 0 1; 0 1 0],'same');

% Interpolation for the blue at the missing points
B = B + conv2(B, [1 1 1; 1 0 1; 1 1 1], 'same');

% Interpolation for the red at the missing points
R = R + conv2(R, [1 1 1; 1 0 1; 1 1 1], 'same');

% Calculate number of surrounding pixels of colour plane used
Rpix = conv2(red_mask, ones(3,3),'same');
Gpix = conv2(green_mask, ones(3,3),'same');
Gpix = Gpix - green_mask.*Gpix + green_mask;
Bpix = conv2(blue_mask, ones(3,3),'same');
% Appropriate division
R=R./Rpix;
G=G./Gpix;
B=B./Bpix;

IDrgb(:,:,1)=R; IDrgb(:,:,2)=G; IDrgb(:,:,3)=B;
% IDrgb = floor(IDrgb);

end