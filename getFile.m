% getFile opens a file and outputs it to an easy to use format
%
% INPUT:
%   filename: filename of image to be captured
%   filetype: filetype currently supports
%       bmp/jpg
%       raw_cmbu
%       raw_cowell
%       raw_ov
%       raw_sony
%   outputFormat: requested output
%       'raw8'
%       'raw10'
%       'bayer8'
%       'bayer10'
%       'rgb8'
%       'rgb10'
%       'interpolated8'
%       'interpolated10'
%
% OUTPUT:
%   ID: image array as requested
%
%  Author: Eugene Lam (modified by)
%          February 2011
%******************************************************************

function [ID] = getFile(filename, inputFormat, outputFormat)
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%	parameters
%%%%%%%%%%%%%%%%%%%%%%%%%
%      subtractPedestal10bit =64;  % 10-bit
%      subtractPedestal12bit =256;  % 12-bit
    subtractPedestal10bit =0;  % 10-bit
    subtractPedestal12bit =0;  % 12-bit
    signedValue = true;     % whether pedestal subtraction should result in signed
    debug = false;
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%	input format to get common raw that is just numbers
%%%%%%%%%%%%%%%%%%%%%%%%%
    switch lower(inputFormat)
        case {'bmp','jpg'}
            RGB = imread(filename);
            bayerFormat = '';
        case {'raw_cmbu'}
            bayerFormat = 'bggr';
            fid = fopen(filename,'r');
            [H V bpp] = getResolution(fid);
            if (bpp == 1)
                raw8 = fread(fid,[H,V],'*uint8','b')';
                raw10 = 4 * raw8;
            else
                raw16 = fread(fid,[H,V],'*uint16','b')';
                raw10 = bitshift(raw16, -6);
                raw12 = bitshift(raw16, -4);
            end
            fclose(fid);
        case {'raw_cowell'}
            bayerFormat = 'bggr';
            fid = fopen(filename,'r');
            [H V bpp] = getResolution(fid);          
            if (bpp == 1)
                raw8 = fread(fid, [H V], '*uint8','b')';
                raw10 = 4 * raw8;
            else
                raw10 = fread(fid, [H V], '*uint16','l')';
            end
            fclose(fid);
        case {'raw_ov'}
            bayerFormat = 'bggr';
            fid = fopen(filename,'r');
            [H V bpp] = getResolution(fid);
            if (bpp == 1)
                raw8 = fread(fid, [H V], '*uint8','b')';
                raw10 = 4 * raw8;
            else
                raw10 = fread(fid, [H V], '*uint16','b')';
            end
            fclose(fid);
        case {'raw_sony'}
            bayerFormat = 'rggb';
            fid = fopen(filename,'r');
            [H V] = fread(fid, 2, 'uint16', 0, 'l');    % first two bytes indicate H and V
            raw10 = fread(fid,[H V],'*uint16','l')';
        case {'raw_vpt'}
            bayerFormat = 'bggr';
            fid = fopen(filename,'r');
            [H V bpp] = getResolution(fid);          
            if (bpp == 1)
                raw8 = fread(fid, [H V], '*uint8','b')';
                raw10 = 4 * raw8;
            else
                raw10 = fread(fid, [H V], '*uint16','l')';
            end
            fclose(fid);
    end
    raw8 = raw10 / 4;
    
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%	pedestal formatting
%%%%%%%%%%%%%%%%%%%%%%%%%
    if signedValue
        raw10 = int16(raw10) - subtractPedestal10bit;
        raw12 = int16(raw12) - subtractPedestal12bit;
    else
        raw10 = raw10 - subtractPedestal10bit;
        raw12 = raw12 - subtractPedestal12bit;
    end
	raw8 = raw10 / 4;
        
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%	convert from raw to required format
%%%%%%%%%%%%%%%%%%%%%%%%%
    switch lower(outputFormat)
        case {'raw8'}
            ID = raw8;
        case {'raw10'}
            ID = raw10;
        case {'raw12'}
            ID = raw12;
        case {'bayer8'}
            ch1 = raw10(1:2:end,1:2:end);
            ch2 = raw10(1:2:end,2:2:end);
            ch3 = raw10(2:2:end,1:2:end);
            ch4 = raw10(2:2:end,2:2:end); 
            switch (bayerFormat)
                case 'rggb'
                    R=ch1; Gr=ch2; Gb=ch3; B=ch4;
                case 'bggr'
                    B=ch1; Gb=ch2; Gr=ch3; R=ch4;
                case 'gbrg'
                    Gb=ch1; B=ch2; R=ch3; Gr=ch4;
                case 'grbg'
                    Gr=ch1; R=ch2; B=ch3; Gb=ch4;
            end
            ID(:,:,1) = R;
            ID(:,:,2) = Gr;
            ID(:,:,3) = Gb;
            ID(:,:,4) = B;
            ID = ID/4;
            
%             % added by Cathy for cowell raw
%             ID = ID/2^6;
            
        case {'bayer10'}
            ch1 = raw10(1:2:end,1:2:end);
            ch2 = raw10(1:2:end,2:2:end);
            ch3 = raw10(2:2:end,1:2:end);
            ch4 = raw10(2:2:end,2:2:end); 
            switch (bayerFormat)
                case 'rggb'
                    R=ch1; Gr=ch2; Gb=ch3; B=ch4;
                case 'bggr'
                    B=ch1; Gb=ch2; Gr=ch3; R=ch4;
                case 'gbrg'
                    Gb=ch1; B=ch2; R=ch3; Gr=ch4;
                case 'grbg'
                    Gr=ch1; R=ch2; B=ch3; Gb=ch4;
            end
            ID(:,:,1) = R;
            ID(:,:,2) = Gr;
            ID(:,:,3) = Gb;
            ID(:,:,4) = B;
        case {'rgb8'}
            ID = bilinear_demosaic_3(raw8,bayerFormat);
        case {'rgb10'}
            ID = bilinear_demosaic_3(raw10,bayerFormat);
        case {'interpolated8'}
            RGB8 = bilinear_demosaic_3(double(raw8),bayerFormat);
            RGB8 = uint8(RGB8); % add by Cathy for cowell raw
            ID =rgb2ycbcr2(RGB8);
        case {'interpolated10'}
            RGB10 = bilinear_demosaic3(raw10,bayerFormat);
            ID = rgb2ycbcr2(RGB10);
    end
    clear raw8 raw10 raw12 raw16;
end


function [H V bpp] = getResolution(fid)

    fseek(fid,0,'eof');%go to end of file
    filesize=ftell(fid); %current position is length of file
    fseek(fid,0,'bof'); % go back to beginning of file

    %           VGA     720P,                1.3, 2,   1080P 3,   5 ,  5   5     5   8    
    widths   = [640 672 1104 1104 1280 1292 1296 1280 1376 1376 1600 1920 2048 2752 2608 2624 3264 3840 3856];
    heights  = [480 488 1312 1314 800  808  808   960 960  968  1200 1162 1536 1936 1952 1952 2448 2324 2340];
    
    %widths   = [640 672 1280 1280 1292 1296 1280 1376 1600 1920 2048 2752 2592 2608 2624 3264 3280 3296]; 
    %heights  = [480 488 720  800  808  808   960 960  1200 1080 1536  968 1944 1952 1952 2448 2464 2464];
    
    %widths   = [640 672 1280 1280 1292 1296 1280 1376 1600 1920 2048 2752 2592 2608 2624 3264 3280 3296];
    %heights  = [480 488 720  800  808  808   960 960  1200 1080 1536  968 1944 1952 1952 2448 2464 2464];

    resolutions = heights.*widths;
    find_bpp = filesize ./ resolutions;
    resolution_idx = [find(find_bpp == 1.0,1) find(find_bpp == 1.5,1) find(find_bpp == 2.0,1)];

    if ~isempty(resolution_idx)
        bpp = (1.0 * ~isempty(find(find_bpp == 1.0,1))) + ...
              (1.5 * ~isempty(find(find_bpp == 1.5,1))) + ...
              (2.0 * ~isempty(find(find_bpp == 2.0,1)));
        H = widths(resolution_idx);
        V = heights(resolution_idx);
    else
        bpp = 0; % Sam, 220424
        H=0;
        V=0;
    end
    
end