% function output = LCBRaw_Mono_test(IDFrameSet, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold,field_effective, option_lsc, ls_mat)
clc;
clear;
close all;

bayerFormat = 'rggb';
pedestal = -64;
bitDepth = 10;

threshold = 16;
field_effective = 0.93;
option_lsc = 3;
roiSize = [7 7];

% get all raw images's path
% file_folder = 'D:\Documents\Polight\Temp\2022\0608\DP_220616_1254pcs_Special\LF';
% file_folder = 'D:\Documents\Polight\svn\Banff2\Banff2_EOL_RX\trunk\Doc\RX_LCB_Attenuation\LCB_Particle221107\Particle';
% file_folder = 'D:\Documents\Polight\Temp\2022\1114\LCB_MP_1112Subgraph\LF';
% file_folder = 'D:\Documents\Polight\Temp\2022\1111\LCB_MP_1111Subgraph\1102\HF\False_Check\Low';

% file_folder = 'D:\Documents\Polight\Temp\2022\1117\LCB_MP_1117_Subgraph\LF\raw';
file_folder = 'D:\Documents\Polight\Temp\2022\1118\LCB_FACode_1118_Subgraph\LF\Cube_Undetected\Raw';
% file_folder = 'D:\Documents\Polight\Temp\2022\1117\1122_4_corner_have_particle\LF';

% file_folder = 'D:\Documents\Polight\Temp\2022\1117\LCB_MP_1117_Subgraph\LF\anomaly_detection\False_check\edg_4corner';
file_list = GetFiles(file_folder, ".raw");

for i = 1:length(file_list)
% for i = 2:2
    pathname = strcat( file_list(i).folder, '\' );
    filename = file_list(i).name;
    
    disp(filename);

    % HF/LF
%     if convertCharsToStrings(file_last_folder) == "HF"
    if contains(filename,'HF') % use filename to judge HF or LF, Sam, 220512
%        fid = fid_HF; % sperate summary file, Sam, 220427
        filterWidth = 5;
%         continue;
%     elseif convertCharsToStrings(file_last_folder) == "LF"
    elseif contains(filename,'LF') % use filename to judge HF or LF, Sam, 220512
%        fid = fid_LF; % sperate summary file, Sam, 220427
        filterWidth = 9;
%         continue;
    end
    
    IDFrameSet = ImportRawFiles(filename, pathname);
    
    ls_mat = [];
    LCBraw_monochrome_LSC_DOE_DPBU_Classify(0, sprintf('%04d_', i), filename, IDFrameSet, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold, field_effective, option_lsc, ls_mat);
%     LCBraw_monochrome_LSC_DOE_ABC(filename, IDFrameSet, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold,field_effective, option_lsc, ls_mat);
end



