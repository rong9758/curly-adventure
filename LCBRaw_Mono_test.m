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
% roiSize = [5 5]; % preserve more information to help judge particle or noise, Sam, 220519
% roiSize = [6 6]; % preserve more information to help judge particle or noise, Sam, 220519
roiSize = [7 7];

% get all raw images's path
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\1.POR_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\2.+30(848)_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\3.+60(878)_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\4.+90(908)_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\5.-30(788)_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\6.-60(758)_APC'
%file_folder = 'd:\Documents\Polight\svn\Sphinx\Sphinx_EOL_RX\trunk\Doc\Requirement\221102_LCBraw_AE_Target_Change\7.-90(728)_APC'
% file_folder = 'd:\Documents\Polight\svn\DpbuSW\18_Staff_Meeting\2022\220823\LCB_RAW_Defect';
% file_folder = 'd:\Documents\Polight\svn\DpbuSW\18_Staff_Meeting\2022\220823\LCB_RAW_Defect_2';
%file_folder = 'd:\Documents\Polight\svn\DpbuSW\18_Staff_Meeting\2022\220823\LCB_RAW_Defect_3';
%file_folder = 'd:\Documents\Polight\Temp\2022\0608\DP_220608_1254pcs';
file_folder = 'C:\Users\G1652925\Desktop\LCB\Have_Particle';
% file_folder = 'd:\Documents\Polight\Temp\2022\0608\DP_220616_1254pcs_Special';
% file_folder = 'd:\Documents\Polight\Temp\2022\0608\DP_220610_1254pcs_Special';
% file_folder = 'D:\Documents\Polight\Temp\2022\0608\DP_220608_1254pcs_Special';
% file_folder = 'd:\Documents\Polight\Temp\2022\0608\DP_220608_1254pcs_Test';
% file_folder = 'D:\Documents\Polight\Temp\2022\0608\DP_220608_1254pcs';
% file_folder = 'd:\Documents\Polight\Temp\2022\0526\DP_220526_939pcs_Special';
% file_folder = 'd:\Documents\Polight\Temp\2022\0525\DP_Special';
% file_folder = 'D:\Documents\Polight\Temp\2022\0527\DP_220527_320pcs';
% file_folder = 'd:\Documents\Polight\Temp\2022\0527\DP_220527_320pcs_Special';
% file_folder = 'D:\Documents\Polight\Temp\2022\0526\DP_220526_939pcs';
% file_folder = 'D:\Documents\Polight\Temp\2022\0525\DP_220525_320pcs';
% file_folder = 'D:\Documents\Polight\Temp\2022\0518\DP_220518_618pcs';
% file_folder = 'D:\Documents\Polight\Temp\2022\0518\DP_EOL_240pcs';
% file_folder = 'D:\Documents\Polight\Temp\2022\0518\DP_OQC_179pcs';
% file_folder = 'D:\Documents\Polight\Temp\2022\0509\DP_200pcs';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_10pcs';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Pass_3pcs';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\DP_20pcs';

% debug, Sam, 2205125
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\Debug';

% 1st location, hard
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_10pcs\F8P209404XD15F2AC';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_10pcs\F8P209405C115F2AC';

% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Pass_3pcs\F8P209101PP15F2AY';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_10pcs\F8P2091018115F2AX';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_4pcs_Hard';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_4pcs_Hard\F8P209101PP15F2AY';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_4pcs_Hard\F8P2091003M15F2AT';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_4pcs_Hard\F8P2094090R15F2AL';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Margin_4pcs_Hard\F8P2091018115F2AX';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\LG_Pass_3pcs\F8P209100HZ15F2A7';
% file_folder = 'd:\Documents\Polight\Temp\2022\0422\DP_20pcs\F8P203603KS15F27T';
file_list = GetFiles(file_folder, ".raw");

slash_idx = strfind(file_folder, '\'); % Sam, 220414
file_last_folder = file_folder(slash_idx(end)+1:end); % the last folder name, Sam, 220414
% output_summary = strcat(file_last_folder, '_results.csv');
output_summary_LF = strcat(file_last_folder, '_LF.csv');
output_summary_HF = strcat(file_last_folder, '_HF.csv');
fid_LF = fopen(output_summary_LF, 'w');
fid_HF = fopen(output_summary_HF, 'w');

title_1 = sprintf('%s,%s,%s,%s,%s,', 'filename', 'col', 'row', 'col_', 'row_');
title_2 = sprintf('%s,%s,%s,%s,', ...
                  'IDbin_o_roi_par_avg', 'IDbin_o_roi_bor_avg', 'IDbin_hist_roi_par_avg', 'IDbin_hist_roi_bor_avg');
%         LCBraw_center LCBraw_center_x LCBraw_center_y,...         % original LCB scores
%         LCBraw_edge   LCBraw_edge_x   LCBraw_edge_y,...
%         LCBraw_corner LCBraw_corner_x LCBraw_corner_y,...
title_3 = sprintf('%s,%s,', 'LCBraw_center_x_DP', 'LCBraw_center_y_DP');
title_4 = sprintf('%s,%s,', 'LCBraw_edge_x_DP', 'LCBraw_edge_y_DP');
title_5 = sprintf('%s,%s,', 'LCBraw_corner_x_DP', 'LCBraw_corner_y_DP');
title_6 = sprintf('%s,', 'attenuation_DP');
title_7 = sprintf('%s,%s,', 'LCBraw_center_x_ABC', 'LCBraw_center_y_ABC');
title_8 = sprintf('%s,%s,', 'LCBraw_edge_x_ABC', 'LCBraw_edge_y_DP');
title_9 = sprintf('%s,%s,', 'LCBraw_corner_x_ABC', 'LCBraw_corner_y_ABC');
title_10 = sprintf('%s\r\n', 'attenuation_ABC');
% fprintf(fid_LF, '%s%s%s%s%s%s%s%s', title_1, title_2, title_3, title_4, title_5, title_6, title_7, title_8);
% fprintf(fid_HF, '%s%s%s%s%s%s%s%s', title_1, title_2, title_3, title_4, title_5, title_6, title_7, title_8);
fprintf(fid_LF, '%s%s%s%s%s%s%s%s%s%s', title_1, title_2, title_3, title_4, title_5, title_6, title_7, title_8, title_9, title_10);
fprintf(fid_HF, '%s%s%s%s%s%s%s%s%s%s', title_1, title_2, title_3, title_4, title_5, title_6, title_7, title_8, title_9, title_10);
% fclose(fid_LF); % sperate summary file, Sam, 220427
% fclose(fid_HF); % sperate summary file, Sam, 220427

% fclose(fid_LF); % sperate summary file, Sam, 220427
for i = 1:length(file_list)
% for i = 2:2
    pathname = strcat( file_list(i).folder, '\' );
    filename = file_list(i).name;
    
    disp(filename);
    
%     if convertCharsToStrings(filename) == "BAQ_EOLR_01_C3066_0x2704E091_F8P2091003M15F2AT_20220410154957_LCBRawHF_BB_APC.raw"
%         filename = filename;
%     end
    
    % HF/LF
    slash_idx = strfind(file_list(i).folder, '\'); % Sam, 220414
    file_last_folder = file_list(i).folder(slash_idx(end)+1:end); % the last folder name, Sam, 220414
%     if convertCharsToStrings(file_last_folder) == "HF"
    if contains(filename,'HF') % use filename to judge HF or LF, Sam, 220512
        fid = fid_HF; % sperate summary file, Sam, 220427
        filterWidth = 5;
%         continue;
%     elseif convertCharsToStrings(file_last_folder) == "LF"
    elseif contains(filename,'LF') % use filename to judge HF or LF, Sam, 220512
        fid = fid_LF; % sperate summary file, Sam, 220427
        filterWidth = 9;
%         continue;
    end
    
    IDFrameSet = ImportRawFiles(filename, pathname);
    
    ls_mat = [];
    [output_DPBU] = LCBraw_monochrome_LSC_DOE_DPBU(fid, sprintf('%04d_', i), filename, IDFrameSet, bayerFormat, pedestal, bitDepth, roiSize, filterWidth, threshold, field_effective, option_lsc, ls_mat);
    [output_ABC]  = LCBraw_monochrome_LSC_DOE_ABC(filename, IDFrameSet, bayerFormat, pedestal, bitDepth, [7 7], filterWidth, threshold, field_effective, option_lsc, ls_mat);
%     fprintf(fid, '%s,%s,%s,%s,%s,', filename, num2str(col), num2str(row), num2str(col_), num2str(row_));
%     fprintf(fid, '%d,', length(stats));
%     fprintf(fid, '%d,%d,', area_th_05_num, area_th_10_num);
%     fprintf(fid, '%d,%d,', area_th_15_num, area_th_20_num);
%     fprintf(fid, '%d,%d,', area_th_25_num, area_th_30_num);
%     fprintf(fid, '%d,%.8f,%.8f', hROI*wROI, area_average, area_std); % area of roi, average of area, std of area
    
%         LCBraw_center LCBraw_center_x LCBraw_center_y,...         % original LCB scores
%         LCBraw_edge   LCBraw_edge_x   LCBraw_edge_y,...
%         LCBraw_corner LCBraw_corner_x LCBraw_corner_y,...
    fprintf(fid, '%d,%d,', output_DPBU.testLog(2), output_DPBU.testLog(3));
    fprintf(fid, '%d,%d,', output_DPBU.testLog(5), output_DPBU.testLog(6));
    fprintf(fid, '%d,%d,', output_DPBU.testLog(8), output_DPBU.testLog(9));
%     fprintf(fid, '%d\r\n', output_DPBU.testLog(12));
    fprintf(fid, '%f,', output_DPBU.testLog(12));
    fprintf(fid, '%d,%d,', output_ABC.testLog(2), output_ABC.testLog(3));
    fprintf(fid, '%d,%d,', output_ABC.testLog(5), output_ABC.testLog(6));
    fprintf(fid, '%d,%d,', output_ABC.testLog(8), output_ABC.testLog(9));
%     fprintf(fid, '%d\r\n', output_DPBU.testLog(12));
    fprintf(fid, '%f\r\n', output_ABC.testLog(12));
end

fclose(fid_LF); % sperate summary file, Sam, 220427
fclose(fid_HF); % sperate summary file, Sam, 220427
