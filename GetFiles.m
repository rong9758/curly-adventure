function r_file_list = GetFiles(file_folder, file_extension)

%file_path = 'D:\Documents\Polight\Temp\2022\0407\BAQ_EOLR_01_C3066_0x26EFE7A4_F8P209402GB15F2AT_20220330101656'
file_path_list = dir( strcat(file_folder, '\*.*') );
file_num = length(file_path_list);
r_file_list = [];

if file_num > 0
    for i = 1:file_num
        if file_path_list(i).name=="." || file_path_list(i).name==".."
            continue;
        end
        
        file_path = file_path_list(i).folder + "\" + file_path_list(i).name;
        if contains(file_path_list(i).name, file_extension)==0 % the file without extension, Sam, 220408
            if contains(file_path_list(i).name, ".")==0
                r_file_list = [r_file_list; GetFiles(file_path, file_extension)];
                continue;
            else
                continue;
            end
        end
        
%         if contains(file_path_list(i).name, 'LCBRawHF')>0
%             i=i;
%         end

        r_file_list = [r_file_list; struct('folder', file_path_list(i).folder, 'name', file_path_list(i).name, 'path', file_path)];
    end
end

%r_file_list = file_list;