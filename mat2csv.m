function mat2csv(path, A)
[row col] = size(A);
fid = fopen(path, 'w');

for i = 1:row
    for j = 1:col-1
        fprintf(fid, '%.8f\t', A(i, j));
    end
    fprintf(fid, '%.8f\r\n', A(i, col));
end

fclose(fid);
end
