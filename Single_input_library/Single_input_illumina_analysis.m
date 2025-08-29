clear
clc
%The illumina data contains bin/sample specific barcodes in reads 1 and 2,
%which we will extract, followed by extracting the variant specific
%barcodes and their counts in each of the sample bins

barcodes = ['AAGG'; 'AACC'; 'TTGG'; 'TTCC'; 'GGAA'; 'GGTT'; 'CCAA'; 'CCTT'];% 'GAGA'; 'TCTC'];

%Initialize variables
forward_ID = zeros(10^9, 1); reverse_ID = forward_ID;
counter = 0;
bc1 = uint8(zeros(10^9, 6)); bc2 = bc1;

%Open files (provided upon request)
fid = fopen('100k-Spect_S1_L002_R1_001.fastq');
fid2 = fopen('100k-Spect_S1_L002_R2_001.fastq');

tline = fgetl(fid); tline = fgetl(fid);
tline2 = fgetl(fid2); tline2 = fgetl(fid2);

while ischar(tline)
    if sum(tline(1) == 'ATGCN') > 0
        counter = counter + 1;
            if mod(counter, 10^6) == 0
               disp(counter)
            end
        x = find(sum((tline2(3:6) == barcodes(1:5, :))') == 4);
        if isempty(x)
            y = find(sum((tline2(4:7) == barcodes(6:8, :))') == 4);
                if isempty(y)
                   forward_ID(counter) = 0;
                else
                   forward_ID(counter) = y + 5;
                end
        else
            forward_ID(counter) = x;
        end

        x = find(sum((tline(5:8) == barcodes([1, 2, 6, 7], :))') == 4);
        if isempty(x)
            reverse_ID(counter) = 0;
        else
            reverse_ID(counter) = x;
        end
        tline = seqrcomplement(tline);
        %Find Barcode 1 sequence
        y = strfind(tline, 'AAACG'); z = strfind(tline, 'CAGTCA');
        if ~isempty(y) & ~isempty(z)
            y = y(1); z = z(1);
            l = z - y - 5;
            if l == 18
                b = tline(y+5:z-1);
                b = b([1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17]);
                if contains(b, 'A')
                    b = 'X';
                else
                    b(b == 'T') = '1';
                    b(b == 'G') = '2';
                    b(b == 'C') = '3';
                
                    bc1(counter, 1) = uint8(str2double(b(1:2)));
                    bc1(counter, 2) = uint8(str2double(b(3:4)));
                    bc1(counter, 3) = uint8(str2double(b(5:6)));
                    bc1(counter, 4) = uint8(str2double(b(7:8)));
                    bc1(counter, 5) = uint8(str2double(b(9:10)));
                    bc1(counter, 6) = uint8(str2double(b(11:12)));
                end
            else
                bc1(counter, :) = 0;
            end
        end

        %Find Barcode 2 sequence
        y = strfind(tline, 'GGTAC'); z = strfind(tline, 'GATAC');
        if ~isempty(y) & ~isempty(z)
            y = y(1); z = z(1);
            l = z - y - 5;
            if l == 18
                b = tline(y+5:z-1);
                b = b([1, 2, 4, 5, 7, 8, 10, 11, 13, 14, 16, 17]);
                if contains(b, 'C')
                    b = 'X';
                else
                    b(b == 'T') = '1';
                    b(b == 'G') = '2';
                    b(b == 'A') = '3';
                
                    bc2(counter, 1) = uint8(str2double(b(1:2)));
                    bc2(counter, 2) = uint8(str2double(b(3:4)));
                    bc2(counter, 3) = uint8(str2double(b(5:6)));
                    bc2(counter, 4) = uint8(str2double(b(7:8)));
                    bc2(counter, 5) = uint8(str2double(b(9:10)));
                    bc2(counter, 6) = uint8(str2double(b(11:12)));
                end
            else
                bc2(counter, :) = 0;
            end
        end

    end
    tline = fgetl(fid); tline = fgetl(fid); tline = fgetl(fid); tline = fgetl(fid);
    tline2 = fgetl(fid2); tline2 = fgetl(fid2); tline2 = fgetl(fid2); tline2 = fgetl(fid2);
end

fclose(fid); fclose(fid2);

%% Convert barcodes to consolidated 32-bit numeric (data compression)

%Trim variables
bc1 = bc1(1:counter, :);
bc2 = bc2(1:counter, :);
forward_ID = forward_ID(1:counter);
reverse_ID = reverse_ID(1:counter);

%Convert each dinucleotide into a single integer (1-9), for BCs 1 and 2
bc1(bc1 == 11) = 1; bc1(bc1 == 12) = 2; bc1(bc1 == 13) = 3; 
bc1(bc1 == 21) = 4; bc1(bc1 == 22) = 5; bc1(bc1 == 23) = 6; 
bc1(bc1 == 31) = 7; bc1(bc1 == 32) = 8; bc1(bc1 == 33) = 9; 

bc2(bc2 == 11) = 1; bc2(bc2 == 12) = 2; bc2(bc2 == 13) = 3; 
bc2(bc2 == 21) = 4; bc2(bc2 == 22) = 5; bc2(bc2 == 23) = 6; 
bc2(bc2 == 31) = 7; bc2(bc2 == 32) = 8; bc2(bc2 == 33) = 9; 

%Trim any incorrect barcodes
bc1_analyze = bc1((sum(bc1, 2) > 0) & (sum(bc2, 2) > 0) & forward_ID > 0 & reverse_ID > 0, :);
bc2_analyze = bc2((sum(bc1, 2) > 0) & (sum(bc2, 2) > 0) & forward_ID > 0 & reverse_ID > 0, :);
fid_analyze = forward_ID((sum(bc1, 2) > 0) & (sum(bc2, 2) > 0) & forward_ID > 0 & reverse_ID > 0);
rid_analyze = reverse_ID((sum(bc1, 2) > 0) & (sum(bc2, 2) > 0) & forward_ID > 0 & reverse_ID > 0);

%Convert forward and reverse IDs into a singular sample bin
read_id = (rid_analyze - 1).*8 + fid_analyze;
read_id = uint8(read_id);

%Combine BCs 1 and 2 into a single BC1+2 number
combined_bc = zeros(size(bc1_analyze, 1), 1);
bc1_conv = uint32(zeros(size(bc1_analyze, 1), 1));
bc2_conv = uint32(zeros(size(bc1_analyze, 1), 1));

for i = 1:size(bc1_analyze, 1)
    if mod(i, 10^6) == 0
        disp(i)
    end
    x = sprintf('%d', bc1_analyze(i, :));
    y = sprintf('%d', bc2_analyze(i, :));
    combined_bc(i) = str2double([x, y]); bc1_conv(i) = uint32(str2double(x)); bc2_conv(i) = uint32(str2double(y));
end

%Extract unique barcode 1 counts, sample specific counts, and normalize for
%read depth differences
[bc1_unique, ~, b] = unique(bc1_conv);
bc1_counts = accumarray(b,1);
bc1_dist = zeros(length(bc1_unique), 32);

for i = 1:32
    disp(i)
    x = accumarray(b(read_id == i), 1);
    bc1_dist(1:length(x), i) = x;
end

%% Compute average expression for each barcode

AB_bin_percs = [3.9 28.8 24.6 20.7 10 3.9 2.6 1.65; 2.6 14.6 15.4 19.2 16.42 9.3 12 7];
BC_bin_percs = [5.8 35.2 25.2 18.7 6.4 2.8 2.3 0.37; 3 14 14 20 17.6 10 13.5 4.7];

bc1_dist = 10^8.*bc1_dist./sum(bc1_dist);

bc1_aboff_average = sum(bc1_dist(:, 1:8).*repmat(1:8, length(bc1_counts), 1).*repmat(AB_bin_percs(1, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 1:8).*AB_bin_percs(1, :), 2);
bc1_abon_average = sum(bc1_dist(:, 17:24).*repmat(1:8, length(bc1_counts), 1).*repmat(AB_bin_percs(2, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 17:24).*AB_bin_percs(2, :), 2);

bc1_bcoff_average = sum(bc1_dist(:, 9:16).*repmat(1:8, length(bc1_counts), 1).*repmat(BC_bin_percs(1, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 9:16).*BC_bin_percs(1, :), 2);
bc1_bcon_average = sum(bc1_dist(:, 25:32).*repmat(1:8, length(bc1_counts), 1).*repmat(BC_bin_percs(2, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 25:32).*BC_bin_percs(2, :), 2);

%Convert bins to expression
%This is done using a fit line obtained by plotting the mid points of all
%sorted bins against the bin number. The output from this function
%corresponds to the log10 expression level obtained from CLASSIC
bin_to_expression = @(x)(10.^(0.463.*x + 1.592));

bc1_aboff_exp = bin_to_expression(bc1_aboff_average);% - 100;
bc1_abon_exp = bin_to_expression(bc1_abon_average);% - 100;
bc1_bcoff_exp = bin_to_expression(bc1_bcoff_average);% - 100;
bc1_bcon_exp = bin_to_expression(bc1_bcon_average);% - 100;

%All Barcodes (BC1 + BC2 combined)

[bcall_unique, ~, b] = unique(combined_bc);
bcall_counts = accumarray(b,1);
bcall_dist = zeros(length(bcall_unique), 32);

for i = 1:32
    disp(i)
    x = accumarray(b(read_id == i), 1);
    bcall_dist(1:length(x), i) = x;
end

bcall_dist = 10^8.*bcall_dist./sum(bcall_dist);

bcall_aboff_average = sum(bcall_dist(:, 1:8).*repmat(1:8, length(bcall_counts), 1).*repmat(AB_bin_percs(1, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 1:8).*AB_bin_percs(1, :), 2);
bcall_abon_average = sum(bcall_dist(:, 17:24).*repmat(1:8, length(bcall_counts), 1).*repmat(AB_bin_percs(2, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 17:24).*AB_bin_percs(2, :), 2);
bcall_bcoff_average = sum(bcall_dist(:, 9:16).*repmat(1:8, length(bcall_counts), 1).*repmat(BC_bin_percs(1, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 9:16).*BC_bin_percs(1, :), 2);
bcall_bcon_average = sum(bcall_dist(:, 25:32).*repmat(1:8, length(bcall_counts), 1).*repmat(BC_bin_percs(2, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 25:32).*BC_bin_percs(2, :), 2);

bcall_aboff_exp = bin_to_expression(bcall_aboff_average);% - 100;
bcall_abon_exp = bin_to_expression(bcall_abon_average);% - 100;
bcall_bcoff_exp = bin_to_expression(bcall_bcoff_average);% - 100;
bcall_bcon_exp = bin_to_expression(bcall_bcon_average);% - 100;

bcall_off_exp = 10.^(0.4514*bcall_off_average + 1.9673);% - 100;
bcall_on_exp = 10.^(0.4514*bcall_on_average + 1.9673);% - 100;

%Now do the same for BC2

[bc2_unique, ~, b] = unique(bc2_conv);
bc2_counts = accumarray(b,1);
bc2_dist = zeros(length(bc2_unique), 32);

for i = 1:32
    disp(i)
    x = accumarray(b(read_id == i), 1);
    bc2_dist(1:length(x), i) = x;
end

bc2_dist = 10^8.*bc2_dist./sum(bc2_dist);

bc2_aboff_average = sum(bc2_dist(:, 1:8).*repmat(1:8, length(bc2_counts), 1).*repmat(AB_bin_percs(1, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 1:8).*AB_bin_percs(1, :), 2);
bc2_abon_average = sum(bc2_dist(:, 17:24).*repmat(1:8, length(bc2_counts), 1).*repmat(AB_bin_percs(2, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 17:24).*AB_bin_percs(2, :), 2);
bc2_bcoff_average = sum(bc2_dist(:, 9:16).*repmat(1:8, length(bc2_counts), 1).*repmat(BC_bin_percs(1, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 9:16).*BC_bin_percs(1, :), 2);
bc2_bcon_average = sum(bc2_dist(:, 25:32).*repmat(1:8, length(bc2_counts), 1).*repmat(BC_bin_percs(2, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 25:32).*BC_bin_percs(2, :), 2);

bc2_aboff_exp = bin_to_expression(bc2_aboff_average);% - 100;
bc2_abon_exp = bin_to_expression(bc2_abon_average);% - 100;
bc2_bcoff_exp = bin_to_expression(bc2_bcoff_average);% - 100;
bc2_bcon_exp = bin_to_expression(bc2_bcon_average);% - 100;

%% Compare with nanopore data - AB first

%"From 100k_SpectR_BC_Nanopore_Export.mat", Load "BC1_conv", "BC2_conv",
%and "barcoded_variants" (variables provided upon request)
AB_100k_variants = barcoded_variants; clear barcoded_variants
AB_reporter_assign_final = AB_100k_variants(:, 6:8);
AB_synTF_assign_final = AB_100k_variants(:, 1:5);
BC1_abfinal_final = BC1_conv; BC2_abfinal_final = BC2_conv;

reporter_assign_final_final = AB_reporter_assign_final(sum(AB_reporter_assign_final > 0, 2) == size(AB_reporter_assign_final, 2) & sum(AB_synTF_assign_final > 0, 2) == size(AB_synTF_assign_final, 2), :);
synTF_assign_final_final = AB_synTF_assign_final(sum(AB_reporter_assign_final > 0, 2) == size(AB_reporter_assign_final, 2) & sum(AB_synTF_assign_final > 0, 2) == size(AB_synTF_assign_final, 2), :);
BC1_conv = BC1_abfinal_final(sum(AB_reporter_assign_final > 0, 2) == size(AB_reporter_assign_final, 2) & sum(AB_synTF_assign_final > 0, 2) == size(AB_synTF_assign_final, 2), :);
BC2_conv = BC2_abfinal_final(sum(AB_reporter_assign_final > 0, 2) == size(AB_reporter_assign_final, 2) & sum(AB_synTF_assign_final > 0, 2) == size(AB_synTF_assign_final, 2), :);

BC1_abfinal_final = BC1_conv; BC2_abfinal_final = BC2_conv;
AB_reporter_assign_final = reporter_assign_final_final;
AB_synTF_assign_final = synTF_assign_final_final;
clear BC1_conv; clear BC2_conv; clear reporter_assign_final_final;
clear synTF_assign_final_final

eu_ABassign_final = AB_synTF_assign_final(:, 1) + 4.*(AB_synTF_assign_final(:, 2) - 1) + 16.*(AB_synTF_assign_final(:, 3) - 1) + 80.*(AB_synTF_assign_final(:, 4) - 1) + 240.*(AB_synTF_assign_final(:, 5) - 1) + 2880*(AB_reporter_assign_final(:, 1) - 1) + 11520*(AB_reporter_assign_final(:, 2) - 1) + 34560*(AB_reporter_assign_final(:, 3) - 1);

[un_nanoporeabBC1, ~, b] = unique(BC1_abfinal_final);
BC1ab_counts_mapping = accumarray(b, 1);

for i = 1:length(BC1ab_counts_mapping)
    j = eu_ABassign_final(b == i);
    BC1ab_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

[un_nanoporeabBC2, ~, b] = unique(BC2_abfinal_final);
BC2ab_counts_mapping = accumarray(b, 1);

for i = 1:length(BC2ab_counts_mapping)
    j = eu_ABassign_final(b == i);
    BC2ab_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode 2 is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

BCallab_final_final = zeros(length(BC2_abfinal_final), 1);

for i = 1:length(BCallab_final_final)
    if mod(i, 10000) == 0
        disp(i)
    end
    BCallab_final_final(i) = str2double([num2str(BC1_abfinal_final(i)) num2str(BC2_abfinal_final(i))]);
end

[un_nanoporeabBCall, ~, b] = unique(BCallab_final_final);
BCall_abcounts_mapping = accumarray(b, 1);

for i = 1:length(BCall_abcounts_mapping)
    j = eu_ABassign_final(b == i);
    BCall_abcounts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

AB_off_nanopore = zeros(size(AB_reporter_assign_final, 1), 1);
AB_on_nanopore = zeros(size(AB_reporter_assign_final, 1), 1);
AB_read_char = zeros(size(AB_reporter_assign_final, 1), 3);

for i = 1:length(AB_off_nanopore)
    if mod(i, 10000) == 0
       disp(i)
    end
    x = bcall_aboff_exp2(bcall_unique2 == BCallab_final_final(i), :);
    y = bcall_abon_exp2(bcall_unique2 == BCallab_final_final(i), :);
    if ~isempty(x)
        AB_read_char(i, 1) = bcall_counts2(bcall_unique2 == BCallab_final_final(i));
        AB_read_char(i, 2) = 12;
        AB_read_char(i, 3) = BCall_abcounts_mapping(un_nanoporeabBCall == BCallab_final_final(i), 2);
    else
        x = bc1_aboff_exp(bc1_unique == BC1_abfinal_final(i), :);
        y = bc1_abon_exp(bc1_unique == BC1_abfinal_final(i), :);
        if ~isempty(x)
            AB_read_char(i, 1) = bc1_counts(bc1_unique == BC1_abfinal_final(i));
            AB_read_char(i, 2) = 1;
            AB_read_char(i, 3) = BC1ab_counts_mapping(un_nanoporeabBC1 == BC1_abfinal_final(i), 2);
        else
            x = bc2_aboff_exp(bc2_unique == BC2_abfinal_final(i), :);
            y = bc2_abon_exp(bc2_unique == BC2_abfinal_final(i), :);
            if ~isempty(x)
                AB_read_char(i, 1) = bc2_counts(bc2_unique == BC2_abfinal_final(i));
                AB_read_char(i, 2) = 2;
                AB_read_char(i, 3) = BC2ab_counts_mapping(un_nanoporeabBC2 == BC2_abfinal_final(i), 2);
            else
                x = 0; y = 0;
            end
        end
    end
    AB_off_nanopore(i) = x; AB_on_nanopore(i) = y;
end

%% Compare with nanopore data - BC

%"From 100k_SpectR_BC_Nanopore_Export.mat", Load "BC1_conv", "BC2_conv",
%and "barcoded_variants" (variables provided upon request)
BC_100k_variants = barcoded_variants; clear barcoded_variants
BC_reporter_assign_final = BC_100k_variants(:, 6:8);
BC_synTF_assign_final = BC_100k_variants(:, 1:5);
BC1_bcfinal_final = BC1_conv; clear BC1_conv
BC2_bcfinal_final = BC2_conv; clear BC2_conv

reporter_assign_final_final = BC_reporter_assign_final(sum(BC_reporter_assign_final > 0, 2) == size(BC_reporter_assign_final, 2) & sum(BC_synTF_assign_final > 0, 2) == size(BC_synTF_assign_final, 2), :);
synTF_assign_final_final = BC_synTF_assign_final(sum(BC_reporter_assign_final > 0, 2) == size(BC_reporter_assign_final, 2) & sum(BC_synTF_assign_final > 0, 2) == size(BC_synTF_assign_final, 2), :);
BC1_conv = BC1_bcfinal_final(sum(BC_reporter_assign_final > 0, 2) == size(BC_reporter_assign_final, 2) & sum(BC_synTF_assign_final > 0, 2) == size(BC_synTF_assign_final, 2), :);
BC2_conv = BC2_bcfinal_final(sum(BC_reporter_assign_final > 0, 2) == size(BC_reporter_assign_final, 2) & sum(BC_synTF_assign_final > 0, 2) == size(BC_synTF_assign_final, 2), :);

BC1_bcfinal_final = BC1_conv; BC2_bcfinal_final = BC2_conv;
BC_reporter_assign_final = reporter_assign_final_final;
BC_synTF_assign_final = synTF_assign_final_final;
clear BC1_conv; clear BC2_conv; clear reporter_assign_final_final;
clear synTF_assign_final_final

eu_BCassign_final = BC_synTF_assign_final(:, 1) + 4.*(BC_synTF_assign_final(:, 2) - 1) + 16.*(BC_synTF_assign_final(:, 3) - 1) + 80.*(BC_synTF_assign_final(:, 4) - 1) + 240.*(BC_synTF_assign_final(:, 5) - 1) + 2880*(BC_reporter_assign_final(:, 1) - 1) + 11520*(BC_reporter_assign_final(:, 2) - 1) + 34560*(BC_reporter_assign_final(:, 3) - 1);

[un_nanoporebcBC1, ~, b] = unique(BC1_bcfinal_final);
BC1bc_counts_mapping = accumarray(b, 1);

for i = 1:length(BC1bc_counts_mapping)
    j = eu_BCassign_final(b == i);
    BC1bc_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

[un_nanoporebcBC2, ~, b] = unique(BC2_bcfinal_final);
BC2bc_counts_mapping = accumarray(b, 1);

for i = 1:length(BC2bc_counts_mapping)
    j = eu_BCassign_final(b == i);
    BC2bc_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode 2 is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

BCallbc_final_final = zeros(length(BC2_bcfinal_final), 1);

for i = 1:length(BCallbc_final_final)
    if mod(i, 10000) == 0
        disp(i)
    end
    BCallbc_final_final(i) = str2double([num2str(BC1_bcfinal_final(i)) num2str(BC2_bcfinal_final(i))]);
end

[un_nanoporebcBCall, ~, b] = unique(BCallbc_final_final);
BCall_bccounts_mapping = accumarray(b, 1);

for i = 1:length(BCall_bccounts_mapping)
    j = eu_BCassign_final(b == i);
    BCall_bccounts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

BC_off_nanopore = zeros(size(BC_reporter_assign_final, 1), 1);
BC_on_nanopore = zeros(size(BC_reporter_assign_final, 1), 1);
BC_read_char = zeros(size(BC_reporter_assign_final, 1), 3);

for i = 1:length(BC_off_nanopore)
    if mod(i, 10000) == 0
       disp(i)
    end
    x = bcall_bcoff_exp(bcall_unique == BCallbc_final_final(i), :);
    y = bcall_bcon_exp(bcall_unique == BCallbc_final_final(i), :);
    if ~isempty(x)
        BC_read_char(i, 1) = bcall_counts(bcall_unique == BCallbc_final_final(i));
        BC_read_char(i, 2) = 12;
        BC_read_char(i, 3) = BCall_bccounts_mapping(un_nanoporebcBCall == BCallbc_final_final(i), 2);
    else
        x = bc1_bcoff_exp(bc1_unique == BC1_bcfinal_final(i), :);
        y = bc1_bcon_exp(bc1_unique == BC1_bcfinal_final(i), :);
        if ~isempty(x)
            BC_read_char(i, 1) = bc1_counts(bc1_unique == BC1_bcfinal_final(i));
            BC_read_char(i, 2) = 1;
            BC_read_char(i, 3) = BC1bc_counts_mapping(un_nanoporebcBC1 == BC1_bcfinal_final(i), 2);
        else
            x = bc2_bcoff_exp(bc2_unique == BC2_bcfinal_final(i), :);
            y = bc2_bcon_exp(bc2_unique == BC2_bcfinal_final(i), :);
            if ~isempty(x)
                BC_read_char(i, 1) = bc2_counts(bc2_unique == BC2_bcfinal_final(i));
                BC_read_char(i, 2) = 2;
                BC_read_char(i, 3) = BC2bc_counts_mapping(un_nanoporebcBC2 == BC2_bcfinal_final(i), 2);
            else
                x = 0; y = 0;
            end
        end
    end
    BC_off_nanopore(i) = x; BC_on_nanopore(i) = y;
end

%% Map illumina data to variants using the barcodes from illumina and nanopore indexing

ordered_eu_expAB = zeros(103680, 3);
ordered_eu_expBC = zeros(103680, 3);

ABeu_ordered_assign_split = zeros(103680, 8);
ABeu_ordered_assign_split(:, 1) = repmat((1:4)', 103680/4, 1);
ABeu_ordered_assign_split(:, 2) = repmat([ones(1, 4), 2*ones(1, 4), 3*ones(1, 4), 4*ones(1, 4)]', 103680/16, 1);
ABeu_ordered_assign_split(:, 3) = repmat((1:4)', 103680/4, 1);
ABeu_ordered_assign_split(:, 4) = repmat((1:4)', 103680/4, 1);
ABeu_ordered_assign_split(:, 5) = repmat((1:4)', 103680/4, 1);

for i = 1:3
    for j = 1:3
        for k = 1:4
            for l = 1:12
                for m = 1:3
                    for n = 1:5
                        for o = 1:4
                            for p = 1:4
                                ABeu_ordered_assign_split(p + (o-1)*4 + (n-1)*16 + (m-1)*80 + (l-1)*240 + (k-1)*2880 + (j-1)*11520 + (i-1)*34560, :) = [p o n m l k j i];
                            end
                        end
                    end
                end
            end
        end
    end
end

for i = 1:size(ordered_eu_expAB, 1)
    if mod(i, 1000) == 0
        disp(i)
    end
    c = AB_off_nanopore(eu_ABassign_final == i & AB_read_char(:, 3) == 1);
    d = AB_on_nanopore(eu_ABassign_final == i & AB_read_char(:, 3) == 1);
    a = c(~isnan(c)); b = d(~isnan(d));
    if ~isempty(a) && ~isempty(b)
            a = mean(a);
            b = mean(b);
            if isempty(b)
                b = 0;
            end
            c = b/a;
            ordered_eu_expAB(i, :) = [a b c];
    end    
end

ordered_eu_expAB(ordered_eu_expAB(:, 1) > 0, 3) = (ordered_eu_expAB(ordered_eu_expAB(:, 1) > 0, 2) - 200)./(ordered_eu_expAB(ordered_eu_expAB(:, 1) > 0, 1) - 200);

BCeu_ordered_assign_split = ABeu_ordered_assign_split;

for i = 1:size(ordered_eu_expBC, 1)
    if mod(i, 1000) == 0
        disp(i)
    end
    c = BC_off_nanopore(eu_BCassign_final == i & BC_read_char(:, 3) == 1);
    d = BC_on_nanopore(eu_BCassign_final == i & BC_read_char(:, 3) == 1);
    a = c(~isnan(c) & ~isnan(d)); b = d(~isnan(c) & ~isnan(d));
    if ~isempty(a) && ~isempty(b)
            a = mean(a);
            b = mean(b);
            if isempty(b)
                b = 0;
            end
            c = b/a;
            ordered_eu_expBC(i, :) = [a b c];
    end    
end

ordered_eu_expBC(ordered_eu_expBC(:, 1) > 0, 3) = (ordered_eu_expBC(ordered_eu_expBC(:, 1) > 0, 2) - 200)./(ordered_eu_expBC(ordered_eu_expBC(:, 1) > 0, 1) - 200);
