clear
clc
%% Extract barcodes

%The illumina data contains bin/sample specific barcodes in reads 1 and 2,
%which we will extract, followed by extracting the variant specific
%barcodes and their counts in each of the sample bins

%Initialize variables
barcodes = ['AAGG'; 'AACC'; 'TTGG'; 'TTCC'; 'GGAA'; 'GGTT'; 'CCAA'; 'CCTT'];% 'GAGA'; 'TCTC'];
forward_ID = zeros(10^9, 1); %reverse_ID = forward_ID;
reverse_ID = zeros(10^9, 1); %reverse_ID = forward_ID;
counter = 0;
bc1 = uint8(zeros(10^9, 6)); bc2 = bc1;

%Open files (provided upon request)
fid = fopen('Stuff_S1_L001_R1_001.fastq');
fid2 = fopen('Stuff_S1_L001_R2_001.fastq');

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
        y = strfind(tline, 'AAACG'); z = strfind(tline, 'CAGTTA');
        if ~isempty(y) && ~isempty(z)
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
        y = strfind(tline, 'AGGTAC'); z = strfind(tline, 'GATAC');
        if ~isempty(y) && ~isempty(z)
            y = y(1); z = z(1);
            l = z - y - 6;
            if l == 18
                b = tline(y+6:z-1);
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

bin_percs = [3.28 40.34 23.26 8.50 6.76 5.01 4.94 4.11; 3.2	42.55 26.57 9.94 6.12 3.66 2.62 1.66; 4.35 39.59 28.12 9.4 5.72 3.8 3.47 2.33; 6.14 46.71 28.23 8.67 3.93 1.38 0.79 0.41];

bc1_dist = 10^8.*bc1_dist./sum(bc1_dist);

bc1_off_average = sum(bc1_dist(:, 25:32).*repmat(1:8, length(bc1_counts), 1).*repmat(bin_percs(4, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 25:32).*bin_percs(4, :), 2);
bc1_oht_average = sum(bc1_dist(:, 17:24).*repmat(1:8, length(bc1_counts), 1).*repmat(bin_percs(3, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 17:24).*bin_percs(3, :), 2);
bc1_gzv_average = sum(bc1_dist(:, 9:16).*repmat(1:8, length(bc1_counts), 1).*repmat(bin_percs(2, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 9:16).*bin_percs(2, :), 2);
bc1_on_average = sum(bc1_dist(:, 1:8).*repmat(1:8, length(bc1_counts), 1).*repmat(bin_percs(1, :), length(bc1_counts), 1), 2)./sum(bc1_dist(:, 1:8).*bin_percs(1, :), 2);

%Convert bins to expression
%This is done using a fit line obtained by plotting the mid points of all
%sorted bins against the bin number. The output from this function
%corresponds to the log10 expression level obtained from CLASSIC
bin_to_expression = @(x)(10.^(0.463.*x + 1.592));

bc1_off_exp = bin_to_expression(bc1_off_average);
bc1_oht_exp = bin_to_expression(bc1_oht_average);
bc1_gzv_exp = bin_to_expression(bc1_gzv_average);
bc1_on_exp = bin_to_expression(bc1_on_average);

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

bcall_off_average = sum(bcall_dist(:, 25:32).*repmat(1:8, length(bcall_counts), 1).*repmat(bin_percs(4, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 25:32).*bin_percs(4, :), 2);
bcall_oht_average = sum(bcall_dist(:, 17:24).*repmat(1:8, length(bcall_counts), 1).*repmat(bin_percs(3, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 17:24).*bin_percs(3, :), 2);
bcall_gzv_average = sum(bcall_dist(:, 9:16).*repmat(1:8, length(bcall_counts), 1).*repmat(bin_percs(2, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 9:16).*bin_percs(2, :), 2);
bcall_on_average = sum(bcall_dist(:, 1:8).*repmat(1:8, length(bcall_counts), 1).*repmat(bin_percs(1, :), length(bcall_counts), 1), 2)./sum(bcall_dist(:, 1:8).*bin_percs(1, :), 2);

bcall_off_exp = bin_to_expression(bcall_off_average);% - 100;
bcall_oht_exp = bin_to_expression(bcall_oht_average);% - 100;
bcall_gzv_exp = bin_to_expression(bcall_gzv_average);% - 100;
bcall_on_exp = bin_to_expression(bcall_on_average);% - 100;

%BC 2

[bc2_unique, ~, b] = unique(bc2_conv);
bc2_counts = accumarray(b,1);
bc2_dist = zeros(length(bc2_unique), 32);

for i = 1:32
    disp(i)
    x = accumarray(b(read_id == i), 1);
    bc2_dist(1:length(x), i) = x;
end

bc2_dist = 10^8.*bc2_dist./sum(bc2_dist);

bc2_off_average = sum(bc2_dist(:, 25:32).*repmat(1:8, length(bc2_counts), 1).*repmat(bin_percs(4, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 25:32).*bin_percs(4, :), 2);
bc2_oht_average = sum(bc2_dist(:, 17:24).*repmat(1:8, length(bc2_counts), 1).*repmat(bin_percs(3, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 17:24).*bin_percs(3, :), 2);
bc2_gzv_average = sum(bc2_dist(:, 9:16).*repmat(1:8, length(bc2_counts), 1).*repmat(bin_percs(2, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 9:16).*bin_percs(2, :), 2);
bc2_on_average = sum(bc2_dist(:, 1:8).*repmat(1:8, length(bc2_counts), 1).*repmat(bin_percs(1, :), length(bc2_counts), 1), 2)./sum(bc2_dist(:, 1:8).*bin_percs(1, :), 2);

bc2_off_exp = bin_to_expression(bc2_off_average);% - 100;
bc2_oht_exp = bin_to_expression(bc2_oht_average);% - 100;
bc2_gzv_exp = bin_to_expression(bc2_gzv_average);% - 100;
bc2_on_exp = bin_to_expression(bc2_on_average);% - 100;

%% Compare with nanopore data

%from "Reads_all_combined.mat", Load "BC1_conv_all", "BC2_conv_all",
%and "barcoded_variants_all_all" (variables provided upon request)
Sub_1M_variants = barcoded_variants_all_all; clear barcoded_variants_all_all
reporter_ass_final = Sub_1M_variants(:, 1:3);
synTF_ass_final = Sub_1M_variants(:, 4:11);
orientation_final = Sub_1M_variants(:, 12);
BC1_final_final = BC1_conv_all; BC2_final_final = BC2_conv_all;

clear BC1_conv_all; clear BC2_conv_all;

eu_ass_final = Sub_1M_variants(:, 1) + 3.*(Sub_1M_variants(:, 2) - 1) + 12.*(Sub_1M_variants(:, 3) - 1) + 24.*(Sub_1M_variants(:, 4) - 1) + 72.*(Sub_1M_variants(:, 5) - 1) + 288.*(Sub_1M_variants(:, 6) - 1) + 576.*(Sub_1M_variants(:, 7) - 1) + 2304.*(Sub_1M_variants(:, 8) - 1) + 6912.*(Sub_1M_variants(:, 9) - 1) + 27648.*(Sub_1M_variants(:, 10) - 1) + 55296.*(Sub_1M_variants(:, 11) - 1) + 221184.*(Sub_1M_variants(:, 12) - 1);

[un_nanoporeBC1, ~, b] = unique(BC1_final_final);
BC1_counts_mapping = accumarray(b, 1);

for i = 1:length(BC1_counts_mapping)
    j = eu_ass_final(b == i);
    BC1_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

[un_nanoporeBC2, ~, b] = unique(BC2_final_final);
BC2_counts_mapping = accumarray(b, 1);

for i = 1:length(BC2_counts_mapping)
    j = eu_ass_final(b == i);
    BC2_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode 2 is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

BCall_final_final = zeros(length(BC2_final_final), 1);

for i = 1:length(BCall_final_final)
    if mod(i, 10000) == 0
        disp(i)
    end
    BCall_final_final(i) = str2double([num2str(BC1_final_final(i)) num2str(BC2_final_final(i))]);
end

[un_nanoporeBCall, ~, b] = unique(BCall_final_final);
BCall_counts_mapping = accumarray(b, 1);

for i = 1:length(BCall_counts_mapping)
    j = eu_ass_final(b == i);
    BCall_counts_mapping(i, 2) = length(unique(j)); %Number of EUs the barcode is linked to
    if mod(i, 10000) == 0
        disp(i)
    end
end

off_nanopore = zeros(size(reporter_ass_final, 1), 1);
oht_nanopore = zeros(size(reporter_ass_final, 1), 1);
gzv_nanopore = zeros(size(reporter_ass_final, 1), 1);
on_nanopore = zeros(size(reporter_ass_final, 1), 1);
read_char = zeros(size(reporter_ass_final, 1), 3);

for i = 1:length(off_nanopore)
    if mod(i, 1000) == 0
       disp(i)
    end
    x = bcall_off_exp(bcall_unique == BCall_final_final(i), :);
    y = bcall_oht_exp(bcall_unique == BCall_final_final(i), :);
    z = bcall_gzv_exp(bcall_unique == BCall_final_final(i), :);
    w = bcall_on_exp(bcall_unique == BCall_final_final(i), :);
    if ~isempty(x)
        read_char(i, 1) = bcall_counts(bcall_unique == BCall_final_final(i));
        read_char(i, 2) = 12;
        read_char(i, 3) = BCall_counts_mapping(un_nanoporeBCall == BCall_final_final(i), 2);
    else
        x = bc1_off_exp(bc1_unique == BC1_final_final(i), :);
        y = bc1_oht_exp(bc1_unique == BC1_final_final(i), :);
        z = bc1_gzv_exp(bc1_unique == BC1_final_final(i), :);
        w = bc1_on_exp(bc1_unique == BC1_final_final(i), :);
        if ~isempty(x)
            read_char(i, 1) = bc1_counts(bc1_unique == BC1_final_final(i));
            read_char(i, 2) = 1;
            read_char(i, 3) = BC1_counts_mapping(un_nanoporeBC1 == BC1_final_final(i), 2);
        else
            x = bc2_off_exp(bc2_unique == BC2_final_final(i), :);
            y = bc2_oht_exp(bc2_unique == BC2_final_final(i), :);
            z = bc2_gzv_exp(bc2_unique == BC2_final_final(i), :);
            w = bc2_on_exp(bc2_unique == BC2_final_final(i), :);
            if ~isempty(x)
                read_char(i, 1) = bc2_counts(bc2_unique == BC2_final_final(i));
                read_char(i, 2) = 2;
                read_char(i, 3) = BC2_counts_mapping(un_nanoporeBC2 == BC2_final_final(i), 2);
            else
                x = 0; y = 0; z = 0; w = 0;
            end
        end
    end
    off_nanopore(i) = x; oht_nanopore(i) = y; gzv_nanopore(i) = z; on_nanopore(i) = w;
end

%% Map illumina data to variants using the barcodes from illumina and nanopore indexing

ordered_eu_exp = zeros(1327104, 12);
eu_ordered_ass_split = zeros(1327104, 12);

for i = 1:6
    for j = 1:4
        for k = 1:2
            for l = 1:4
                for m = 1:3
                    for n = 1:4
                        for o = 1:2
                            for p = 1:4
                                for q = 1:3
                                    for r = 1:2
                                        for s = 1:4
                                            for t = 1:3
                                                eu_ordered_ass_split(t + (s-1)*3 + (r-1)*12 + (q-1)*24 + (p-1)*72 + (o-1)*288 + (n-1)*576 + (m-1)*2304 + (l-1)*6912 + (k-1)*27648 + (j-1)*55296 + (i-1)*221184, :) = [t s r q p o n m l k j i];
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

for i = 1:size(ordered_eu_exp, 1)
    if mod(i, 1000) == 0
        disp(i)
    end
    e = off_nanopore(eu_ass_final == i & read_char(:, 2) == 12);
    f = oht_nanopore(eu_ass_final == i & read_char(:, 2) == 12);
    g = gzv_nanopore(eu_ass_final == i & read_char(:, 2) == 12);
    h = on_nanopore(eu_ass_final == i & read_char(:, 2) == 12);
    a = e(~isnan(e)); b = f(~isnan(f)); c = g(~isnan(g)); d = h(~isnan(h));
    w = zeros(1, 4);
    if ~isempty(a)
        [f, xi] = ksdensity(log10(a)); fmaxa = 10.^xi(f == max(f)); fmaxa = mean(fmaxa);
        w(1) = length(a);
        a = mean(a);
    else
        a = 0; fmaxa = 0;
    end
    if ~isempty(b)
        [f, xi] = ksdensity(log10(b)); fmaxb = 10.^xi(f == max(f)); fmaxb = mean(fmaxb);
        w(2) = length(b);
        b = mean(b);
    else
        b = 0; fmaxb = 0;
    end
    if ~isempty(c)
        [f, xi] = ksdensity(log10(c)); fmaxc = 10.^xi(f == max(f)); fmaxc = mean(fmaxc);
        w(3) = length(c);
        c = mean(c);
    else
        c = 0; fmaxc = 0;
    end
    if ~isempty(d)
        [f, xi] = ksdensity(log10(d)); fmaxd = 10.^xi(f == max(f)); fmaxd = mean(fmaxd);
        w(4) = length(d);
        d = mean(d);
    else
        d = 0; fmaxd = 0;
    end
    ordered_eu_exp(i, :) = [a b c d fmaxa fmaxb fmaxc fmaxd w];
end
