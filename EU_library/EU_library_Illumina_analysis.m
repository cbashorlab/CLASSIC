clear
clc

%The illumina data contains bin/sample specific barcodes in reads 1 and 2,
%which we will extract, followed by extracting the variant specific
%barcodes and their counts in each of the sample bins

%Load sample specific barcodes, read in fastq data files
barcodes = ['AAGG'; 'AACC'; 'TTGG'; 'TTCC'; 'GGAA'; 'GGTT'; 'CCAA'; 'CCTT'; 'GAGA'; 'TCTC'];
[~, Sequence, ~] = fastqread('Sony-MA900_S1_L001_R1_001.fastq');
[~, Sequence2, ~] = fastqread('Sony-MA900_S1_L001_R2_001.fastq');

%Initialize variables to read sample specific barcodes
forward_ID = zeros(size(Sequence)); reverse_ID = forward_ID;

%Extract sample barcodes
for i = 1:length(Sequence)
    if mod(i, 100000) == 0
        disp(i)
    end
    a = cell2mat(Sequence(i));
    b = cell2mat(Sequence2(i));
    for j = 1:5%size(barcodes, 1)
        if contains(a(5:8), barcodes(j, :))
            reverse_ID(i) = j;
        end
        if contains(b(3:6), barcodes(j, :))
            forward_ID(i) = j;
        end
    end
    for j = 6:8
        if contains(a(5:8), barcodes(j, :))
            reverse_ID(i) = j;
        end
        if contains(b(4:7), barcodes(j, :))
            forward_ID(i) = j;
        end
    end
end

%Create a read count table for each sample barcode (QC)
primer_table = zeros(10, 10);

for i = 1:10
    for j = 1:10
        primer_table(i, j) = sum(forward_ID == j & reverse_ID == i);
    end
end

%Only retain barcodes that were included in this run
reverse_ID(reverse_ID == 6) = 3; reverse_ID(reverse_ID == 7) = 4;
fid_keep = forward_ID(forward_ID > 0 & reverse_ID > 0 & reverse_ID < 5);
rid_keep = reverse_ID(forward_ID > 0 & reverse_ID > 0 & reverse_ID < 5);
s_keep = Sequence(forward_ID > 0 & reverse_ID > 0 & reverse_ID < 5);
read_id = (rid_keep - 1)*8 + fid_keep;

%Read 1 from illumina runs always has better quality data, so we will use
%that to extract our variant barcodes. First, we reverse complement the
%reads so all barcodes are in the top strand and consistent with our
%nanopore indexing

for i = 1:length(s_keep)
    if mod(i, 100000) == 0 
        disp(i)
    end
    s_keep(i) = cellstr(seqrcomplement(cell2mat(s_keep(i))));
end

%Extract barcodes now
barcodes = cell(length(s_keep), 2);
barcode_lengths = zeros(size(s_keep));

for i = 1:length(s_keep)
    if mod(i, 100000) == 0
        disp(i)
    end
    a = cell2mat(s_keep(i));
    x = strfind(a, 'GAAACG');
    y = strfind(a, 'CAGTTC');

    if ~isempty(x) && ~isempty(y)
        barcodes(i, 1) = cellstr(a(x(1)+6:y(1)-1));
        barcode_lengths(i) = y(1) - x(1) - 6;
        if contains(a, 'GACATCAGC')
            barcodes(i, 2) = cellstr('1');
        elseif contains(a, 'TGCTTCTGC')
            barcodes(i, 2) = cellstr('2');
        else
            barcodes(i, 2) = cellstr('0');
        end
    else
        barcodes(i, 1) = cellstr('X');
        barcodes(i, 2) = cellstr('0');
    end
end

%Trim barcodes that have basecalling errors
barcodes_trimmed = barcodes(~contains(barcodes(:, 1), 'X') & barcode_lengths' == 18 & ~contains(barcodes(:, 2), '0'), :);
read_id_trimmed = read_id(~contains(barcodes(:, 1), 'X') & barcode_lengths' == 18 & ~contains(barcodes(:, 2), '0'));

bc1 = barcodes_trimmed(contains(barcodes_trimmed(:, 2), '1'), 1); bc2 = barcodes_trimmed(contains(barcodes_trimmed(:, 2), '2'), 1);
read_id_trimmed_bc1 = read_id_trimmed(contains(barcodes_trimmed(:, 2), '1'));
read_id_trimmed_bc2 = read_id_trimmed(contains(barcodes_trimmed(:, 2), '2'));

%Extract unique barcodes and compute their counts with accummarray (only
%the barcode 1 piece first)
[bc1_unique, ~, b] = unique(bc1);
bc1_dist = zeros(length(bc1_unique), 30);
bc1_counts = accumarray(b,1);

%Compute read counts for each barcode in each bin
for i = 1:30
    disp(i)
    x = accumarray(b(read_id_trimmed_bc1 == i), 1);
    bc1_dist(1:length(x), i) = x;
end

%Repeat the same steps as above, for barcode piece 2
[bc2_unique, ~, b] = unique(bc2);
bc2_dist = zeros(length(bc2_unique), 30);
bc2_counts = accumarray(b,1);

for i = 1:30
    disp(i)
    x = accumarray(b(read_id_trimmed_bc2 == i), 1);
    bc2_dist(1:length(x), i) = x;
end

%Trim barcodes that had < 10 reads
bc1_un_trim = bc1_unique(bc1_counts > 10); bc2_un_trim = bc2_unique(bc2_counts > 10);
bc1_dist_trim = bc1_dist(bc1_counts > 10, :); bc2_dist_trim = bc2_dist(bc2_counts > 10, :);
bc1_counts_trim = bc1_counts(bc1_counts > 10); bc2_counts_trim = bc2_counts(bc2_counts > 10);

%Normalize for read depth differences across different sample bins
bc1_dist_trim = 1000000*bc1_dist_trim./sum(bc1_dist_trim);
bc2_dist_trim = 1000000*bc2_dist_trim./sum(bc2_dist_trim);

%Normalize by fraction of library in each bin, and compute an average bin
%for each barcode (weighted average of counts in each bin)

%For library 1, sort 1 (4021), BC1 first
bc1_bin_averages_4021 = zeros(length(bc1_un_trim), 1);
bin1_percs = [1.47 9.15 30.75 18.34 9.31 7.25 5.88 5.03 5.89 1.32];

for i = 1:length(bc1_bin_averages_4021)
    disp(i)
    a = bc1_dist_trim(i, 1:10);
    if sum(a > 0) > 1
        bc1_bin_averages_4021(i) = sum(a.*bin1_percs.*(1:1:10))/sum(a.*bin1_percs);
    end
end

%Now library 1, sort 1 (4021), BC2
bc2_bin_averages_4021 = zeros(length(bc2_un_trim), 1);

for i = 1:length(bc2_bin_averages_4021)
    disp(i)
    a = bc2_dist_trim(i, 1:10);
    if sum(a > 0) > 1
        bc2_bin_averages_4021(i) = sum(a.*bin1_percs.*[1:1:10])/sum(a.*bin1_percs);
    end
end

%Now library 1, sort 2 (4022)
bc1_bin_averages_4022 = zeros(length(bc1_un_trim), 1);

for i = 1:length(bc1_bin_averages_4022)
    disp(i)
    a = bc1_dist_trim(i, 11:20);
    if sum(a > 0) > 1
        bc1_bin_averages_4022(i) = sum(a.*bin1_percs.*[1:1:10])/sum(a.*bin1_percs);
    end
end

bc2_bin_averages_4022 = zeros(length(bc2_un_trim), 1);

for i = 1:length(bc2_bin_averages_4022)
    disp(i)
    a = bc2_dist_trim(i, 11:20);
    if sum(a > 0) > 1
        bc2_bin_averages_4022(i) = sum(a.*bin1_percs.*[1:1:10])/sum(a.*bin1_percs);
    end
end

%Now library 2 (404)
bc1_bin_averages_404 = zeros(length(bc1_un_trim), 1);
bin1_percs_404 = [0.53 3.66 25.53 29.78 11 8.88 6.3 5.94 6.94 0.95];

for i = 1:length(bc1_bin_averages_404)
    disp(i)
    a = bc1_dist_trim(i, 21:30);
    if sum(a > 0) > 1
        bc1_bin_averages_404(i) = sum(a.*bin1_percs_404.*[1:1:10])/sum(a.*bin1_percs_404);
    end
end

bc2_bin_averages_404 = zeros(length(bc2_un_trim), 1);

for i = 1:length(bc2_bin_averages_404)
    disp(i)
    a = bc2_dist_trim(i, 21:30);
    if sum(a > 0) > 1
        bc2_bin_averages_404(i) = sum(a.*bin1_percs_404.*[1:1:10])/sum(a.*bin1_percs_404);
    end
end

%% Now we map illumina data to variants using the barcodes from illumina and nanopore indexing

%load('pKR402_nanopore_data.mat') 
%load('pKR402_nanopore_data.mat') %Variables provided separately, upon request

%We will match variants based on their nanopore barcodes (bc1_unique), and
%store all barcodes corresponding to each of the 384 variants

%Start with 402 Sort #1 (4021)
bc_384_binary = pKR402_assign_table > 0;
average_384_expression_4021 = cell(384, 1);

for i = 1:384
    disp(i)
    z = [];
    x = pKR402_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1); %If a barcode is uniquely mapped to the i-th element
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc1_un_trim, cell2mat(x(j)))) > 0
                y = bc1_bin_averages_4021(contains(bc1_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1); %Append its average expression to "z"
            end
        end
    end
    average_384_expression_4021(i) = {z(z > 0)}; %Assign all non 0 barcodes to the i-th element of the variable
end

%Now for variants that had the second BC2 piece
bc_384_binary = pKRsub_ass_table > 0;

for i = 1:384
    disp(i)
    if sum(pKRsub_ass_table(:, i)) > 100
    z = [];
    x = pKRsub_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1);
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc2_un_trim, cell2mat(x(j)))) > 0
                y = bc2_bin_averages_4021(contains(bc2_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1);
            end
        end
    end
        average_384_expression_4021(i) = {z(z > 0)};
    end
end

%Now do 402 sort 2 (4022)
bc_384_binary = pKR402_ass_table > 0;
average_384_expression_4022 = cell(384, 1);

for i = 1:384
    disp(i)
    z = [];
    x = pKR402_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1);
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc1_un_trim, cell2mat(x(j)))) > 0
                y = bc1_bin_averages_4022(contains(bc1_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1);
            end
        end
    end
    average_384_expression_4022(i) = {z(z > 0)};
end

bc_384_binary = pKRsub_ass_table > 0;

for i = 1:384
    disp(i)
    if sum(pKRsub_ass_table(:, i)) > 100
    z = [];
    x = pKRsub_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1);
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc2_un_trim, cell2mat(x(j)))) > 0
                y = bc2_bin_averages_4022(contains(bc2_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1);
            end
        end
    end
        average_384_expression_4022(i) = {z(z > 0)};
    end
end

%Now do library 2 (404)
bc_384_binary = pKR404_ass_table > 0;
average_384_expression_404 = cell(384, 1);

for i = 1:384
    disp(i)
    z = [];
    z2 = [];
    x = pKR404_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1);
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc1_un_trim, cell2mat(x(j)))) > 0
                y = bc1_bin_averages_404(contains(bc1_un_trim, cell2mat(x(j))));
                y2 = bc1_counts_trim(contains(bc1_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1);
                z2(end + 1) = y2(1);
            end
        end
    end
    average_384_expression_404(i) = {z(z > 0)};
    zeta(i) = {z2(z > 0)};
end

bc_384_binary = pKRsub_ass_table > 0;

for i = 1:384
    disp(i)
    if sum(pKRsub_ass_table(:, i)) > 100
    z = [];
    z2 = [];
    x = pKRsub_bc1_unique((sum(bc_384_binary' > 0) == 1)' & bc_384_binary(:, i) == 1);
    if ~isempty(x)
        for j = 1:length(x)
            if sum(contains(bc2_un_trim, cell2mat(x(j)))) > 0
                y = bc2_bin_averages_404(contains(bc2_un_trim, cell2mat(x(j))));
                y2 = bc2_counts_trim(contains(bc2_un_trim, cell2mat(x(j))));
                z(end + 1) = y(1);
                z2(end + 1) = y2(1);
            end
        end
    end
        average_384_expression_404(i) = {z(z > 0)};
        zeta(i) = {z2(z > 0)};
    end
end

%% Convert bins to expression

%This is done using a fit line obtained by plotting the mid points of all
%sorted bins against the bin number. The output from this function
%corresponds to the log10 expression level obtained from CLASSIC

bin_to_expression = @(x) (0.2291*x + 2.1578);

%The next few lines of code will, for each variant, extract their relevant
%barcodes, convert those to expression levels (log10), compute kernel
%density across the expression space, and save a single expression value 
%for each variant.

average_expression = zeros(384, 3);

for i = 1:size(average_expression, 1)
    disp(i)
    a = 10.^bin_to_expression(cell2mat(average_384_expression_4021(i)));
    b = 10.^bin_to_expression(cell2mat(average_384_expression_4022(i)));
    c = 10.^bin_to_expression(cell2mat(average_384_expression_404(i)));
    [f, xi] = ksdensity(a); [f2, xi2] = ksdensity(b); [f3, xi3] = ksdensity(c); 
    average_expression(i, :) = [(xi(f == max(f))) (xi2(f2 == max(f2))) (xi3(f3 == max(f3)))];
end

