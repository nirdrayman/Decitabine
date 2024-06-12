function snpTable = vcfToTable(vcfFilePath) %outputTablePath)
    % Initialize cell arrays to store the extracted information
    locations = {};
    refBases = {};
    altBases = {};
    depths = [];
    vafs = [];

    % Open the VCF file
    fid = fopen(vcfFilePath, 'r');
    if fid == -1
        error('Unable to open VCF file %s.', vcfFilePath);
    end

    % Read the file line by line
    while ~feof(fid)
        line = fgetl(fid);
        if startsWith(line, '#')
            continue;
        end
        
        fields = strsplit(line, '\t');
        chrom = fields{1};
        pos = fields{2};
        ref = fields{4};
        alt = fields{5};
        info = fields{8};

        % Only consider SNPs (both REF and ALT are single nucleotides)
        if length(ref) == 1 && length(alt) == 1
            % Convert position to numeric value
            location = str2double(pos);

            % Parse INFO field
            infoFields = strsplit(info, ';');
            infoMap = containers.Map();
            for i = 1:length(infoFields)
                keyValue = strsplit(infoFields{i}, '=');
                if length(keyValue) == 2
                    key = keyValue{1};
                    value = keyValue{2};
                    infoMap(key) = value;
                end
            end

            % Extract relevant fields
            if isKey(infoMap, 'DP')
                depth = str2double(infoMap('DP'));
            else
                depth = NaN;
            end

            if isKey(infoMap, 'DP4')
                dp4 = str2double(strsplit(infoMap('DP4'), ','));
                refCount = dp4(1) + dp4(2);
                altCount = dp4(3) + dp4(4);
                vaf = altCount / (refCount + altCount);
            else
                vaf = NaN;
            end

            % Append to arrays
            locations{end+1, 1} = location;
            refBases{end+1, 1} = ref;
            altBases{end+1, 1} = alt;
            depths(end+1, 1) = depth;
            vafs(end+1, 1) = vaf;
        end
    end

    % Close the VCF file
    fclose(fid);

    % Create table
    snpTable = table(locations, refBases, altBases, depths, vafs, ...
        'VariableNames', {'Location', 'RefBase', 'AltBase', 'DepthOfSequencing', 'VAF'});

%     % Save table to CSV file if output path is given
%     if nargin > 1
%         writetable(snpTable, outputTablePath);
%         fprintf('Data has been successfully exported to %s\n', outputTablePath);
%     end
end

% Example usage:
% snpTable = vcfToTable('input.vcf', 'output_table.csv');