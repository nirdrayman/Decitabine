function snpCountsTable = MapSnpToGene(geneTable, varargin)
    % Initialize array to store SNP counts for each gene
    numGenes = height(geneTable);
    snpCounts = zeros(numGenes, nargin - 1); % Each column for one SNP file

    % Loop through each SNP file
    for k = 1:nargin - 1
        SNP_file = varargin{k};
        
        locations = cell2mat(SNP_file.Location);
        
        % Loop through each SNP in the current file
        for i = 1:length(locations)
            snpLocation = locations(i);
            
            % Find genes that the SNP maps to
            for j = 1:numGenes
                geneStart = geneTable.gene_start(j);
                geneEnd = geneTable.gene_end(j);
                if snpLocation >= geneStart && snpLocation <= geneEnd
                    snpCounts(j, k) = snpCounts(j, k) + 1;
                end
            end
        end
    end

    % Create the output table with SNP counts for each file
    snpCountsTable = geneTable;
    callingWorkspaceVars = arrayfun(@inputname, 2:nargin, 'UniformOutput', false);
    for k = 1:length(varargin)
        snpCountsTable.(callingWorkspaceVars{k}) = snpCounts(:, k);
    end
end
