%% Analyze DEC plaque purified

cd ('C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque purified mutants sequencing')

% Colors
black = [0 0 0];
red = [0.6350 0.0780 0.1840];
green = [0.4660 0.6740 0.1880];
blue = [0 0.4470 0.7410];
purple = [0.4940 0.1840 0.5560];
orange = [0.8500 0.3250 0.0980];

% Read SNPs
P0=vcfToTable('S0_variants_filtered.vcf');
% P0_2=vcfToTable('F0_variants_filtered.vcf');
% P0_3=vcfToTable('CS0_variants_filtered.vcf');
LB11=vcfToTable('LB11_variants_filtered.vcf');
LB12=vcfToTable('LB12_variants_filtered.vcf');
LB21=vcfToTable('LB21_variants_filtered.vcf');
%LB23=vcfToTable('LB23_variants_filtered.vcf'); %removed - looks like duplicated sequencing of LB31
LB31=vcfToTable('LB31_variants_filtered.vcf');
CS33=vcfToTable('CS33_variants_filtered.vcf');

% remove SNPs found in P0 from other

loc_P0=cell2mat(P0.Location);
loc_11=cell2mat(LB11.Location);
loc_12=cell2mat(LB12.Location);
loc_21=cell2mat(LB21.Location);
%loc_23=cell2mat(LB23.Location);
loc_31=cell2mat(LB31.Location);
loc_CS33=cell2mat(CS33.Location);

[~,remove]=ismember(loc_P0,loc_11); remove(remove==0)=[]; LB11(remove,:)=[];
[~,remove]=ismember(loc_P0,loc_12); remove(remove==0)=[]; LB12(remove,:)=[];
[~,remove]=ismember(loc_P0,loc_21); remove(remove==0)=[]; LB21(remove,:)=[];
%[~,remove]=ismember(loc_P0,loc_23); remove(remove==0)=[]; LB23(remove,:)=[];
[~,remove]=ismember(loc_P0,loc_31); remove(remove==0)=[]; LB31(remove,:)=[];
[~,remove]=ismember(loc_P0,loc_CS33); remove(remove==0)=[]; CS33(remove,:)=[];

%% Analyze on substitution level
% Types of SNPs
snpTypes = {'A>C', 'A>G', 'A>T', ...
            'C>A', 'C>G', 'C>T', ...
            'G>A', 'G>C', 'G>T', ...
            'T>A', 'T>C', 'T>G'};

% LB11
LB11_SNP_type=zeros(12,1);
for ii=1:size(LB11,1)
    if LB11.RefBase{ii}=='A' && LB11.AltBase{ii}=='C'
        LB11_SNP_type(1)=LB11_SNP_type(1)+1;
    elseif LB11.RefBase{ii}=='A' && LB11.AltBase{ii}=='G'
        LB11_SNP_type(2)=LB11_SNP_type(2)+1;
    elseif LB11.RefBase{ii}=='A' && LB11.AltBase{ii}=='T'
        LB11_SNP_type(3)=LB11_SNP_type(3)+1;
    elseif LB11.RefBase{ii}=='C' && LB11.AltBase{ii}=='A'
        LB11_SNP_type(4)=LB11_SNP_type(4)+1;
    elseif LB11.RefBase{ii}=='C' && LB11.AltBase{ii}=='G'
        LB11_SNP_type(5)=LB11_SNP_type(5)+1;
    elseif LB11.RefBase{ii}=='C' && LB11.AltBase{ii}=='T'
        LB11_SNP_type(6)=LB11_SNP_type(6)+1;    
    elseif LB11.RefBase{ii}=='G' && LB11.AltBase{ii}=='A'
        LB11_SNP_type(7)=LB11_SNP_type(7)+1;
    elseif LB11.RefBase{ii}=='G' && LB11.AltBase{ii}=='C'
        LB11_SNP_type(8)=LB11_SNP_type(8)+1; 
    elseif LB11.RefBase{ii}=='G' && LB11.AltBase{ii}=='T'
        LB11_SNP_type(9)=LB11_SNP_type(9)+1;
    elseif LB11.RefBase{ii}=='T' && LB11.AltBase{ii}=='A'
        LB11_SNP_type(10)=LB11_SNP_type(10)+1; 
    elseif LB11.RefBase{ii}=='T' && LB11.AltBase{ii}=='C'
        LB11_SNP_type(11)=LB11_SNP_type(11)+1;    
    elseif LB11.RefBase{ii}=='T' && LB11.AltBase{ii}=='G'
        LB11_SNP_type(12)=LB11_SNP_type(12)+1;
    end
end

% LB12
LB12_SNP_type=zeros(12,1);
for ii=1:size(LB12,1)
    if LB12.RefBase{ii}=='A' && LB12.AltBase{ii}=='C'
        LB12_SNP_type(1)=LB12_SNP_type(1)+1;
    elseif LB12.RefBase{ii}=='A' && LB12.AltBase{ii}=='G'
        LB12_SNP_type(2)=LB12_SNP_type(2)+1;
    elseif LB12.RefBase{ii}=='A' && LB12.AltBase{ii}=='T'
        LB12_SNP_type(3)=LB12_SNP_type(3)+1;
    elseif LB12.RefBase{ii}=='C' && LB12.AltBase{ii}=='A'
        LB12_SNP_type(4)=LB12_SNP_type(4)+1;
    elseif LB12.RefBase{ii}=='C' && LB12.AltBase{ii}=='G'
        LB12_SNP_type(5)=LB12_SNP_type(5)+1;
    elseif LB12.RefBase{ii}=='C' && LB12.AltBase{ii}=='T'
        LB12_SNP_type(6)=LB12_SNP_type(6)+1;    
    elseif LB12.RefBase{ii}=='G' && LB12.AltBase{ii}=='A'
        LB12_SNP_type(7)=LB12_SNP_type(7)+1;
    elseif LB12.RefBase{ii}=='G' && LB12.AltBase{ii}=='C'
        LB12_SNP_type(8)=LB12_SNP_type(8)+1; 
    elseif LB12.RefBase{ii}=='G' && LB12.AltBase{ii}=='T'
        LB12_SNP_type(9)=LB12_SNP_type(9)+1;
    elseif LB12.RefBase{ii}=='T' && LB12.AltBase{ii}=='A'
        LB12_SNP_type(10)=LB12_SNP_type(10)+1; 
    elseif LB12.RefBase{ii}=='T' && LB12.AltBase{ii}=='C'
        LB12_SNP_type(11)=LB12_SNP_type(11)+1;    
    elseif LB12.RefBase{ii}=='T' && LB12.AltBase{ii}=='G'
        LB12_SNP_type(12)=LB12_SNP_type(12)+1;
    end
end

% LB21
LB21_SNP_type=zeros(12,1);
for ii=1:size(LB21,1)
    if LB21.RefBase{ii}=='A' && LB21.AltBase{ii}=='C'
        LB21_SNP_type(1)=LB21_SNP_type(1)+1;
    elseif LB21.RefBase{ii}=='A' && LB21.AltBase{ii}=='G'
        LB21_SNP_type(2)=LB21_SNP_type(2)+1;
    elseif LB21.RefBase{ii}=='A' && LB21.AltBase{ii}=='T'
        LB21_SNP_type(3)=LB21_SNP_type(3)+1;
    elseif LB21.RefBase{ii}=='C' && LB21.AltBase{ii}=='A'
        LB21_SNP_type(4)=LB21_SNP_type(4)+1;
    elseif LB21.RefBase{ii}=='C' && LB21.AltBase{ii}=='G'
        LB21_SNP_type(5)=LB21_SNP_type(5)+1;
    elseif LB21.RefBase{ii}=='C' && LB21.AltBase{ii}=='T'
        LB21_SNP_type(6)=LB21_SNP_type(6)+1;    
    elseif LB21.RefBase{ii}=='G' && LB21.AltBase{ii}=='A'
        LB21_SNP_type(7)=LB21_SNP_type(7)+1;
    elseif LB21.RefBase{ii}=='G' && LB21.AltBase{ii}=='C'
        LB21_SNP_type(8)=LB21_SNP_type(8)+1; 
    elseif LB21.RefBase{ii}=='G' && LB21.AltBase{ii}=='T'
        LB21_SNP_type(9)=LB21_SNP_type(9)+1;
    elseif LB21.RefBase{ii}=='T' && LB21.AltBase{ii}=='A'
        LB21_SNP_type(10)=LB21_SNP_type(10)+1; 
    elseif LB21.RefBase{ii}=='T' && LB21.AltBase{ii}=='C'
        LB21_SNP_type(11)=LB21_SNP_type(11)+1;    
    elseif LB21.RefBase{ii}=='T' && LB21.AltBase{ii}=='G'
        LB21_SNP_type(12)=LB21_SNP_type(12)+1;
    end
end

% % LB23
% LB23_SNP_type=zeros(12,1);
% for ii=1:size(LB23,1)
%     if LB23.RefBase{ii}=='A' && LB23.AltBase{ii}=='C'
%         LB23_SNP_type(1)=LB23_SNP_type(1)+1;
%     elseif LB23.RefBase{ii}=='A' && LB23.AltBase{ii}=='G'
%         LB23_SNP_type(2)=LB23_SNP_type(2)+1;
%     elseif LB23.RefBase{ii}=='A' && LB23.AltBase{ii}=='T'
%         LB23_SNP_type(3)=LB23_SNP_type(3)+1;
%     elseif LB23.RefBase{ii}=='C' && LB23.AltBase{ii}=='A'
%         LB23_SNP_type(4)=LB23_SNP_type(4)+1;
%     elseif LB23.RefBase{ii}=='C' && LB23.AltBase{ii}=='G'
%         LB23_SNP_type(5)=LB23_SNP_type(5)+1;
%     elseif LB23.RefBase{ii}=='C' && LB23.AltBase{ii}=='T'
%         LB23_SNP_type(6)=LB23_SNP_type(6)+1;    
%     elseif LB23.RefBase{ii}=='G' && LB23.AltBase{ii}=='A'
%         LB23_SNP_type(7)=LB23_SNP_type(7)+1;
%     elseif LB23.RefBase{ii}=='G' && LB23.AltBase{ii}=='C'
%         LB23_SNP_type(8)=LB23_SNP_type(8)+1; 
%     elseif LB23.RefBase{ii}=='G' && LB23.AltBase{ii}=='T'
%         LB23_SNP_type(9)=LB23_SNP_type(9)+1;
%     elseif LB23.RefBase{ii}=='T' && LB23.AltBase{ii}=='A'
%         LB23_SNP_type(10)=LB23_SNP_type(10)+1; 
%     elseif LB23.RefBase{ii}=='T' && LB23.AltBase{ii}=='C'
%         LB23_SNP_type(11)=LB23_SNP_type(11)+1;    
%     elseif LB23.RefBase{ii}=='T' && LB23.AltBase{ii}=='G'
%         LB23_SNP_type(12)=LB23_SNP_type(12)+1;
%     end
% end

% LB31
LB31_SNP_type=zeros(12,1);
for ii=1:size(LB31,1)
    if LB31.RefBase{ii}=='A' && LB31.AltBase{ii}=='C'
        LB31_SNP_type(1)=LB31_SNP_type(1)+1;
    elseif LB31.RefBase{ii}=='A' && LB31.AltBase{ii}=='G'
        LB31_SNP_type(2)=LB31_SNP_type(2)+1;
    elseif LB31.RefBase{ii}=='A' && LB31.AltBase{ii}=='T'
        LB31_SNP_type(3)=LB31_SNP_type(3)+1;
    elseif LB31.RefBase{ii}=='C' && LB31.AltBase{ii}=='A'
        LB31_SNP_type(4)=LB31_SNP_type(4)+1;
    elseif LB31.RefBase{ii}=='C' && LB31.AltBase{ii}=='G'
        LB31_SNP_type(5)=LB31_SNP_type(5)+1;
    elseif LB31.RefBase{ii}=='C' && LB31.AltBase{ii}=='T'
        LB31_SNP_type(6)=LB31_SNP_type(6)+1;    
    elseif LB31.RefBase{ii}=='G' && LB31.AltBase{ii}=='A'
        LB31_SNP_type(7)=LB31_SNP_type(7)+1;
    elseif LB31.RefBase{ii}=='G' && LB31.AltBase{ii}=='C'
        LB31_SNP_type(8)=LB31_SNP_type(8)+1; 
    elseif LB31.RefBase{ii}=='G' && LB31.AltBase{ii}=='T'
        LB31_SNP_type(9)=LB31_SNP_type(9)+1;
    elseif LB31.RefBase{ii}=='T' && LB31.AltBase{ii}=='A'
        LB31_SNP_type(10)=LB31_SNP_type(10)+1; 
    elseif LB31.RefBase{ii}=='T' && LB31.AltBase{ii}=='C'
        LB31_SNP_type(11)=LB31_SNP_type(11)+1;    
    elseif LB31.RefBase{ii}=='T' && LB31.AltBase{ii}=='G'
        LB31_SNP_type(12)=LB31_SNP_type(12)+1;
    end
end

% CS33
CS33_SNP_type=zeros(12,1);
for ii=1:size(CS33,1)
    if CS33.RefBase{ii}=='A' && CS33.AltBase{ii}=='C'
        CS33_SNP_type(1)=CS33_SNP_type(1)+1;
    elseif CS33.RefBase{ii}=='A' && CS33.AltBase{ii}=='G'
        CS33_SNP_type(2)=CS33_SNP_type(2)+1;
    elseif CS33.RefBase{ii}=='A' && CS33.AltBase{ii}=='T'
        CS33_SNP_type(3)=CS33_SNP_type(3)+1;
    elseif CS33.RefBase{ii}=='C' && CS33.AltBase{ii}=='A'
        CS33_SNP_type(4)=CS33_SNP_type(4)+1;
    elseif CS33.RefBase{ii}=='C' && CS33.AltBase{ii}=='G'
        CS33_SNP_type(5)=CS33_SNP_type(5)+1;
    elseif CS33.RefBase{ii}=='C' && CS33.AltBase{ii}=='T'
        CS33_SNP_type(6)=CS33_SNP_type(6)+1;    
    elseif CS33.RefBase{ii}=='G' && CS33.AltBase{ii}=='A'
        CS33_SNP_type(7)=CS33_SNP_type(7)+1;
    elseif CS33.RefBase{ii}=='G' && CS33.AltBase{ii}=='C'
        CS33_SNP_type(8)=CS33_SNP_type(8)+1; 
    elseif CS33.RefBase{ii}=='G' && CS33.AltBase{ii}=='T'
        CS33_SNP_type(9)=CS33_SNP_type(9)+1;
    elseif CS33.RefBase{ii}=='T' && CS33.AltBase{ii}=='A'
        CS33_SNP_type(10)=CS33_SNP_type(10)+1; 
    elseif CS33.RefBase{ii}=='T' && CS33.AltBase{ii}=='C'
        CS33_SNP_type(11)=CS33_SNP_type(11)+1;    
    elseif CS33.RefBase{ii}=='T' && CS33.AltBase{ii}=='G'
        CS33_SNP_type(12)=CS33_SNP_type(12)+1;
    end
end

f=figure; hold;
barh(1.5:12.5,LB11_SNP_type,'FaceColor',red,'BarWidth',0.3)
barh(1.2:12.2,LB12_SNP_type,'FaceColor',green,'BarWidth',0.3)
barh(0.9:11.9,LB21_SNP_type,'FaceColor',blue,'BarWidth',0.3)
barh(0.6:11.6,LB31_SNP_type,'FaceColor',purple,'BarWidth',0.3)
yticks(1:12)
yticklabels(snpTypes)
xlabel 'Number of SNPs'
ylim ([0.2,12.8])
f.Position=[680   499   481   479];
set(gca,'FontSize',16)
box('off')
legend({'D1.1','D1.2','D2.1','D3.1'},'Location','southeast','Box','off','FontSize',14)
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_substitution_level','emf')

%% VAF histogram

VAF_11=LB11.VAF;
VAF_12=LB12.VAF;
VAF_21=LB21.VAF;
%VAF_23=LB23.VAF;
VAF_31=LB31.VAF;

f=figure; hold;
subplot(2,2,1); histogram(VAF_11,'BinEdges',0:0.1:1,'FaceColor',red); title 'D1.1'
set(gca,'FontSize',16); ylabel '# SNPs'; xlabel 'SNP Frequency'; yticks(0:25:100)
subplot(2,2,2); histogram(VAF_12,'BinEdges',0:0.1:1,'FaceColor',green); title 'D1.2'
set(gca,'FontSize',16); yticks(0:20:60) 
subplot(2,2,3); histogram(VAF_21,'BinEdges',0:0.1:1,'FaceColor',blue); title 'D2.1'
set(gca,'FontSize',16); yticks(0:25:100)
subplot(2,2,4); histogram(VAF_31,'BinEdges',0:0.1:1,'FaceColor',purple); title 'D3.1'
set(gca,'FontSize',16); yticks(0:20:80)
f.Position=[680   499   481   479];
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_VAF_histograms','emf')

%% SNPs / VAF / Location

loc_11=cell2mat(LB11.Location);
loc_12=cell2mat(LB12.Location);
loc_21=cell2mat(LB21.Location);
%loc_23=cell2mat(LB23.Location);
loc_31=cell2mat(LB31.Location);

VAF_11=LB11.VAF;
VAF_12=LB12.VAF;
VAF_21=LB21.VAF;
%VAF_23=LB23.VAF;
VAF_31=LB31.VAF;

f=figure; hold;
scatter(loc_11,VAF_11,80,red,"filled",'MarkerFaceAlpha',0.5)
scatter(loc_12,VAF_12,80,green,"filled",'MarkerFaceAlpha',0.5)
scatter(loc_21,VAF_21,80,blue,"filled",'MarkerFaceAlpha',0.5)
scatter(loc_31,VAF_31,80,purple,"filled",'MarkerFaceAlpha',0.5)
xlim([1,136376])
xticks([]);
ylabel 'SNP Frequency'
yticks(0:0.25:1)
set(gca,'FontSize',20)
box('off')
f.Position=[275,676,1375,302];
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_location_vs_VAF','emf')

%% gene level analysis

%read in HSV-1 genes location from GTF file
geneTable = extractGeneDataFromExcel('C:\Users\Nir\UC Irvine\BioSci-Drayman Lab-Team - Documents\NGS sequencing\Laura\Vervet_ND02_genome\S17_exon.gtf');

%count SNPs in the viral genes
snpCountsTable = MapSnpToGene(geneTable, LB11,LB12,LB21,LB31);

% remove the miR
snpCountsTable(78:93,:)=[];

%find genes hit in all 4 mutants
snps_per_gene=snpCountsTable{:,4:end};
hits=[];
for ii=1:size(snps_per_gene,1)
    if ~ismember(0,snps_per_gene(ii,:))
        hits=[hits;ii];
    end
end

%average # SNPs
f=figure;
bar(1:length(hits),sum(snpCountsTable{hits,4:end},2)./4,'FaceColor',black)
xticks(1:length(hits))
xticklabels(snpCountsTable{hits,1})
%set(gca,'XTickLabelRotation',30)
ylabel 'Average # SNPs'
f.Position=[680   499   481   479];
set(gca,'FontSize',16)
box('off')
title 'Genes mutated in all isolates'
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_gene_level','emf')

% stacked bar # SNPs
f=figure;
b=bar(1:length(hits),snpCountsTable{hits,4:end},'stacked');
b(1).FaceColor=red;
b(2).FaceColor=green;
b(3).FaceColor=blue;
b(4).FaceColor=purple;
xticks(1:length(hits))
xticklabels(snpCountsTable{hits,1})
set(gca,'XTickLabelRotation',60)
ylabel '# SNPs'
f.Position=[680   499   481   479];
set(gca,'FontSize',16)
legend({'D1.1','D1.2','D2.1','D3.1'},'Location','northwest','Box','off','FontSize',16)
box('off')
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_gene_level_stacked','emf')


% SNPs locations in the hits
f=figure; hold;
for ii=1:14
    subplot(3,5,ii)
    scatter(loc_11,VAF_11,80,red,"filled",'MarkerFaceAlpha',0.5); hold;
    scatter(loc_12,VAF_12,80,green,"filled",'MarkerFaceAlpha',0.5)
    scatter(loc_21,VAF_21,80,blue,"filled",'MarkerFaceAlpha',0.5)
    %scatter(loc_23,VAF_23,80,purple,"filled",'MarkerFaceAlpha',0.5)
    scatter(loc_31,VAF_31,80,purple,"filled",'MarkerFaceAlpha',0.5)
    xlim([geneTable{hits(ii),2},geneTable{hits(ii),3}])
    ylim([0,1])
    title(geneTable{hits(ii),1})
end
%saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\plaque_purified_gene_level','emf')

