%% Analyze P145-P148
% download VCF files from GitHub

% P145 - parental virus. P146-P148: Decitabine resistant HSV-1

% Colors
black = [0 0 0];
red = [0.6350 0.0780 0.1840];
green = [0.4660 0.6740 0.1880];
blue = [0 0.4470 0.7410];
purple = [0.4940 0.1840 0.5560];
orange = [0.8500 0.3250 0.0980];

% Read SNPs
P145=vcfToTable('P145_variants_filtered.vcf');
P146=vcfToTable('P146_variants_filtered.vcf');
P147=vcfToTable('P147_variants_filtered.vcf');
P148=vcfToTable('P148_variants_filtered.vcf');

% remove SNPs found in P145 from P146-P148

loc_P145=cell2mat(P145.Location);
loc_P146=cell2mat(P146.Location);
loc_P147=cell2mat(P147.Location);
loc_P148=cell2mat(P148.Location);

[~,remove]=ismember(loc_P145,loc_P146);
remove(remove==0)=[];
P146(remove,:)=[];

[~,remove]=ismember(loc_P145,loc_P147);
remove(remove==0)=[];
P147(remove,:)=[];

[~,remove]=ismember(loc_P145,loc_P148);
remove(remove==0)=[];
P148(remove,:)=[];

%% Analyze on substitution level
% Types of SNPs
snpTypes = {'A>C', 'A>G', 'A>T', ...
            'C>A', 'C>G', 'C>T', ...
            'G>A', 'G>C', 'G>T', ...
            'T>A', 'T>C', 'T>G'};

% P145
P145_SNP_type=zeros(12,1);
for ii=1:size(P145,1)
    if P145.RefBase{ii}=='A' && P145.AltBase{ii}=='C'
        P145_SNP_type(1)=P145_SNP_type(1)+1;
    elseif P145.RefBase{ii}=='A' && P145.AltBase{ii}=='G'
        P145_SNP_type(2)=P145_SNP_type(2)+1;
    elseif P145.RefBase{ii}=='A' && P145.AltBase{ii}=='T'
        P145_SNP_type(3)=P145_SNP_type(3)+1;
    elseif P145.RefBase{ii}=='C' && P145.AltBase{ii}=='A'
        P145_SNP_type(4)=P145_SNP_type(4)+1;
    elseif P145.RefBase{ii}=='C' && P145.AltBase{ii}=='G'
        P145_SNP_type(5)=P145_SNP_type(5)+1;
    elseif P145.RefBase{ii}=='C' && P145.AltBase{ii}=='T'
        P145_SNP_type(6)=P145_SNP_type(6)+1;    
    elseif P145.RefBase{ii}=='G' && P145.AltBase{ii}=='A'
        P145_SNP_type(7)=P145_SNP_type(7)+1;
    elseif P145.RefBase{ii}=='G' && P145.AltBase{ii}=='C'
        P145_SNP_type(8)=P145_SNP_type(8)+1; 
    elseif P145.RefBase{ii}=='G' && P145.AltBase{ii}=='T'
        P145_SNP_type(9)=P145_SNP_type(9)+1;
    elseif P145.RefBase{ii}=='T' && P145.AltBase{ii}=='A'
        P145_SNP_type(10)=P145_SNP_type(10)+1; 
    elseif P145.RefBase{ii}=='T' && P145.AltBase{ii}=='C'
        P145_SNP_type(11)=P145_SNP_type(11)+1;    
    elseif P145.RefBase{ii}=='T' && P145.AltBase{ii}=='G'
        P145_SNP_type(12)=P145_SNP_type(12)+1;
    end
end

% P146
P146_SNP_type=zeros(12,1);
for ii=1:size(P146,1)
    if P146.RefBase{ii}=='A' && P146.AltBase{ii}=='C'
        P146_SNP_type(1)=P146_SNP_type(1)+1;
    elseif P146.RefBase{ii}=='A' && P146.AltBase{ii}=='G'
        P146_SNP_type(2)=P146_SNP_type(2)+1;
    elseif P146.RefBase{ii}=='A' && P146.AltBase{ii}=='T'
        P146_SNP_type(3)=P146_SNP_type(3)+1;
    elseif P146.RefBase{ii}=='C' && P146.AltBase{ii}=='A'
        P146_SNP_type(4)=P146_SNP_type(4)+1;
    elseif P146.RefBase{ii}=='C' && P146.AltBase{ii}=='G'
        P146_SNP_type(5)=P146_SNP_type(5)+1;
    elseif P146.RefBase{ii}=='C' && P146.AltBase{ii}=='T'
        P146_SNP_type(6)=P146_SNP_type(6)+1;    
    elseif P146.RefBase{ii}=='G' && P146.AltBase{ii}=='A'
        P146_SNP_type(7)=P146_SNP_type(7)+1;
    elseif P146.RefBase{ii}=='G' && P146.AltBase{ii}=='C'
        P146_SNP_type(8)=P146_SNP_type(8)+1; 
    elseif P146.RefBase{ii}=='G' && P146.AltBase{ii}=='T'
        P146_SNP_type(9)=P146_SNP_type(9)+1;
    elseif P146.RefBase{ii}=='T' && P146.AltBase{ii}=='A'
        P146_SNP_type(10)=P146_SNP_type(10)+1; 
    elseif P146.RefBase{ii}=='T' && P146.AltBase{ii}=='C'
        P146_SNP_type(11)=P146_SNP_type(11)+1;    
    elseif P146.RefBase{ii}=='T' && P146.AltBase{ii}=='G'
        P146_SNP_type(12)=P146_SNP_type(12)+1;
    end
end

% P147
P147_SNP_type=zeros(12,1);
for ii=1:size(P147,1)
    if P147.RefBase{ii}=='A' && P147.AltBase{ii}=='C'
        P147_SNP_type(1)=P147_SNP_type(1)+1;
    elseif P147.RefBase{ii}=='A' && P147.AltBase{ii}=='G'
        P147_SNP_type(2)=P147_SNP_type(2)+1;
    elseif P147.RefBase{ii}=='A' && P147.AltBase{ii}=='T'
        P147_SNP_type(3)=P147_SNP_type(3)+1;
    elseif P147.RefBase{ii}=='C' && P147.AltBase{ii}=='A'
        P147_SNP_type(4)=P147_SNP_type(4)+1;
    elseif P147.RefBase{ii}=='C' && P147.AltBase{ii}=='G'
        P147_SNP_type(5)=P147_SNP_type(5)+1;
    elseif P147.RefBase{ii}=='C' && P147.AltBase{ii}=='T'
        P147_SNP_type(6)=P147_SNP_type(6)+1;    
    elseif P147.RefBase{ii}=='G' && P147.AltBase{ii}=='A'
        P147_SNP_type(7)=P147_SNP_type(7)+1;
    elseif P147.RefBase{ii}=='G' && P147.AltBase{ii}=='C'
        P147_SNP_type(8)=P147_SNP_type(8)+1; 
    elseif P147.RefBase{ii}=='G' && P147.AltBase{ii}=='T'
        P147_SNP_type(9)=P147_SNP_type(9)+1;
    elseif P147.RefBase{ii}=='T' && P147.AltBase{ii}=='A'
        P147_SNP_type(10)=P147_SNP_type(10)+1; 
    elseif P147.RefBase{ii}=='T' && P147.AltBase{ii}=='C'
        P147_SNP_type(11)=P147_SNP_type(11)+1;    
    elseif P147.RefBase{ii}=='T' && P147.AltBase{ii}=='G'
        P147_SNP_type(12)=P147_SNP_type(12)+1;
    end
end

% P148
P148_SNP_type=zeros(12,1);
for ii=1:size(P148,1)
    if P148.RefBase{ii}=='A' && P148.AltBase{ii}=='C'
        P148_SNP_type(1)=P148_SNP_type(1)+1;
    elseif P148.RefBase{ii}=='A' && P148.AltBase{ii}=='G'
        P148_SNP_type(2)=P148_SNP_type(2)+1;
    elseif P148.RefBase{ii}=='A' && P148.AltBase{ii}=='T'
        P148_SNP_type(3)=P148_SNP_type(3)+1;
    elseif P148.RefBase{ii}=='C' && P148.AltBase{ii}=='A'
        P148_SNP_type(4)=P148_SNP_type(4)+1;
    elseif P148.RefBase{ii}=='C' && P148.AltBase{ii}=='G'
        P148_SNP_type(5)=P148_SNP_type(5)+1;
    elseif P148.RefBase{ii}=='C' && P148.AltBase{ii}=='T'
        P148_SNP_type(6)=P148_SNP_type(6)+1;    
    elseif P148.RefBase{ii}=='G' && P148.AltBase{ii}=='A'
        P148_SNP_type(7)=P148_SNP_type(7)+1;
    elseif P148.RefBase{ii}=='G' && P148.AltBase{ii}=='C'
        P148_SNP_type(8)=P148_SNP_type(8)+1; 
    elseif P148.RefBase{ii}=='G' && P148.AltBase{ii}=='T'
        P148_SNP_type(9)=P148_SNP_type(9)+1;
    elseif P148.RefBase{ii}=='T' && P148.AltBase{ii}=='A'
        P148_SNP_type(10)=P148_SNP_type(10)+1; 
    elseif P148.RefBase{ii}=='T' && P148.AltBase{ii}=='C'
        P148_SNP_type(11)=P148_SNP_type(11)+1;    
    elseif P148.RefBase{ii}=='T' && P148.AltBase{ii}=='G'
        P148_SNP_type(12)=P148_SNP_type(12)+1;
    end
end

f=figure; hold;
barh(1.3:12.3,P146_SNP_type,'FaceColor',red,'BarWidth',0.3)
barh(1:12,P147_SNP_type,'FaceColor',blue,'BarWidth',0.3)
barh(0.7:11.7,P148_SNP_type,'FaceColor',green,'BarWidth',0.3)
yticks(1:12)
yticklabels(snpTypes)
xlabel 'Number of SNPs'
ylim ([0.2,12.8])
f.Position=[680   499   481   479];
set(gca,'FontSize',16)
box('off')
legend({'DEC-1','DEC-2','DEC-3'},'Location','northeast','Box','off','FontSize',14)
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\P146_P148_substitution_level','emf')

% percent G>C / C>G
transversions_P146=(P146_SNP_type(5)+P146_SNP_type(8))./sum(P146_SNP_type)*100
transversions_P147=(P147_SNP_type(5)+P148_SNP_type(8))./sum(P147_SNP_type)*100
transversions_P148=(P148_SNP_type(5)+P148_SNP_type(8))./sum(P148_SNP_type)*100

%% SNPs / VAF / Location

loc_P146=cell2mat(P146.Location);
loc_P147=cell2mat(P147.Location);
loc_P148=cell2mat(P148.Location);

VAF_P146=P146.VAF;
VAF_P147=P147.VAF;
VAF_P148=P148.VAF;

f=figure; hold;
scatter(loc_P146,VAF_P146,80,red,"filled",'MarkerFaceAlpha',0.7)
scatter(loc_P147,VAF_P147,80,blue,"filled",'MarkerFaceAlpha',0.7)
scatter(loc_P148,VAF_P148,80,green,"filled",'MarkerFaceAlpha',0.7)
xlim([1,136376])
xticks([]);
ylabel 'Variant Frequency'
yticks(0:0.25:1)
set(gca,'FontSize',20)
box('off')
f.Position=[275,676,1375,302];
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\P146_P148_location_vs_VAF','emf')

%% Number of SNPs
f=figure; hold;
bar(0.6,length(loc_P146),'FaceColor',red,'BarWidth',0.3)
bar(1,length(loc_P147),'FaceColor',blue,'BarWidth',0.3)
bar(1.4,length(loc_P148),'FaceColor',green,'BarWidth',0.3)
xticks(0.6:0.4:1.4)
xticklabels({'DEC-1','DEC-2','DEC-3'})
ylabel 'Number of SNPs'
yticks(0:50:200)
xlim ([0.3,1.7])
f.Position=[680   499   240   240];
set(gca,'FontSize',14)
box('off')
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\Num_SNPs_pool','emf')

%% Concordance
common_P146_P147=intersect(loc_P146,loc_P147);
common_P147_P148=intersect(loc_P147,loc_P148);
common_P146_P148=intersect(loc_P146,loc_P148);

f=figure; hold;
bar(0.6,length(common_P146_P147),'FaceColor',black,'BarWidth',0.3)
%bar(1,length(common_P146_P148),'FaceColor',black,'BarWidth',0.3)
bar(1,0,'FaceColor',black,'BarWidth',0.3)
bar(1.4,length(common_P147_P148),'FaceColor',black,'BarWidth',0.3)
xticks(0.6:0.4:1.4)
xticklabels({'DEC-1&2','DEC-1&3','DEC-2&3'})
ylabel 'Shared SNPs'
xlim ([0.3,1.7])
f.Position=[680   499   240   240];
set(gca,'FontSize',14)
box('off')
saveas(f,'C:\Users\Nir\OneDrive - UC Irvine\Writing\Manuscripts\24 HSV1 drug repurposing screen\Num_sahred_SNPs_P146-P148','emf')

