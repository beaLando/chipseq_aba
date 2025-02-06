# GFF formatting
## First need to create a matrix with scores based on the read density values in the bigWig files for each gene/genomic region (in BED format!!):
### Going from GFF to BED is surprisingly difficult, because GFF files are nasty
#### First, I fix original GFF file, then I convert it to GTF, using agat:

conda activate gff_gtf

agat_convert_sp_gxf2gxf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff
agat_convert_sp_gff2gtf.pl --gff ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff -o ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gtf

#### I also save bed file with only +/-1000bp around TSS genes

awk 'BEGIN{OFS="\t"}{if($3=="gene"){$3="tss";if($7 == "+"){$5=$4}else{$4=$5};print $0}}' ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat.gff | awk {'print $1"\t"$4"\t"$5"\t"$9"\t"$6"\t"$7'} | bedops --range -1000:1000 --everything - > ara_refs/GCF_000001735.4_TAIR10.1_genomic_agat_TSS.bed

conda deactivate

