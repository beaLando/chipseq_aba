# This generates files that can be used instead of files generated at step 5a
# In this case I downsample all samples to a minimum value, instead of downsampling one of two samples in pair of replicates to equalize read depth between replicates
# This step is to do a quality check and see if similar peaks are identified at low depths, given we have samples with very low depths

conda activate chippy

in_dir="mapped"
out_dir="mapped/downsampled_ramp"
seed=42
mkdir -p ${out_dir}

# Downsample at three different depths (check which are the three shallowest samples and take those three depths
## Sample matrix
dos2unix ./qc_out/chipqc_all/chip_readsize_fragsize.csv
my_array_samp=( $(cut -d ";" -f1 ./qc_out/chipqc_all/chip_readsize_fragsize.csv | sort -u) ) #I should have actually removed inputs from here, the non-downsampled inputs in ./mapped can be used

## The three samples with the least reads were 8, 6, and 4 for me, in increasing order for number of reads mapped
dpth1=`samtools view -c ./mapped/Unknown_CE349-003R0008.sorted.bam`
dpth2=`samtools view -c ./mapped/Unknown_CE349-003R0006.sorted.bam`
dpth3=`samtools view -c ./mapped/Unknown_CE349-003R0004.sorted.bam`

my_array_dpth=( ${dpth1} ${dpth2} ${dpth3} )

for dpth in "${my_array_dpth[@]}"; do
    echo "subsampling ${dpth} random reads"

    out_dir_reads=$(echo ${out_dir}"/"${dpth}"_reads")
    mkdir -p ${out_dir_reads}

    for samp in "${my_array_samp[@]}"; do
        echo "in genotype ${samp}"

    biostar145820 --seed ${seed} -n ${dpth} -o ${out_dir_reads}/${samp}.dwnsmp.bam ${in_dir}/${samp}.sorted.bam

	samtools sort -o ${out_dir_reads}/${samp}.dwnsmp.sorted.bam ${out_dir_reads}/${samp}.dwnsmp.bam && rm ${out_dir_reads}/${samp}.dwnsmp.bam

	samtools index ${out_dir_reads}/${samp}.dwnsmp.sorted.bam

    done
done


## Also subsample to have in total between the two cocit samples the same number of mapped reads as the number of mapped reads for the two cocitaba samples together
dpth4=2677917 #(dpth1+dpth2)/2
samp1="Unknown_CE349-003R0002"
samp2="Unknown_CE349-003R0004"
my_array_samp=( ${samp1} ${samp2} )

for samp in "${my_array_samp[@]}"; do
        echo "in genotype ${samp}"

    out_dir_reads=$(echo ${out_dir}"/"${dpth4}"_reads")
    mkdir -p ${out_dir_reads}

    biostar145820 --seed ${seed} -n ${dpth4} -o ${out_dir_reads}/${samp}.dwnsmp.bam ${in_dir}/${samp}.sorted.bam

	samtools sort -o ${out_dir_reads}/${samp}.dwnsmp.sorted.bam ${out_dir_reads}/${samp}.dwnsmp.bam && rm ${out_dir_reads}/${samp}.dwnsmp.bam

	samtools index ${out_dir_reads}/${samp}.dwnsmp.sorted.bam
done


conda deactivate

rm mapped/downsampled_ramp/*/.sorted.bam

