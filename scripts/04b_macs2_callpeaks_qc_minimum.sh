# Other subsamplings

# Sum of mapped reads of rep1+rep2 for cocit samples equals sum of mapped reads of rep1+rep2 of cocitaba samples, with reads distributed evenly or unevenly between replicates


## CASE1: subsampling distributes reads equally to subsampled rep1 and rep2 of cocit samples
mkdir -p ./macs2_out/2677917_reads
mkdir -p ./macs2_out/2677917_reads/gsMask

raw_bam_dir="mapped"
downsamp_bam_dir1="${raw_bam_dir}/downsampled_ramp/2677917_reads"
downsamp_bam_dir2="${raw_bam_dir}/downsampled"
macs_out="./macs2_out/2677917_reads"
gs_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"

inpt1=Unknown_CE349-003R0001
inpt2=Unknown_CE349-003R0003
inpt="${downsamp_bam_dir2}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir2}/${inpt2}.dwnsmp.sorted.bam"

### One sample at a time
imm1=Unknown_CE349-003R0002
imm2=Unknown_CE349-003R0004
imm=(${imm1} ${imm2})

for rep in "${imm[@]}"; do

    if [[ ! -f "${macs_out}/${rep}_peaks.narrowPeak" ]]; then
         macs2 callpeak -t ${downsamp_bam_dir1}/${rep}.dwnsmp.sorted.bam -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n ${rep}
    fi

    # remove all peaks that overlap greenscreen
    bedtools intersect -v -wa -a ${macs_out}/${rep}_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/${rep}_peaks.narrowPeak

done

### Pooled samples
imm="${downsamp_bam_dir1}/${imm1}.dwnsmp.sorted.bam ${downsamp_bam_dir1}/${imm2}.dwnsmp.sorted.bam"

if [[ ! -f "${macs_out}/cocit_peaks.narrowPeak" ]]; then
     macs2 callpeak -t ${imm} -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n "cocit"
fi

# remove all peaks that overlap greenscreen
bedtools intersect -v -wa -a ${macs_out}/cocit_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/cocit_peaks.narrowPeak


## CASE2: subsampling distributes reads differently to subsampled rep1 and rep2 of cocit samples
### NB: I can use previously downsampled files
#mkdir -p ./mapped/downsampled_ramp/r1_4405927_r2_949907_reads
#cp ./mapped/downsampled_ramp/4405927_reads/Unknown_CE349-003R0002.dwnsmp.sorted.bam ./mapped/downsampled_ramp/r1_4405927_r2_949907_reads/
#cp ./mapped/downsampled_ramp/4405927_reads/Unknown_CE349-003R0002.dwnsmp.sorted.bam.bai ./mapped/downsampled_ramp/r1_4405927_r2_949907_reads/
#cp ./mapped/downsampled_ramp/949907_reads/Unknown_CE349-003R0004.dwnsmp.sorted.bam ./mapped/downsampled_ramp/r1_4405927_r2_949907_reads/
#cp ./mapped/downsampled_ramp/949907_reads/Unknown_CE349-003R0004.dwnsmp.sorted.bam.bai ./mapped/downsampled_ramp/r1_4405927_r2_949907_reads/

mkdir -p ./macs2_out/r1_4405927_r2_949907_reads
mkdir -p ./macs2_out/r1_4405927_r2_949907_reads/gsMask

raw_bam_dir="mapped"
downsamp_bam_dir1="${raw_bam_dir}/downsampled_ramp/r1_4405927_r2_949907_reads"
downsamp_bam_dir2="${raw_bam_dir}/downsampled"
macs_out="./macs2_out/r1_4405927_r2_949907_reads"
gs_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"

inpt1=Unknown_CE349-003R0001
inpt2=Unknown_CE349-003R0003
inpt="${downsamp_bam_dir2}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir2}/${inpt2}.dwnsmp.sorted.bam"

### One sample at a time
imm1=Unknown_CE349-003R0002
imm2=Unknown_CE349-003R0004
imm=(${imm1} ${imm2})

for rep in "${imm[@]}"; do

    if [[ ! -f "${macs_out}/${rep}_peaks.narrowPeak" ]]; then
         macs2 callpeak -t ${downsamp_bam_dir1}/${rep}.dwnsmp.sorted.bam -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n ${rep}
    fi

    # remove all peaks that overlap greenscreen
    bedtools intersect -v -wa -a ${macs_out}/${rep}_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/${rep}_peaks.narrowPeak

done

### Pooled samples
imm="${downsamp_bam_dir1}/${imm1}.dwnsmp.sorted.bam ${downsamp_bam_dir1}/${imm2}.dwnsmp.sorted.bam"

if [[ ! -f "${macs_out}/cocit_peaks.narrowPeak" ]]; then
     macs2 callpeak -t ${imm} -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n "cocit"
fi

# remove all peaks that overlap greenscreen
bedtools intersect -v -wa -a ${macs_out}/cocit_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/cocit_peaks.narrowPeak