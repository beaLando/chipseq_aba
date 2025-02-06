# SET VARIABLES

raw_bam_dir="mapped"
downsamp_bam_dir="${raw_bam_dir}/downsampled"
bam_suffix="dwnsmp.sorted.bam"
macs2_out="macs2_out"
gs_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"
glIDR=1 #-log10(0.1)



# PEAK CALLING on SINGLE REPLICATES

conda activate chippy

# make macs2 output directory
mkdir -p ${macs2_out}
mkdir -p ${macs2_out}/gsMask

dos2unix ./qc_out/chipqc_all/chip_readsize_fragsize.csv
my_array_trt=( $(cut -d ";" -f4 ./qc_out/chipqc_all/chip_readsize_fragsize.csv |  sed 's/Unknown_CE349-003R000*_//' | sed 's/^i//' | sed 's/^gfp//' | sort -u) ) 


# PEAK CALLING on SINGLE CHIP SAMPLES

for trt in "${my_array_trt[@]}"; do
    echo "working with treatment ${trt}"
    
    arrName=$(awk -F ";" -v pat="$trt""$" '$4~pat' ./qc_out/chipqc_all/chip_readsize_fragsize.csv)

    inpt1=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    inpt2=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    inpt="${downsamp_bam_dir}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir}/${inpt2}.dwnsmp.sorted.bam"

    imm1=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    imm2=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    imm=(${imm1} ${imm2})

    echo "input samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
    echo "immunoprecipitated samples are ${imm1} and ${imm2}" #I call peaks on individual samples

    # run MACS2 on the two reps

    for rep in "${imm[@]}"; do

                if [[ ! -f "${macs2_out}/${rep}_peaks.narrowPeak" ]]; then
                   macs2 callpeak -t ${raw_bam_dir}/${rep}.sorted.bam -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs2_out} -n ${rep}
                fi

                bedtools intersect -v -wa -a ${macs2_out}/${rep}_peaks.narrowPeak -b ${gs_regions} > ${macs2_out}/gsMask/${rep}_peaks.narrowPeak

    done
done
    


# PEAK CALLING on POOLED REPLICATES

for trt in "${my_array_trt[@]}"; do
    echo "working with treatment ${trt}"
    
    arrName=$(awk -F ";" -v pat="$trt""$" '$4~pat' ./qc_out/chipqc_all/chip_readsize_fragsize.csv)

    inpt1=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    inpt2=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    inpt="${downsamp_bam_dir}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir}/${inpt2}.dwnsmp.sorted.bam"

    imm1=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    imm2=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    imm="${downsamp_bam_dir}/${imm1}.dwnsmp.sorted.bam ${downsamp_bam_dir}/${imm2}.dwnsmp.sorted.bam"

    echo "input samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
    echo "immunoprecipitated samples are ${imm}" #This time I call peaks on pooled immunoprecipitated samples

    # run MACS2
        
    if [[ ! -f "${macs2_out}/${trt}_peaks.narrowPeak" ]];then
        	macs2 callpeak -t ${imm} -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs2_out} -n ${trt}
    fi

    bedtools intersect -v -wa -a ${macs2_out}/${trt}_peaks.narrowPeak -b ${gs_regions} > ${macs2_out}/gsMask/${trt}_peaks.narrowPeak

done

conda deactivate



# FINDING IDR PEAKS AND MASKING WITH GREENSCREEN

conda activate idr_chip

## Make IDR output directory
mkdir -p macs2_out/IDR

## Sort peak by -log10(p-value)
sort -k8,8nr macs2_out/Unknown_CE349-003R0002_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0002_peaks_sortedbylP.narrowPeak 
sort -k8,8nr macs2_out/Unknown_CE349-003R0004_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0004_peaks_sortedbylP.narrowPeak
sort -k8,8nr macs2_out/Unknown_CE349-003R0006_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0006_peaks_sortedbylP.narrowPeak
sort -k8,8nr macs2_out/Unknown_CE349-003R0008_peaks.narrowPeak > macs2_out/Unknown_CE349-003R0008_peaks_sortedbylP.narrowPeak

## IDR
idr --samples macs2_out/Unknown_CE349-003R0002_peaks_sortedbylP.narrowPeak macs2_out/Unknown_CE349-003R0004_peaks_sortedbylP.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file macs2_out/IDR/cocit_idr.peaks \
--plot \
--log-output-file macs2_out/IDR/cocit.idr.log

idr --samples macs2_out/Unknown_CE349-003R0006_peaks_sortedbylP.narrowPeak macs2_out/Unknown_CE349-003R0008_peaks_sortedbylP.narrowPeak \
--input-file-type narrowPeak \
--rank p.value \
--output-file macs2_out/IDR/cocitaba_idr.peaks \
--plot \
--log-output-file macs2_out/IDR/cocitaba.idr.log

### I filter at global idr >= 0.1, but see:
#### https://groups.google.com/g/idr-discuss/c/Z4X_x2ajTe4
#### https://groups.google.com/g/idr-discuss/c/PZ6tiZ0C-kg
awk -F"\t" -v glIDR=${glIDR} 'BEGIN{OFS="\t"} $12>=glIDR{print}' macs2_out/IDR/cocit_idr.peaks > macs2_out/IDR/cocit_idr01.peaks
awk -F"\t" -v glIDR=${glIDR} 'BEGIN{OFS="\t"} $12>=glIDR{print}' macs2_out/IDR/cocitaba_idr.peaks > macs2_out/IDR/cocitaba_idr01.peaks

intersectBed -v -wa -a macs2_out/IDR/cocit_idr01.peaks -b $gs_regions > macs2_out/IDR/cocit_idr01_gsMask.peaks
intersectBed -v -wa -a macs2_out/IDR/cocitaba_idr01.peaks -b $gs_regions > macs2_out/IDR/cocitaba_idr01_gsMask.peaks

conda deactivate