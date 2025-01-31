raw_bam_dir="mapped"
downsamp_bam_dir="${raw_bam_dir}/downsampled"
bam_suffix="dwnsmp.sorted.bam"
macs2_out="macs2_out"
gs_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"
# average basepair q-value threshold (log10)
q=10

# make macs2 output directory
mkdir -p ${macs2_out}/noMask_qval${q}
mkdir -p ${macs2_out}/gsMask_qval${q}

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

    #chip_fsize=fragment size for each experiment, but I can let macs2 guess it, since it's all same experiment

    echo "input samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
    echo "immunoprecipitated samples are ${imm1} and ${imm2}" #I call peaks on individual samples

    # run MACS2 on the two reps

    for rep in "${imm[@]}"; do

                if [[ ! -f "${macs2_out}/${rep}_peaks.narrowPeak" ]]; then
                   macs2 callpeak -t ${raw_bam_dir}/${rep}.sorted.bam -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 150 -g 119482012 --outdir ${macs2_out} -n ${rep}
                fi

                # remove all peaks that do not have an average base pair q-value <=10^(-${q})
                awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q && $1!="ChrC" && $1!="ChrM"{print}' ${macs2_out}/${rep}_peaks.narrowPeak > ${macs2_out}/noMask_qval${q}/${rep}_peaks.narrowPeak

                # remove all peaks that overlap greenscreen
                bedtools intersect -v -wa -a ${macs2_out}/noMask_qval${q}/${rep}_peaks.narrowPeak -b ${gs_regions} > ${macs2_out}/gsMask_qval${q}/${rep}_peaks.narrowPeak
    done
done
    


# PEAK CALLING on POOLED CHIP SAMPLES

for trt in "${my_array_trt[@]}"; do
    echo "working with treatment ${trt}"
    
    arrName=$(awk -F ";" -v pat="$trt""$" '$4~pat' ./qc_out/chipqc_all/chip_readsize_fragsize.csv)

    inpt1=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    inpt2=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    inpt="${downsamp_bam_dir}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir}/${inpt2}.dwnsmp.sorted.bam"

    imm1=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==1{print $0}')
    imm2=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==2{print $0}')
    imm="${downsamp_bam_dir}/${imm1}.dwnsmp.sorted.bam ${downsamp_bam_dir}/${imm2}.dwnsmp.sorted.bam"

    echo "immunoprecipitated samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
    echo "gfp samples are ${imm}" #This time I call peaks on pooled immunoprecipitated samples

    # run MACS2
        
    if [[ ! -f "${macs2_out}/${trt}_peaks.narrowPeak" ]];then
        	macs2 callpeak -t ${imm} -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 150 -g 119482012 --outdir ${macs2_out} -n ${trt}
	fi
        
        awk -F"\t" -v q=${q} 'BEGIN{OFS="\t"} $9>=q && $1!="ChrC" && $1!="ChrM"{print}' ${macs2_out}/${trt}_peaks.narrowPeak > ${macs2_out}/noMask_qval${q}/${trt}_peaks.narrowPeak

        bedtools intersect -v -wa -a ${macs2_out}/noMask_qval${q}/${trt}_peaks.narrowPeak -b ${gs_regions} > ${macs2_out}/gsMask_qval${q}/${trt}_peaks.narrowPeak

done
