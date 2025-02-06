# Immunoprecipitated sample files must be from step 03a, while I can keep on using inputs from step 03

conda activate chippy

raw_bam_dir="mapped"
downsamp_bam_dir1="${raw_bam_dir}/downsampled_ramp"
downsamp_bam_dir2="${raw_bam_dir}/downsampled"
bam_suffix="dwnsmp.sorted.bam"
gs_regions="ara_refs/arabidopsis_greenscreen_20inputs.bed"


dos2unix ./qc_out/chipqc_all/chip_readsize_fragsize.csv
my_array_trt=( $(cut -d ";" -f4 ./qc_out/chipqc_all/chip_readsize_fragsize.csv |  sed 's/Unknown_CE349-003R000*_//' | sed 's/^i//' | sed 's/^gfp//' | sort -u) ) 


# PEAK CALLING on SINGLE CHIP SAMPLES

for d in $downsamp_bam_dir1/*; do
    if [ -d "$d" ]; then

    d_base="$(basename $d)"
    echo "working with immunoprecipitated samples downsampled to $d_base"
    
    macs_out=./macs2_out/${d_base}
    mkdir -p ${macs_out}
    mkdir -p ${macs_out}/gsMask

        for trt in "${my_array_trt[@]}"; do
        echo "working with treatment ${trt}"
    
        arrName=$(awk -F ";" -v pat="$trt""$" '$4~pat' ./qc_out/chipqc_all/chip_readsize_fragsize.csv)

        inpt1=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==1{print $0}')
        inpt2=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==2{print $0}')
        inpt="${downsamp_bam_dir2}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir2}/${inpt2}.dwnsmp.sorted.bam"

        imm1=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==1{print $0}')
        imm2=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==2{print $0}')
        imm=(${imm1} ${imm2})

        echo "input samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
        echo "immunoprecipitated samples are ${imm1} and ${imm2}" #I call peaks on individual samples

        # run MACS2 on the two reps

        for rep in "${imm[@]}"; do

                    if [[ ! -f "${macs_out}/${rep}_peaks.narrowPeak" ]]; then
                       macs2 callpeak -t ${d}/${rep}.dwnsmp.sorted.bam -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n ${rep}
                    fi

                    # remove all peaks that overlap greenscreen
                    bedtools intersect -v -wa -a ${macs_out}/${rep}_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/${rep}_peaks.narrowPeak
            done
        done
    fi
done



# PEAK CALLING on POOLED CHIP SAMPLES

for d in $downsamp_bam_dir1/*; do
        if [ -d "$d" ]; then

        d_base="$(basename $d)"
        echo "working with immunoprecipitated samples downsampled to $d_base"
    
        macs_out=./macs2_out/${d_base} 

        for trt in "${my_array_trt[@]}"; do
            echo "working with treatment ${trt}"
    
            arrName=$(awk -F ";" -v pat="$trt""$" '$4~pat' ./qc_out/chipqc_all/chip_readsize_fragsize.csv)

            inpt1=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==1{print $0}')
            inpt2=$(echo "$arrName" | awk -F ";" -v pat="_i" '$4~pat {print $1}' | awk 'NR==2{print $0}')
            inpt="${downsamp_bam_dir2}/${inpt1}.dwnsmp.sorted.bam ${downsamp_bam_dir2}/${inpt2}.dwnsmp.sorted.bam"

            imm1=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==1{print $0}')
            imm2=$(echo "$arrName" | awk -F ";" -v pat="_gfp" '$4~pat {print $1}' | awk 'NR==2{print $0}')
            imm="${d}/${imm1}.dwnsmp.sorted.bam ${d}/${imm2}.dwnsmp.sorted.bam"

            echo "input samples are ${inpt}" #I always use non-immunoprecipitated samples together as control
            echo "immunoprecipitated samples are ${imm}" #This time I call peaks on pooled immunoprecipitated samples

            # run MACS2
        
            if [[ ! -f "${macs_out}/${trt}_peaks.narrowPeak" ]];then
        	    macs2 callpeak -t ${imm} -c ${inpt} -f BAM --keep-dup auto --nomodel --extsize 170 -g 101274395 -p 1e-3 --outdir ${macs_out} -n ${trt}
	    fi
        
            bedtools intersect -v -wa -a ${macs_out}/${trt}_peaks.narrowPeak -b ${gs_regions} > ${macs_out}/gsMask/${trt}_peaks.narrowPeak

        done
    fi
done


conda deactivate