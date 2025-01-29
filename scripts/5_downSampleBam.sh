in_dir="mapped"
out_dir="mapped/downsampled"
seed=42
mkdir -p ${out_dir}

# function to find the minimum
# value in an array
minIndex(){
   arr=("$@")
   min_val=${arr[0]}
   min_idx=0
   for i in ${!arr[@]}; do
        cur_val=${arr[${i}]}
        if [[ ${cur_val} -lt ${min_val} ]]; then
                min_val=${arr[$i]}
                min_idx=${i}
        fi
   done

}


## With a single replicate: there is no need to downsample anything 
#single_rep=("FD_S_2019_Mock" 
#"FD_ft10_tsf1_S_2019_Mock")
#for samp in "${single_rep[@]}"; do
#  cp ${in_dir}/${samp}_R1.dupmark.sorted.bam \
#    ${out_dir}/${samp}_R1.dupmark.sorted.bam
#  cp ${in_dir}/${samp}_R1.dupmark.sorted.bam.bai \
#    ${out_dir}/${samp}_R1.dupmark.sorted.bam.bai
#done

# Downsample given two reps
dos2unix ./qc_out/chipqc_all/chip_readsize_fragsize.csv
my_array=( $(cut -d ";" -f4 ./qc_out/chipqc_all/chip_readsize_fragsize.csv | sort -u) )

for samp in "${my_array[@]}"; do
    echo "working with genotype ${samp}"
    
    rep1=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk 'NR==1{print $0}')
    rep2=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk 'NR==2{print $0}')

     echo "replicate1 is ${rep1}"
     echo "replicate2 is ${rep2}"

     depth1=`samtools view -c ${in_dir}/${rep1}.sorted.bam`
     depth2=`samtools view -c ${in_dir}/${rep2}.sorted.bam`

     arrName=(${depth1} ${depth2})
     minIndex "${arrName[@]}"
     
     let "min_idx = $min_idx + 1"

     for (( rep=1; rep<=2; rep++ )); do

         if [[ ${min_idx} -eq ${rep} ]]; then

            # select name of replicate with smallest read depth
            rep_min=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk -v pat="$rep" 'NR==pat{print $0}')
            echo "replicate with lower read number is $rep_min"

            # copy replicate with smallest read depth
            cp ${in_dir}/${rep_min}.sorted.bam ${out_dir}/${rep_min}.dwnsmp.sorted.bam
            cp ${in_dir}/${rep_min}.sorted.bam.bai ${out_dir}/${rep_min}.dwnsmp.sorted.bam.bai

          else

            # select name of replicate with smallest read depth
            rep_max=$(awk -F ";" -v pat="$samp""$" '$4~pat {print $1}' ./qc_out/chipqc_all/chip_readsize_fragsize.csv | awk -v pat="$rep" 'NR==pat{print $0}')
            echo "replicate with higher read number is $rep_max"
        
            # downsample these replicates
            biostar145820 --seed ${seed} -n ${min_val} -o ${out_dir}/${rep_max}.dwnsmp.bam ${in_dir}/${rep_max}.sorted.bam

	samtools sort -o ${out_dir}/${rep_max}.dwnsmp.sorted.bam ${out_dir}/${rep_max}.dwnsmp.bam && rm ${out_dir}/${rep_max}.dwnsmp.bam

	samtools index ${out_dir}/${rep_max}.dwnsmp.sorted.bam

        fi
    done
done

