
#!/bin/bash

bamtopairs(){
    ${SAMTOOLS} view ${1} | \
    awk -v OFS="\t" '{print $3, $4, int($2%128%64%32/16), int($2/128)}' | \
    awk 'NR%2{printf $0 "\t";next;}1' | \
    awk -v OFS="\t" '{
        sub("1", "-", $3)sub("0", "+", $3);
        sub("1", "-", $7)sub("0", "+", $7);
        print substr($1, 4), $2, $3, substr($5, 4), $6, $7}' | \
    awk '{
        sub("-", "0", $3)sub("+", "1", $3);
        sub("-", "0", $6)sub("+", "1", $6);
        print $0}' \
    > ${1##*/}.pairs
}

source ${1}
mkdir -p ${OUTPUT}
cd ${OUTPUT}
echo "" | tee ${LOGFILE}
echo "Welcome to run HiSIF Linux Shell Scripts..." | tee -a ${LOGFILE}
echo "Version: 1.0.0" | tee -a ${LOGFILE}
if [ ${PREPROCESSING}x = "True"x ]
then
    echo "Transfer BAM file to HiSIF input 6-column pairs" | tee -a ${LOGFILE}
    if [ ${HICPROBAMPATH} ]
    then
        BAMFILE=($(ls ${HICPROBAMPATH}/*_hg19.bwt2pairs.bam))
    elif [ ${BAMPATH} ]
    then
        BAMFILE=($(ls ${BAMPATH}/*.bam))
    fi
    for i in ${!BAMFILE[@]}
    do
        echo ${BAMFILE[i]} | tee -a ${LOGFILE}
        bamtopairs ${BAMFILE[i]}
    done
    echo "Combine all pairs to ${SAMPLE}.pairs" | tee -a ${LOGFILE}
    cat *.pairs > ${SAMPLE}.pairs
    echo "Make the directory split" | tee -a ${LOGFILE}
    mkdir -p split
    line=($(wc -l ${SAMPLE}.pairs))
    length=($(echo ${line[0]} | wc -L))
    echo "split files to split" | tee -a ${LOGFILE}
    split -l $[10**$length/100] -d ${SAMPLE}.pairs split/${SAMPLE}_HiC
    echo "Make the directory ${SAMPLE}" | tee -a ${LOGFILE}
    mkdir -p ${SAMPLE}
    echo "Preprocessing with proc..." | tee -a ${LOGFILE}
    ${HISIF}/bin/proc split ${SAMPLE} -t
fi

echo "Run HiSIF..." | tee -a ${LOGFILE}
for tvalue in $(seq ${THRESHOLD_FIRST} ${THRESHOLD_LAST})
do
    echo "Running HiSIF with the threshold: ${tvalue}" | tee -a ${LOGFILE}
    ${HISIF}/bin/HiSIF \
        -g ${REFERENCE_GENOME} \
        -c ${HISIF}/resources/${CUTTING_FRAGMENTS} \
        -w ${READ_LENGTH} ${CUTTING_SITE_EXTENT} ${BIN_SIZE} \
        -p ${POISSON_MIXTURE1} ${POISSON_MIXTURE2} \
        -t ${tvalue} \
        -i ${ITERATIONS} \
        -f ${FRAGMENT_SIZE} \
        ${SAMPLE} >> ${LOGFILE} 2>&1
done
echo "...Running HiSIF is done..." | tee -a ${LOGFILE}
echo "...END..." | tee -a ${LOGFILE}

