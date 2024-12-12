PISA parse -rule 'CR,R2:4-10,barcodes.txt,CB,1;UR,R1:1-10' -1 test.R1.fq -2 test.R2.fq test.R1.fastq.gz test.R2.fastq.gz

split_barcode(){
        read1=$1
        read2=$2
        out=${1%/*}/01.rmdup_fq
  mkdir $out
  barcode=barcode.txt
  /home/yufeng/anaconda3/envs/R4.0/bin/python3.9 3demultiplexer.py $read1 $read2 $barcode $out
#  ls $out/*R[1,2].fastq | parallel -j 4 -N 2 rmdupfq
}
export -f split_barcode
#ls 01.m6a_seq2/0[1,2]*/*fastq.gz | parallel --tmpdir . -N 2 split_barcode
#ls 01.m6a_seq2/03*/*fastq.gz | parallel --tmpdir . -N 2 split_barcode


star_2_genomic(){
  index=/mnt/12/yuan_jianwen/index/01.hg38_star_index
  read1=$1
  read2=$2
  star_out=${read1%_R1*}_star_out/
  STAR \
  --alignEndsType EndToEnd \
  --genomeDir $index \
  --genomeLoad LoadAndKeep \
  --outBAMcompression 10 \
  --outFileNamePrefix $star_out \
  --outFilterMultimapNmax 1 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterType BySJout \
  --outReadsUnmapped Fastx \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMmode Full \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped Within \
  --outStd Log \
  --readFilesIn $read1 $read2 \
  --runMode alignReads \
  --readFilesCommand zcat \
  --runThreadN 20
}
export -f star_2_genomic

star_multi_map(){
  index=/mnt/12/yuan_jianwen/index/01.hg38_star_index
  read1=$1
  read2=$2
  star_out=${read1%/*}/01.multi_map/
  STAR \
  --alignEndsType EndToEnd \
  --genomeDir $index \
  --genomeLoad LoadAndKeep \
  --outBAMcompression 10 \
  --outFileNamePrefix $star_out \
  --outFilterMultimapNmax 100 \
  --winAnchorMultimapNmax 100 \
  --outSAMmultNmax 1 \
  --outMultimapperOrder Random \
  --runRNGseed 777 \
  --outFilterMultimapScoreRange 1 \
  --outFilterScoreMin 10 \
  --outFilterType BySJout \
  --outReadsUnmapped Within \
  --outSAMattrRGline ID:foo \
  --outSAMattributes All \
  --outSAMmode Full \
  --outSAMtype BAM Unsorted \
  --outSAMunmapped None \
  --outStd Log \
  --readFilesIn $read1 $read2 \
  --runMode alignReads \
  --runThreadN 4
  outbam=$star_out/Aligned.out.bam
  pre=$star_out/Aligned.out
  filterbam=$pre.filter.bam
  rmdB=$pre.randomer_rmdup.bam
  sambamba view -F "not (unmapped or mate_is_unmapped or secondary_alignment)" -f bam -o $filterbam $star_out/Aligned.out.bam
  python /mnt/13/yuan_jianwen/01.cell_stat_m6a/barcodecollasepe.py -b $filterbam -o $rmdB -m $pre.randomer_rmdup.metric && rm $filterbam
}
export -f star_multi_map

export -f star_multi_map
trim_alian(){
        read1=$1
        pre=${1%_R1.fastq*}
        read2=${pre}_R2.fastq
        sample=${pre##*/}
  barcode=barcode.txt

        R1_ad=`grep $sample $barcode | cut -f 3`
        R2_ad="NNNNNNNNNAGATCGGAAGAGCGTCGTGT"
        echo $R1_ad
        out1=${pre}_R1.trim.fq.gz
        out2=${pre}_R2.trim.fq.gz
  parse1=${pre}_R1.parse.fq
  parse2=${pre}_R2.parse.fq

###1 cut adapter
#  gunzip -d $read1.gz
#  gunzip -d $read2.gz
#  cutadapt -a $R1_ad -A "$R2_ad;min_overlap=14" -e 0.1 --no-indels -o $out1 -p $out2 -U 9 -m 15 $read1 $read2 && gzip $read1 && gzip $read2
#  /home/jia_wenqi/anaconda3/bin/fastp -i $out1 -o $out1.gz -I $out2 -O $out2.gz && mv $out1.gz $out1 && mv $out2.gz $out2
###2 alian
#  star_2_genomic $out1 $out2
###2.1 multi mapping reads for TE
#  star_multi_map ${out1%_R1*}_star_out/Unmapped.out.mate* && rm ${out1%_R1*}_star_out/Unmapped.out.mate*

###3
  outbam=${out1%_R1*}_star_out/Aligned.out.bam
  pre1=${outbam%.bam}
  filterbam=$pre1.filter.bam
  rmdB=$pre1.randomer_rmdup.bam
  subrmdB=${rmdB/.bam/.sub.bam}
#  sambamba view -F "not (unmapped or mate_is_unmapped or secondary_alignment)" -f bam -o $filterbam $outbam

#  python /mnt/13/yuan_jianwen/01.cell_stat_m6a/barcodecollasepe.py -b $filterbam -o $rmdB -m $pre1.randomer_rmdup.metric
#  sambamba sort -o $rmdB.tmp $rmdB && mv $rmdB.tmp $rmdB && rm $filterbam
#  samtools index $rmdB
#  subBam $rmdB 2000000 $subrmdB
###3.1 Overlap reads clipped
  clipbam=${rmdB/.bam/.clipOverReads.bam}
#  bam clipOverlap --in $rmdB --out $clipbam --stats > ${clipbam/bam/stat}

###3.2 read1Only bam
  read1B=${rmdB/.bam/.R1.bam}
#  sambamba view -F "first_of_pair" -f bam -o $read1B $rmdB


###4 m6a site depth stat
  m6aR=/mnt/1/lv_yuan/jingxia/workspace/02.m6A_seq2_M_109_NGN2/m6A_seq2_yjw/00.ref/hg38_m6a_site_from_m6aseq2.flank25bp.bed
  rRNABed=/mnt/1/lv_yuan/jingxia/workspace/02.m6A_seq2_M_109_NGN2/m6A_seq2_yjw/00.ref/gencode.v38.annotation.rRNA.bed
  depthout=$pre1.region.depth.bed
  subdepthout=${depthout/.bed/.sub.bed}
  depthoutR1=${depthout/.bed/.R1only.bed}
  subdepthoutR1=${subdepthout/.bed/.R1only.bed}
#  sambamba depth region  -L $m6aR -t 4 -o $depthout $rmdB
#  sambamba depth region  -L $m6aR -t 4 -o $subdepthout $subrmdB
#  sambamba depth region -F "first_of_pair" -L $m6aR -t 4 -o $depthoutR1 $rmdB
#  sambamba depth region -F "first_of_pair" -L $m6aR -t 4 -o $subdepthoutR1 $subrmdB
#  cut -f1-8 $depthout | intersectBed -a - -b $rRNABed -v -header > ${depthout/bed/rRNArm.bed}
#  cut -f1-8 $subdepthout | intersectBed -a - -b $rRNABed -v -header > ${subdepthout/bed/rRNArm.bed}
#  cut -f1-8 $depthoutR1 | intersectBed -a - -b $rRNABed -v -header > ${depthoutR1/bed/rRNArm.bed}
#  cut -f1-8 $subdepthoutR1 | intersectBed -a - -b $rRNABed -v -header > ${subdepthoutR1/bed/rRNArm.bed}

###5 gene level stat
#  gtf=/mnt/12/yuan_jianwen/hg38/gencode.v41.annotation.gtf
#  featureCounts -a $gtf -o ${rmdB/bam/count} -g gene_id --extraAttributes gene_name,gene_type \
#         -Q 10 -s 1 -O -d 30 --fraction \
#         --minOverlap 5 \
#         --largestOverlap -p -T 5 $rmdB
#  sambamba index $rmdB
#  sambamba depth region -L $m6aR -t 5 -F "mapping_quality > 1 and first_of_pair" -o ${depthout/bed/highq.bed} $rmdB

###TE element
TEout=${rmdB%/*}/02.TE
mkdir $TEout


mbam=${rmdB%/*}/01.multi_map/Aligned.out.randomer_rmdup.bam
mergeBam=${rmdB%.bam}.mulin.bam
#sambamba sort -o $mbam.tmp $mbam && mv $mbam.tmp $mbam
#sambamba merge ${rmdB%.bam}.mulin.bam $mbam $rmdB

/home/wu_xinyu/miniconda3/bin/TEcount --mode multi \
  -b $mergeBam \
  --GTF /mnt/12/yuan_jianwen/hg38/gencode.v32.annotation.gtf \
  --TE /mnt/12/yuan_jianwen/hg38/repeak_element/hg38_rmsk_TE.gtf \
  --project $TEout/TEcount \
  --stranded forward \
  --sortByPos

}
export -f trim_alian
#ls 0[1,2]*/01*/*[0,8]h*R1.fastq | parallel -j 4 trim_alian
ls 0[1]*/01*/*[0,8]h*R1.fastq | parallel -j 4 trim_alian
#ls 01.m6a_seq2/03*/01*/*R1.fastq | parallel -j 6 trim_alian


f1(){
        gtf="/mnt/12/yuan_jianwen/hg38/gencode.v34.annotation.gtf.gz"
        bam=$1
        pre=${1%.bam}
        rmdbam=$pre.rmdup.bam
#       sambamba markdup -r -t 3 $bam $pre.rmdup.bam
#       featureCounts --fraction -T 5 -p -t gene -g gene_name -a $gtf -M -O -o ${rmdbam%.bam}.count $1
        #featureCounts -T 5 -p -t gene -g gene_name -a $gtf -o ${pre}.count $bam
#  samtools sort -o $pre.sorted.bam     $bam
#  samtools view -b $pre.sorted.bam chr{1..22} chrX chrY chrM > $pre.temp && mv $pre.temp $pre.sorted.bam
#  samtools index ${pre}.sorted.bam
}
export -f f1

# step 3 count
#ls 0[3-4]/*new.sorted.bam | parallel f1
#allc="all.new.rmdup.count.txt"
#sed '1d' 01*/ERL.new.sorted.count | cut -f1 > $allc
#for i in 0*/*new.sorted.rmdup.count; do sed '1d' $i | cut -f7 | paste $allc - > ${allc}_1 && mv ${allc}_1 ${allc}; done


### step 4 call peak
# docker run -it --mount src=/mnt,target=/mnt,type=bind --workdir="${PWD}" brianyee/clipper:61d5456 bash
# reinstall clipper in /mnt/12/yuan_jianwen/tools/clipper-2.0.1 to avoid bug of chr missing
# pip install .
callpeak(){
        bam=$1
  prefix=${1%.bam}
        clipper -o ${prefix}.clipper.peak -b $bam -s hg38 --processors 3
}
export -f callpeak
#ls 01*/*bam | parallel --tmpdir . callpeak
