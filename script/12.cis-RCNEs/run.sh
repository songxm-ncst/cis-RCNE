##Edited by songht@20201222  cmb@bnu Email: songhongtao@bnu.edu.cn##
#####depandency tools#####
# PHAST v1.4 #
# bedtools v2.31.0 #
# conda activate wgd

export wkdir="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs"
export maf="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs/01.maf/new1_ref_species_mulitway.maf"
export scripts="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs/scripts"
export coll_segment="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs/02.collinear_segments/Ath_29way_segments.bed" # 12way collinear segments with bed file
export CDS="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs/03.cds_bed/Ath.cds.bed" #coding sequence bed file in cucumber genome
export mostCons="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/15.cis_rcne_hpj/12.cis-RCNEs/04.mostCons/merge_filter_Ath_mulitway_most-cons.bed"
export len_cut=20 #set the minimal length of elements with base pair as unit #
export iden_cut=0.6 #set the minimal alignment identitiy cutoff for screening elements #

## Step1 Retain the noncoding elements from mostCons elements ##
bedtools intersect -a $mostCons -b $CDS -v > $wkdir/mostCons_nocoding.bed 

## Step2 Retain the cis-regulatory elements that located in collinear segments ##
bedtools intersect -a $wkdir/mostCons_nocoding.bed -b $coll_segment -wa |sort -u > $wkdir/mostCons_nocoding_inColl.bed

## Step3 Retain >= minimal length elements from mostCons_nocoding_inColl elements ##
awk -v len_cut=$len_cut '{if (($3-$2+1)>= len_cut){print $_}}' $wkdir/mostCons_nocoding_inColl.bed |sort -u > $wkdir/mostCons_nocoding_inColl_$len_cut.bed

## Step4 Using maf_parse to extract elements alignments based on 16way alignments ##
maf_parse  -g $wkdir/mostCons_nocoding_inColl_$len_cut.bed $maf  > $wkdir/elements.submaf

## Step5 Using maf_parse to extract elements alignments based on 16way alignments ##
perl $scripts/maf_parse.pl  $wkdir/elements.submaf $len_cut $iden_cut  > $wkdir/alternate_cisRCNEs.maf

## Step6 Retain elements with minimal alignment identity from elements ##
#maf_parse $maf -g $wkdir/mostCons_nocoding_inColl_$len_cut.bed > $wkdir/new.bed
#maf_parse  -g $wkdir/mostCons_nocoding_$len_cut.bed $maf  > $wkdir/elements.submaf

## Step6 Choose 29 genome segment ### 只保留每一块等于29个物种的区块 #
python 1.py
## 结果 cisRCNEs.maf out2.maf

###  保守原件 ./04.mostCons/allMostCons.bed


# 根据 mostCons_nocoding_inColl.bed 给 cisRCNEs_pro.bed 重新命名
grep Ath cisRCNEs.maf |sed 's/[ ][ ]*/ /g' |awk -F " " '{print$2,$3,$3+$4,"0",$5}' |sed 's/ /\t/g' > cisRCNEs_pro.bed
python name.py mostCons_nocoding_inColl_20.bed cisRCNEs_pro.bed |sort -u > cisRCNEs.bed

# diff mostCons_nocoding_inColl_20.bed cisRCNEs.bed 检查两个文件是否一致




