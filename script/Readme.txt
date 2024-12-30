# conda activate anno
# 1. 基因组比对首先需要屏蔽基因组上的重复序列。
mkdir 02.RepeatMasker

sh 01.run_RepeatMasker.sh
"""
#conda activate anno
RepeatMasker -pa 12 -species Brassicaceae -e rmblast -xsmall -html -s -gff Ath.genome.fa -dir repeat_result_Ath 1>Ath_log.o.txt 2>Ath_log.e.txt &
RepeatMasker -pa 12 -species Brassicaceae -e rmblast -xsmall -html -s -gff Bca.genome.fa -dir repeat_result_Bca 1>Bca_log.o.txt 2>Bca_log.e.txt &
RepeatMasker -pa 12 -species Brassicaceae -e rmblast -xsmall -html -s -gff Bju.genome.fa -dir repeat_result_Bju 1>Bju_log.o.txt 2>Bju_log.e.txt &
RepeatMasker -pa 12 -species Brassicaceae -e rmblast -xsmall -html -s -gff Bna.genome.fa -dir repeat_result_Bna 1>Bna_log.o.txt 2>Bna_log.e.txt &
"""

# 将fa文件转换生成2bit文件，计算每条染色体的长度，用于后续的分析。
## 根据名字把基因组文件拆分
faSplit byName Ath.genome.sm.fa Ath_fa
## 统计基因组信息(conoda activate anno)
faToTwoBit Ath.genome.sm.fa Ath.genome.sm.fa.2bit
faSize Ath.genome.sm.fa -detailed > Ath.genome.sm.fa.faSize

# 2. 使用 LASTZ 对基因组进行两两比对(首先对去除基因组重复区域后的基因组保留染色体基因)
mkdir 03.lastz

sh 01.lastz.sh
"""
lastz Chr1.fa ../02.RepeatMasker/repeat_result_Aar/Aar.genome.sm.fa O=400 E=30 T=1 H=2000 X=9400 L=3000 K=3000 --format=axt --ambiguous=n --ambiguous=iupac > Aar_Chr1.axt &
lastz Chr2.fa ../02.RepeatMasker/repeat_result_Aar/Aar.genome.sm.fa O=400 E=30 T=1 H=2000 X=9400 L=3000 K=3000 --format=axt --ambiguous=n --ambiguous=iupac > Aar_Chr2.axt &
lastz Chr3.fa ../02.RepeatMasker/repeat_result_Aar/Aar.genome.sm.fa O=400 E=30 T=1 H=2000 X=9400 L=3000 K=3000 --format=axt --ambiguous=n --ambiguous=iupac > Aar_Chr3.axt &
lastz Chr4.fa ../02.RepeatMasker/repeat_result_Aar/Aar.genome.sm.fa O=400 E=30 T=1 H=2000 X=9400 L=3000 K=3000 --format=axt --ambiguous=n --ambiguous=iupac > Aar_Chr4.axt &
lastz Chr5.fa ../02.RepeatMasker/repeat_result_Aar/Aar.genome.sm.fa O=400 E=30 T=1 H=2000 X=9400 L=3000 K=3000 --format=axt --ambiguous=n --ambiguous=iupac > Aar_Chr5.axt &
"""

# 3. Chaining 从axt文件得到chain文件
mkdir 04.chain 
	mkdir 2bit ; cd 2bit ;ln -s ../../02.RepeatMasker/repeat_result_*/*.2bit ./ ; cd ..
	mkdir axt ; cd axt ;ln -s ../../03.lastz/*axt  ./ ; cd ..
	mkdir chain

	sh chr1.sh ; sh chr2.sh ; sh chr3.sh ; sh chr4.sh ; chr5.sh
"""
nohup axtChain axt/Aal_Chr1.axt ../02.RepeatMasker/repeat_result_Ath/Ath.genome.sm.fa.2bit 2bit/Aal.genome.sm.fa.2bit chain/Aal_Chr1.chain -minScore=3000 -linearGap=medium &
nohup axtChain axt/Aal_Chr2.axt ../02.RepeatMasker/repeat_result_Ath/Ath.genome.sm.fa.2bit 2bit/Aal.genome.sm.fa.2bit chain/Aal_Chr2.chain -minScore=3000 -linearGap=medium &
nohup axtChain axt/Aal_Chr3.axt ../02.RepeatMasker/repeat_result_Ath/Ath.genome.sm.fa.2bit 2bit/Aal.genome.sm.fa.2bit chain/Aal_Chr3.chain -minScore=3000 -linearGap=medium &
nohup axtChain axt/Aal_Chr4.axt ../02.RepeatMasker/repeat_result_Ath/Ath.genome.sm.fa.2bit 2bit/Aal.genome.sm.fa.2bit chain/Aal_Chr4.chain -minScore=3000 -linearGap=medium &
nohup axtChain axt/Aal_Chr5.axt ../02.RepeatMasker/repeat_result_Ath/Ath.genome.sm.fa.2bit 2bit/Aal.genome.sm.fa.2bit chain/Aal_Chr5.chain -minScore=3000 -linearGap=medium &
"""

# 4. 合并chain文件，并过滤
mkdir 05.chainMergeSort

sh 01.chainMergeSort.sh
"""
chainMergeSort ../04.chain/chain/Aal*.chain > Aal.chain &
chainMergeSort ../04.chain/chain/Aar*.chain > Aar.chain &
chainMergeSort ../04.chain/chain/Aly*.chain > Aly.chain &
chainMergeSort ../04.chain/chain/Ane*.chain > Ane.chain &
chainMergeSort ../04.chain/chain/Aru*.chain > Aru.chain &
"""

# 5. 连接chain，并为下一步做准备
mkdir 06.chainPreNet

ln -s ../02.RepeatMasker/*/*genome.sm.fa.faSize ./
sh 01.chainPreNet.sh
"""
chainPreNet ../05.chainMergeSort/Aal.chain Ath.genome.sm.fa.faSize Aal.genome.sm.fa.faSize Aal.pre.chain &
chainPreNet ../05.chainMergeSort/Aar.chain Ath.genome.sm.fa.faSize Aar.genome.sm.fa.faSize Aar.pre.chain &
chainPreNet ../05.chainMergeSort/Aly.chain Ath.genome.sm.fa.faSize Aly.genome.sm.fa.faSize Aly.pre.chain &
chainPreNet ../05.chainMergeSort/Ane.chain Ath.genome.sm.fa.faSize Ane.genome.sm.fa.faSize Ane.pre.chain &
chainPreNet ../05.chainMergeSort/Aru.chain Ath.genome.sm.fa.faSize Aru.genome.sm.fa.faSize Aru.pre.chain &
"""

# 6. chainNet 排序
mkdir 07.chainNet
ln -s ../02.RepeatMasker/*/*genome.sm.fa.faSize ./
ln -s ../06.chainPreNet/*pre.chain ./

sh 01.chainNet.sh
"""
chainNet Aal.pre.chain -minSpace=1 Ath.genome.sm.fa.faSize Aal.genome.sm.fa.faSize stdout /dev/null | netSyntenic stdin Aal_noClass.net
chainNet Aar.pre.chain -minSpace=1 Ath.genome.sm.fa.faSize Aar.genome.sm.fa.faSize stdout /dev/null | netSyntenic stdin Aar_noClass.net
chainNet Aly.pre.chain -minSpace=1 Ath.genome.sm.fa.faSize Aly.genome.sm.fa.faSize stdout /dev/null | netSyntenic stdin Aly_noClass.net
chainNet Ane.pre.chain -minSpace=1 Ath.genome.sm.fa.faSize Ane.genome.sm.fa.faSize stdout /dev/null | netSyntenic stdin Ane_noClass.net
chainNet Aru.pre.chain -minSpace=1 Ath.genome.sm.fa.faSize Aru.genome.sm.fa.faSize stdout /dev/null | netSyntenic stdin Aru_noClass.net
"""

# 7. Maffing 并生成 maf 文件
mkdir 08.Maffing ;cd 08.Maffing
ln -s ../02.RepeatMasker/*/*.genome.sm.fa.2bit ./
ln -s ../02.RepeatMasker/*/*genome.sm.fa.faSize ./
ln -s ../06.chainPreNet/*pre.chain ./
ln -s ../07.chainNet/*_noClass.net ./

sh 01.Maffing.sh
"""
netToAxt Aal_noClass.net Aal.pre.chain Ath.genome.sm.fa.2bit Aal.genome.sm.fa.2bit stdout | axtSort stdin Ath_Aal.axt
netToAxt Aar_noClass.net Aar.pre.chain Ath.genome.sm.fa.2bit Aar.genome.sm.fa.2bit stdout | axtSort stdin Ath_Aar.axt
netToAxt Aly_noClass.net Aly.pre.chain Ath.genome.sm.fa.2bit Aly.genome.sm.fa.2bit stdout | axtSort stdin Ath_Aly.axt
netToAxt Ane_noClass.net Ane.pre.chain Ath.genome.sm.fa.2bit Ane.genome.sm.fa.2bit stdout | axtSort stdin Ath_Ane.axt
netToAxt Aru_noClass.net Aru.pre.chain Ath.genome.sm.fa.2bit Aru.genome.sm.fa.2bit stdout | axtSort stdin Ath_Aru.axt
"""
sh 02.maf.sh
"""
axtToMaf Ath_Aal.axt Ath.genome.sm.fa.faSize Aal.genome.sm.fa.faSize Ath_Aal.maf -tPrefix=Ath. -qPrefix=Aal.
axtToMaf Ath_Aar.axt Ath.genome.sm.fa.faSize Aar.genome.sm.fa.faSize Ath_Aar.maf -tPrefix=Ath. -qPrefix=Aar.
axtToMaf Ath_Aly.axt Ath.genome.sm.fa.faSize Aly.genome.sm.fa.faSize Ath_Aly.maf -tPrefix=Ath. -qPrefix=Aly.
axtToMaf Ath_Ane.axt Ath.genome.sm.fa.faSize Ane.genome.sm.fa.faSize Ath_Ane.maf -tPrefix=Ath. -qPrefix=Ane.
axtToMaf Ath_Aru.axt Ath.genome.sm.fa.faSize Aru.genome.sm.fa.faSize Ath_Aru.maf -tPrefix=Ath. -qPrefix=Aru.
"""

# 8. MULTIZ 合并
mkdir 09.MULTIZ

ln -s ../08.Maffing/*.maf ./
roast - T=`pwd` R=30 M=20 E=Ath "((Aal Aar Aly Ane Aru Bca Bju Bna Bni Bol Bra Cen Chi Cru Csa Dni Ech Iin Lma Mpy Ovi Raq Rbr Rsa Sal Sar Spa Tar) Ath)" *.*.maf ref_species_mulitway.maf > roast.sh
sh roast.sh

# 9. i-Adhore 共线性比对
mkdir 10.i-Adhore

ln -s ../09.MULTIZ/table/new_ref_species_mulitway.table ./
sh 1.1.gff.tmp.sh
"""
grep Ath new_ref_species_mulitway.table | sort -k3,3 -k5n > tmp/Ath_gff.tmp
grep Aal new_ref_species_mulitway.table | sort -k3,3 -k5n > tmp/Aal_gff.tmp
grep Aar new_ref_species_mulitway.table | sort -k3,3 -k5n > tmp/Aar_gff.tmp
grep Aly new_ref_species_mulitway.table | sort -k3,3 -k5n > tmp/Aly_gff.tmp
grep Ane new_ref_species_mulitway.table | sort -k3,3 -k5n > tmp/Ane_gff.tmp
"""
sh 1.2.families.csv.sh
"""
perl scripts/write_families.pl new_ref_species_mulitway.table > ./families.csv
"""
sh 1.3.mkdir_file.sh
"""
perl scripts/write_scaffold_lst_files.pl tmp/Aal_gff.tmp Aal &
perl scripts/write_scaffold_lst_files.pl tmp/Aar_gff.tmp Aar &
perl scripts/write_scaffold_lst_files.pl tmp/Aly_gff.tmp Aly &
perl scripts/write_scaffold_lst_files.pl tmp/Ane_gff.tmp Ane &
perl scripts/write_scaffold_lst_files.pl tmp/Aru_gff.tmp Aru &
"""
sh 1.4.write_the_data.ini
# ini file indicating the running parameters used in i-Adhore
# conda activate wgd
# ~/00.Software/01.conda/Miniconda3/envs/wgd/bin/i-adhore
"""
cat genome.lst datasetII.ini > data.ini
"""
sh 1.5.i-adhore.sh
"""
i-adhore data.ini
"""
sh 2.segments_table.sh
"""
echo 'Step 2 Start at'
date
export reference="Ath"
export nway_out="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/segments_table/n-way"
export segments="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/output/segments.txt"
export Step2="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/scripts/Step2.pl"
echo 'preparing n-way [n=3~17] segments'
mkdir ${nway_out}
for n in {3..17}
        do
                perl $Step2 $segments $reference `seq 2 $n` > ${nway_out}/${n}-way_segments.table
        done
echo 'Step 2 Finished at'
date
"""
sh 3.segments_gff3.sh
"""
echo 'Step 3 Start at'
date
export MAAs_table="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/new_ref_species_mulitway.table"
export segments_table="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/segments_table/n-way/17-way_segments.table"
export nway_out="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/segments_gff3/n-way"
export Step3="/data01/masterHome/lchunjin/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/09.i-Adhore/scripts/Step3.pl"

for n in {3..17}
        do
                perl scripts/Step3.pl segments_table/n-way/${n}-way_segments.table $MAAs_table | sort -u > ${nway_out}/${n}way_segments.gff3
        done
"""



# 10. 运行phastCons
mkdir 10.phastCons; cd 10.phastCons
sh 1.phyloFit.sh
即
phyloFit -i MAF Ath_mulitway.maf --tree new_SpeciesTree_rooted_str.txt

sh 2.mafSplit.sh
mkdir maf
mafSplit _.bed maf/ Ath_mulitway.maf -byTarget -useFullSequenceName

# 运行phastCons计算保守性分数。
mkdir wig
mkdir bed
for i in maf/*.maf; do x=`basename $i .maf`; phastCons $i phyloFit.mod --target-coverage 0.2 --expected-length 80 --rho 0.4 --msa-format MAF --seqname $x --most-conserved bed/$x.most-cons.bed > wig/$x.wig; done
cat wig/*.wig >> Ath_mulitway.wig
cat bed/*.most-cons.bed >> Ath_mulitway_most-cons.bed


sh 3.wigToBigWig.sh
1.使用wigToBigWig将wig文件转换为二进制的bigWig文件后，就能在Genome Browser中展示了。
wigToBigWig Ath_mulitway.wig Ath.genome.sm.fa.faSize Ath.bw

sh 4.bigWigAverageOverBed.sh
2.假设给定一条序列比如lncRNA，我们想计算它的平均phastCons分数，可以使用bigWigAverageOverBed实现。
bigWigAverageOverBed Ath.bw Ath_mulitway_most-cons.bed Ath.tab
3.假设给定一条序列比如lncRNA，我们想计算它的平均phastCons分数，可以使用bigWigAverageOverBed实现。



# 10.cis-regulatory_conserved_noncoding
mkdir cis-regulatory_conserved_noncoding ; cd cis-regulatory_conserved_noncoding

1
mkdir maf ; cd maf
ln -s ../../MULTIZ/ref_species_mulitway.maf ./
grep -v "#" ref_species_mulitway.maf > 16_way_mulitway.maf

2
# conda activate rnaseq
mkdir collinear_segments ; cd collinear_segments
ln -s ../../i-Adhore/segments_gff3/n-way/17way_segments.gff3 ./
grep Ath 17way_segments.gff3 > Ath_17way_segments.gff3
convert2bed -i gff -o bed <Ath_17way_segments.gff3> Ath_17way_segments.bed
convert2bed -i gff -o bed <17way_segments.gff3> 17way_segments.bed
#awk -F '\t' '{print$1"\t"$2"\t"$3}' Ath_17way_segments.bed |sort -u |sort -k2n > new_Ath_17way_segments.bed

3
mkdir cds_bed ; cd cds_bed
conda activate jcvi
ln -s ../../../01.data/01.Brassicaceae/30.Arabidopsis_thaliana/Arabidopsis_thaliana.gff ./
sh 1.jcvi_cds_bed.sh
 即 python -m jcvi.formats.gff bed --type=CDS --key=ID Arabidopsis_thaliana.gff > Ath.cds.bed
#awk -F '\t' '{print$1"\t"$2"\t"$3}' Ath.cds.bed |sort -u |sort -k2n > new_Ath.cds.bed

4
mkdir mostCons ; cd mostCons
ln -s ../../phastCons/Ath_mulitway_most-cons.bed ./allMostCons.bed
mkdir chr ; cd chr
ln -s ~/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/phastCons/bed/* ./
ln -s ~/03.Pipeline/01.Comparative_genomes/01.Brassicaceae/12.cis_rcne/phastCons/wig/* ./

5
cd ..
sh run.sh


