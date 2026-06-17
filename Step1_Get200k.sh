cat ../../../ATACPeaks/Monkey_dTAC.bed | awk -F'\t' -v OFS='\t' '{print $1,$2,$3,$1"_"$2"_"$3}' | sortBed > region.txt

cat Monkey_TSS.txt | sed '1d' | awk 'BEGIN{OFS="\t"}{print $3,$4-202000,$4+202000,$1}'| awk '{$2=$2<0?1:$2}1'|tr ' ' '\t' > Monkey_Promoter_200k.bed
cat Monkey_Promoter_200k.bed | awk -v OFS='\t' '{print $1,$2,$3,$4}' | sortBed > a;mv a Monkey_Promoter_200k.bed
bedtools intersect -a region.txt -b Monkey_Promoter_200k.bed -wa -wb | cut -f 4,8 > peak_gene.txt
bedtools intersect -a region.txt -b Monkey_Promoter_200k.bed -wa -wb | awk '{print $3-$7+100000}'|awk '{$1=$1<0?-$1:$1}1' > dis
