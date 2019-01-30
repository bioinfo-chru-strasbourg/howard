
#!/bin/bash
#################################
## Extract contig 
#################################


# $1 VCF
# $2 Ref.fa.fai
# $3 force/replace contig 

existContig=$($BCFTOOLS view -h $1 | grep "##contig=" -c)
existCHRVCF=$($BCFTOOLS view -H $1 | head -n1 | cut -f1 | grep chr -c)
existCHRFAI=$(cat $2 | head -n1 | cut -f1 | grep chr -c)
if [ ! -z $3 ] || !(($existContig)); then needContig=1; fi;

#echo "needConfig=$needContig"
#echo "existCHRVCF=$existCHRVCF"
#echo "existCHRFAI=$existCHRFAI"

#head $1

if (($needContig)); then 

	# Header without contig
	$BCFTOOLS view -h $1 | grep ^##contig= -v | grep ^## > $1.header
	
	# contig section

	if (($existCHRVCF)); then
		if (($existCHRFAI)); then
			cat $2 | awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}'  >> $1.header
		else
			cat $2 | awk '{printf("##contig=<ID=chr%s,length=%d>\n",$1,$2);}'  >> $1.header
		fi
	else
		if (($existCHRFAI)); then
			cat $2 | sed s/chr//gi | awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}'  >> $1.header
		else
			cat $2 | awk '{printf("##contig=<ID=%s,length=%d>\n",$1,$2);}' >> $1.header
		fi
	fi;

	# CHROM...
	$BCFTOOLS view -h $1 | grep ^##contig= -v | grep ^#CHROM  >> $1.header

	# Varaints
	#$BCFTOOLS view -H $1
	$BCFTOOLS reheader $1 -h $1.header

fi;



