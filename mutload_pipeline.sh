#!/bin/sh
ref="<file path to reference genome>"
bed="<file path to bed file>"
picard="<file path to picard>" # alternatively, gatk4 can be used
bwa="<file path to bwa>"
python2="<file path to python2>"
UCM="<file path to UnifiedConsensusMaker.py>"
samtools="<file path to samtools>"
bcftools="<file path to bcftools>"
fastq="fastq" # if applicable, specify a different directory containing source fastq-files
ext1="_R1.fastq.gz" # adjust to the general part of the filename of read1
ext2="_R2.fastq.gz" # adjust to the general part of the filename of read2
# note that output will be saved to the current working directory


mkdir tmp
mkdir bam_unaligned
mkdir bam
mkdir bamdump
mkdir DCSdump
mkdir DCS_output
mkdir doc

if [ ! -d $fastq ]
then
  echo ERROR: fastq directory not found.
  exit
fi

for i in $fastq/*$ext1
	 	
	do
	
	bn=`basename $i`
	fname=${bn/$ext1/}
        r1=$fastq/$fname$ext1
        r2=$fastq/$fname$ext2
	
	echo `date` " [$$] Trimming for: " $fname	
	paste -- <(zgrep "" $r1) <(zgrep "" $r2) | perl -ne \
		'
		++$count;
		if ($count == 1) {
			open(OUT1, ">>out1.tmp");
			open(OUT2, ">>out2.tmp");
		}
		if ($count % 1000000 == 0) {
			print "processed $count lines\n";
			open(OUT1, ">>out1.tmp");
			open(OUT2, ">>out2.tmp");
			print OUT1 join("\n", @file1),"\n";
			print OUT2 join("\n", @file2),"\n";
			@file1 = ();
			@file2 = ();
		}
		if ($count % 4 == 2) {
			$_ =~ /([ACTGN]+)\s([ACTGN]+)/i;
			$fwd = $1;
			$rev = $2;
		# The 12 in the code below specifies how much is required to recognize the start of fragment. 
		# To save as much info as possible, I copy the barcode!
			$b1 = substr($fwd, 0, 12);
			$b2 = substr($rev, 0, 12);
			$bc1 = reverse($b1);
			$bc2 = reverse($b2);
			$bc1 =~ tr/ACTGactg/TGACtgac/;
			$bc2 =~ tr/ACTGactg/TGACtgac/;
			if ($rev =~ /.+$bc1/ && $fwd =~ /(.+$bc2)/) {
		
				$m1 = $1;
				$rev =~ /(.+$bc1)/;
				$m2 = $1;
				$l1 = length($m1)+length($b1);
				$l2 = length($m2)+length($b2);
				if ($l1 != $l2) {
					++$yikes;
					push @file1, $b1.$fwd;
					$l1 = length($fwd)+length($b1);
					push @file2, $b2.$rev;
					$l2 = length($rev)+length($b2);
				} else {
					push @file1, $b1.$m1;
					push @file2, $b2.$m2;
				}		
			} else {
				push @file1, $b1.$fwd;
				$l1 = length($fwd)+length($b1);
				push @file2, $b2.$rev;
				$l2 = length($rev)+length($b2);
			}
		} elsif ($count % 4 == 0) {
			push @file1, "J" x $l1;
			push @file2, "J" x $l2;
		} else {
			chomp $_;
			@line = split (/\t/,$_);
			push @file1, $line[0];
			push @file2, $line[1];
		}
		if(eof()){	
			print "lost $yikes reads due to discordant length\n"; 
			print OUT1 join("\n", @file1),"\n";
			print OUT2 join("\n", @file2),"\n";
			close(OUT1);close(OUT2);
		}
		' 
	

	echo `date` " [$$] Making unaligned bam-file for: " $fname

	java -jar $picard FastqToSam F1=out1.tmp F2=out2.tmp O=bam_unaligned/$fname.unaligned.bam SM=$fname TMP_DIR=tmp
	rm out*.tmp

	echo `date` " [$$] Doing some really convoluted file splitting for: " $fname
	$samtools view bam_unaligned/$fname.unaligned.bam | awk '{ print > "/bamdump/"length($10)".sam" }'
	rm bam_unaligned/$fname.unaligned.bam

	echo `date` " [$$] Creating DCS fasta-files for: " $fname
	for l in bamdump/*.sam
		do
		bn=`basename $l`
		sn=${bn/.sam/}
		echo `date` " [$$] Creating DCS fasta-file for: " $bn
		$python2 $UCM --input $l --prefix DCSdump/$sn --minmem 2 --spacerlen 0
	done
	zcat DCSdump/*_read1_dcs.fq.gz > DCS_output/$fname"_read1_dcs.fq"
	gzip DCS_output/$fname"_read1_dcs.fq"
	zcat DCSdump/*_read2_dcs.fq.gz > DCS_output/$fname"_read2_dcs.fq"
	gzip DCS_output/$fname"_read2_dcs.fq"
	rm bamdump/*.*
	rm DCSdump/*.*


done


for i in DCS_output/*_read1_dcs.fq.gz
	do

	bn=`basename $i`
	fname=${bn/_read1*/}
        r1="DCS_output/"$fname"_read1_dcs.fq.gz"
        r2="DCS_output/"$fname"_read2_dcs.fq.gz"

	
$bwa aln $ref $r1 > bam/read_1.aln

$bwa aln $ref $r2 > bam/read_2.aln

$bwa sampe -s $ref bam/read_1.aln bam/read_2.aln $r1 $r2 > bam/$fname"_DCS_PE.aln.sam"

$samtools view -Sbu bam/$fname"_DCS_PE.aln.sam" | $samtools sort -o bam/$fname"_DCS_PE.aln.sort.bam"

$samtools index bam/$fname"_DCS_PE.aln.sort.bam"

$samtools view -F 4 -b bam/$fname"_DCS_PE.aln.sort.bam" > bam/$fname"_DCS_PE.filt.bam"

java -jar $picard AddOrReplaceReadGroups INPUT=bam/$fname"_DCS_PE.filt.bam" OUTPUT=bam/$fname"_DCS.bam" RGLB=UW RGPL=Illumina RGPU=ATATAT RGSM=default

$samtools index bam/$fname"_DCS.bam"



rm tmp.vcf
rm bam/read_1.aln
rm bam/read_2.aln

done

rm bam/*PE.aln*

rm bam/*PE.filt.bam


for i in bam/*_DCS.bam
	do
	bn=`basename $i`
	export fname=${bn/_*/}
	echo $fname
	$bcftools mpileup -A -I -q 20 -Q 20 -a AD -f $ref -d 10000 -R $bed $i | perl -ne '
		if ($_ =~ /#CHROM/) {print "Locus\tA\tC\tG\tT\tN";}
		if ($_ =~ /^chr/) {
			chomp $_;
			@line = split (/\t/, $_);
			print "\n" . $line[0] . ":" . $line[1];
			@si = split (/:/, $line[9]);
			@counts = split (/,/, $si[1]);
			@alleles = split (/,/, $line[3].",".$line[4]);
			@pbc = (0, 0, 0, 0, 0);
			foreach $base (0 .. @alleles-1) {
				if ($alleles[$base] eq A) {
					$pbc[0] = $counts[$base];
				} elsif ($alleles[$base] eq C) {
					$pbc[1] = $counts[$base];
				} elsif ($alleles[$base] eq G) {
					$pbc[2] = $counts[$base];
				} elsif ($alleles[$base] eq T) {
					$pbc[3] = $counts[$base];
				} else {
					$pbc[4] = $counts[$base];
				}
			}
			foreach $b (@pbc) {
				print "\t" . $b;
			}
		}
	' > "doc/"$fname"_doc.tsv"
done

# There will also be an Rscript available soon to obtain the mutation load information directly from the pipeline. For now the depth of coverage data has to be analyzed from an R session.


exit


