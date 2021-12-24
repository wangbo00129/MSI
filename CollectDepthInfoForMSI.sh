path_bam=$1
path_output=$2

perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /TTTGA(T+)GAGAA/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr4:55598012-55598436 )  > Bat25.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /TTCTC(A+)TCAAA/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr4:55598012-55598436 )  > Bat25.R 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /CAGGT(A+)GGGTT/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:47641360-47641786 )  > Bat26.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /AACCC(T+)ACCTG/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:47641360-47641786 )  > Bat26.R 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /CAGGA(T+)GAGGC/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:39536490-39536916 )  > Mono27.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /GCCTC(A+)TCCTG/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:39536490-39536916 )  > Mono27.R 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /TTGCT(A+)GGCCA/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr14:23652147-23652567 )  > NR21.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /TGGCC(T+)AGCAA/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr14:23652147-23652567 )  > NR21.R 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /TCCTA(T+)GTGAG/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:95849162-95849585 )  > NR24.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /CTCAC(A+)TAGGA/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr2:95849162-95849585 )  > NR24.R 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /CTGGT(A+)GCCAC/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr11:102193309-102193734 )  > NR27.F 
perl -e 'while(my $tag=<>){chomp($tag);my ($rid,$seq)=(split/\s+/,$tag)[0,9];if($seq =~ /GTGGC(T+)ACCAG/){ $len=length($1);print "$rid\t$len\n";} }'  <(samtools view $path_bam  chr11:102193309-102193734 )  > NR27.R

for i in Bat25 Bat26 Mono27 NR21 NR24 NR27
do
	cat ${i}.*|cut -f2 |sort |uniq -c |sort -rn |awk -v OFS="\t" -v tag=${i} '{print tag,$2,$1}'
	rm ${i}.*
done > $path_output
