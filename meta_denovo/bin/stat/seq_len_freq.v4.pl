#!/usr/bin/perl
#-------------------------------------------------+
#    [APM] This script was generated by amp.pl    |
#    [APM] Created time: 2015-10-21 13:45:34      |
#-------------------------------------------------+
# name: fasta_len_freq.pl
# func: stat the fasta length frequency and some characteristic values
# version: version 4.0
# update: 2017-10-26
#-------------------------------------------------+

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Std;
use File::Basename qw/basename dirname/;
use FindBin qw/$Bin/;

use List::Util qw/sum max/;


# default value is for Meta Genomics, stat the scaftigs length frequency
my %opts = (s=>0,e=>3000,n=>30,t=>"contigs");
getopts('s:e:n:w:t:k:o:f',\%opts);

&usage unless @ARGV == 1;

my $Rscript = "Rscript";

my $file = shift @ARGV;
my $type = $opts{t};

my $sample = $opts{k} ? $opts{k} : basename($file,".$opts{t}.fa");
my $odir   = $opts{o} ? $opts{o} : dirname($file);

# read the sequence file and calculate the sequence length
my $inseq = Bio::SeqIO->new(-file=>$file,-format=>"fasta");
my @lengths;

open my $ofh_len    , ">" , "$odir/$sample.$opts{t}.length" or die $!;
open my $ofh_contig , ">" , "$odir/$sample.$opts{t}.filter.fa" or die $! if ($opts{f});

my $order = "00000000";

while (my$seq = $inseq->next_seq)
{
	my $id = $seq->id;
	my $length = $seq->length;
	my $seqstr = $seq->seq;

	# remove the contig whose length < minmum allowed length
	next if ($length < $opts{s});
	push @lengths , $length;
	$order ++;

	my $trans_len = $length > $opts{e} ? $opts{e} : $length;

	my $seqid = "${sample}_contig_$order";

	print $ofh_len "$seqid\t$length\t$trans_len\n";
	print $ofh_contig ">$seqid\n$seqstr\n" if ($opts{f});
}
close $ofh_contig if ($opts{f});
close $ofh_len;

# do some stat
my $total     = scalar @lengths;
my $total_len = sum(@lengths);
my $max_len   = max(@lengths);
my $avg_len   = sprintf "%0.2f" , $total_len/$total;
my ($n50,$n90)= fetch_N_size(\@lengths,50,90);

open my $ofh_log , ">" , "$odir/$sample.$opts{t}.stat.txt" or die $!;
print $ofh_log <<STAT;
SampleID: $sample
Num: $total
Total len (bp): $total_len
Average len (bp): $avg_len
Max len (bp): $max_len
N50 len (bp): $n50
N90 len (bp): $n90
STAT

close $ofh_log;

# draw frequency figure
my $step = $opts{w} || int(($opts{e}-$opts{s})/$opts{n});

my ($c_hist,$c_line) = ("skyblue","orangeRed");
my $R = <<R;
data = read.table("$odir/$sample.$opts{t}.length",header=F)
cell<-seq(from=$opts{s},to=$opts{e}+$step,by=$step)
nums=data[,3]
r=hist(nums,breaks=cell)
per=c(0,r\$counts/$total*100)
mids=c(0,r\$mids)
ymax=max(r\$counts)*1.1
mycolors <- colorRampPalette(c(\"$c_hist\", \"$c_line\"))(2)
pdf("$odir/$sample.$opts{t}.length_distribution.pdf",width = 8,height = 6)
par(mar=c(5,4,4,8)+0.1)
hist(nums,breaks=cell,xlim=c(0,$opts{e}),ylim=c(0,ymax),col=mycolors[1],ann=FALSE)
legend(x=$opts{e}*0.6,y=ymax*0.90,"Frequence(#)",fill=mycolors[1],bty='n')
legend(x=$opts{e}*0.6,y=ymax*0.98,"Percentage(%)",lty=1,col=mycolors[2],lwd=2,bty='n',seg.len = .8)
par(new=TRUE)
plot(mids,per,type='l',xlim=c(0,$opts{e}),ylim=c(0,max(per)*1.5),axes=F,ann=FALSE,col=mycolors[2],lwd=2)
axis(4,col="black",col.ticks="black",col.axis="black")
mtext("Frequence(#)",side=2,line=2)
mtext("Percentage(%)",side=4,line=2)
title("$type Length Distribution of $sample",xlab="$type length(bp)",line=2)
dev.off()
png("$odir/$sample.$opts{t}.length_distribution.png",width =800,height = 600)
par(mar=c(5,4,4,8)+0.1)
hist(nums,breaks=cell,xlim=c(0,$opts{e}),ylim=c(0,ymax),col=mycolors[1],ann=FALSE)
legend(x=$opts{e}*0.6,y=ymax*0.90,"Frequence(#)",fill=mycolors[1],bty='n')
legend(x=$opts{e}*0.6,y=ymax*0.98,"Percentage(%)",lty=1,col=mycolors[2],lwd=2,bty='n',seg.len = .8)
par(new=TRUE)
plot(mids,per,type='l',xlim=c(0,$opts{e}),ylim=c(0,max(per)*1.5),axes=F,ann=FALSE,col=mycolors[2],lwd=2)
axis(4,col="black",col.ticks="black",col.axis="black")
mtext("Frequence(#)",side=2,line=2)
mtext("Percentage(%)",side=4,line=2)
title("$type Length Distribution of $sample",xlab="$type length(bp)",line=2)
dev.off()
R

open R,">$odir/draw.R" or die $!;
print R $R;
close R;
system("$Rscript $odir/draw.R 1>/dev/null");
system("rm $odir/Rplots.pdf") if (-e "$odir/Rplots.pdf");

sub usage
{
	print <<HELP;
Usage:   perl $0 [options] <*.fatsa>

Options: -s INT    the minmum length to be stat, [0]
         -e INT    the maxmum length to be stat, [3000]
         -t STR    the sequence type, this name will be diaplayed in the figure [contigs]
         -w INT    set the stat windows, default use the -n, is prior to -n
         -n INT    set the sum steps, [30]
         -k STR    set the keyname of output file, default fetch from the filename
         -o STR    set the outpur dir, default use the dir which contains the input file
HELP
	exit;
}

sub fetch_N_size
{
	my $arr = shift;
	my @locis = @_;
	my @locis_sort = sort {$a<=>$b} @locis;

	my $sum = sum(@$arr);
	my $accum_val = 0;
	my %nsizes;
	OUTER:foreach my $val (sort {$b<=>$a} @$arr)
	{
		$accum_val += $val;
		INNER:foreach my $loci (@locis_sort)
		{
			if ($accum_val < $loci * $sum/100)
			{
				next OUTER;
			}
			else
			{
				if($nsizes{$loci})
				{
					next;
				}
				else
				{
					$nsizes{$loci} = $val;
				}
			}
		}
	}

	my @nsizes = map { $nsizes{$_} } @locis;
	return @nsizes;
}
