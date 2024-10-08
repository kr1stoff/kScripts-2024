#!/usr/bin/perl -w
if(@ARGV<1)
{
   print "perl $0 <fasta_seq>[cut_off_len]\n";
   exit;
}

@len=();$t_len=0;$n_100=0;$n_2000=0;
$cut_off=$ARGV[1] || 100;

$len=0;
open IN,"$ARGV[0]";
while($line=<IN>)
{
   chomp($line);
   if($line=~/^>/)
   {
      if($len>=100){$n_100++;}
      if($len>=2000){$n_2000++;}
      if($len>=$cut_off)
      {
         push @len,$len;
         $t_len+=$len;
      }
      $len=0;
   }
   else
   {
      $len+=length($line);
   }
}

if($len>=100){$n_100++;}
if($len>=2000){$n_2000++;}
if($len>=$cut_off)
{
   push @len,$len;
   $t_len+=$len;
}
close IN;

@len = sort {$b <=> $a} @len;
#print "Total_Len(bp)\tN50_Len(bp)\tN90_Len(bp)\tMax_Len(bp)\tNumber>=100bp\tNumber>=2000bp\n";
print "序列总长度(bp)\tN50(bp)\tN90(bp)\t最长的Contig长度(bp)\t>=100bp序列数\t>=2000bp序列数\n";
#print "Max_length\t$len[0]\n";

$nn=10;$sum=0;%n50=();%l50=();
for($i=0;$i<=$#len;$i++)
{
   $sum+=$len[$i];
   if($sum>=$t_len*$nn/100)
   {
      $n50{$nn}=$len[$i];
      $l50{$nn}=$i+1;
      $nn+=10;  
   }
   last if($nn == 100);
}

#print "N90\t$n50{90}\t$l50{90}\n";
#print "N80\t$n50{80}\t$l50{80}\n";
#print "N70\t$n50{70}\t$l50{70}\n";
#print "N60\t$n50{60}\t$l50{60}\n";
#print "N50\t$n50{50}\t$l50{50}\n";
#print "N40\t$n50{40}\t$l50{40}\n";
#print "N30\t$n50{30}\t$l50{30}\n";
#print "N20\t$n50{20}\t$l50{20}\n";
#print "N10\t$n50{10}\t$l50{10}\n";
#print "N90\t$n50{90}\n";
#print "N50\t$n50{50}\n";
#print "Total_length\t$t_len\n";
#print "number>=100bp\t$n_100\n";
#print "number>=2000bp\t$n_2000\n";
print "$t_len\t$n50{50}\t$n50{90}\t$len[0]\t$n_100\t$n_2000\n";
