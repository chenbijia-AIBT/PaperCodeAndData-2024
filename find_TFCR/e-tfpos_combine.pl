use strict;

my $nfile1=1; 
my $dirname1="./d_motif_combine/";
#my $dirname0="./../tf_test/b_DNase_Rk/";
opendir(DIR,$dirname1)|| die "Can't open directory $dirname1"; 
my @filename1 = readdir(DIR); 
foreach   (@filename1){ 
         $nfile1++;
} 
closedir DIR; 
@filename1=sort{$a cmp $b}@filename1;


for(my $nn=0;$nn<$nfile1-1;$nn++)
#for(my $nn=82;$nn<$nfile1-1;$nn++)
{
	#if((@filename1[$nn] ne "Cd4naivewb11970640_UW.txt")){next;}
	if((@filename1[$nn] eq ".")||(@filename1[$nn]eq "..")){next;}
	opendir(DIR2,"$dirname1@filename1[$nn]")|| die "Can't open directory"; 
	my @filename2 = readdir(DIR2);
	my $nft=@filename2;
#	if($nft<494){next;}
	print "@filename1[$nn]\n";
	my $dirname4="./e_tfpos_combine/";
#	if(-e "$dirname4@filename1[$nn]"){next;}
	my @data=();
	my $datai=0;
	foreach my $intf(@filename2)
	{
		if($intf eq "." ||$intf eq ".."){next;}

		open IN,"$dirname1@filename1[$nn]/$intf"or die "read error";
		while(my $line=<IN>)
		{
			chomp $line;
			$line=~s/[\n\r]//;
			my @temp=split /\t/,$line;
			#if((@temp[0] ne "chr2")&&(@temp[0] ne "chr3")&&(@temp[0] ne "chr4")&&(@temp[0] ne "chr5")&&(@temp[0] ne "chr6")&&(@temp[0] ne "chr7")&&(@temp[0] ne "chr8")&&(@temp[0] ne "chr9")&&(@temp[0] ne "chrX")){next;}
			my $tpos=int ((@temp[1] + @temp[2])/2);
			@data[$datai]="@temp[0]\t$tpos";
			$datai++;
		}
		close IN;
	}
	print "total number:\t$datai\n";
	@data=sort{(split /\t/,$a)[1] <=> (split /\t/,$b)[1]}@data;
	print "ok1\n";
	@data=sort{(split /\t/,$a)[0] cmp (split /\t/,$b)[0]}@data;
		
	open OUT,">$dirname4@filename1[$nn]"or die "open error";
	my $prechr="";
	my $prePos=-1000;
	my $num=0;
	foreach(@data)
	{
		my @temp=split /\t/,$_;
		if(($prechr ne @temp[0])||(@temp[1] ne $prePos))
		{
			if(!$prechr)
			{
				$prechr = @temp[0];
				$prePos = @temp[1];
				$num++;			
				next;
			}			
			print OUT "$prechr\t$prePos\t$num\n";				
		
			$prechr = @temp[0];
			$prePos = @temp[1];
			$num=0;		
		}
		$num++;
		#$prePos = @temp[1];			
	}
	print OUT "$prechr\t$prePos\t$num\n";


	close OUT;
	
}