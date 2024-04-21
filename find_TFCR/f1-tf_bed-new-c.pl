use strict;

my $dirname1="./e_tfpos_combine/";
opendir(DIR,$dirname1)|| die "Can't open directory $dirname1"; 
my @filename1 = readdir(DIR); 
closedir DIR; 
my @runfile=();
my $runi=0;
foreach (@filename1)
{ 
	if(($_ eq ".")||($_ eq "..")){next;}
    #if(!($_ =~ /SRX/)){next;}  
    @runfile[$runi]=$_;
    $runi++;
} 
print "$runi\n";

for(my $nn=0;$nn<$runi;$nn++)
#for(my $nn=325;$nn<$runi;$nn++)
{
	
	my $dirname4="./f_tfbed/";
	if(-e "$dirname4@runfile[$nn]"){next;}
	my $command="./guassian_method/guassian_300 $dirname1@runfile[$nn] $dirname4@runfile[$nn]";
	print "$command\n";
	system("$command");
}