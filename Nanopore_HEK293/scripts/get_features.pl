
while(<>)
{
	chomp;
	my @tmp = split(/\t+/,$_);
	my $offset=$#tmp-15;

	my ($contig, $ref_position, $callScore, $Indel, $base, $site, $strand) = ($tmp[0],$tmp[1]."_".$tmp[2],$tmp[10],0,$tmp[$offset+15], $tmp[7], $tmp[11]);
	my($deletionScore, $insertionScore) = (0,0);
	if(/deletion\_score\=([\d\.E-]+)/)
	{ $deletionScore=$1;}

	if(/insertion\_score\=([\d\.E-]+)/)
	{ $insertionScore=$1;}

	$site=$site-$tmp[1]+1 if($strand eq '+');
        $site=$tmp[2]-$site if($strand eq '-');
							  
	printf("%s\t%s\t%s\t%.2f\t%.2f\t%.2f\t%s\t%d\t%s\n",$contig.":".$ref_position.":".$strand, $contig, $ref_position, $callScore, $deletionScore, $insertionScore, $base, $site, $strand);
}

