#*************************************************************************************************************************************************
# This prog caclculates accesibility of the donor groups from all AA

# -------------------------------------------------------------

# I/P :: 1A3Qc1.int       O/P:: 1A3Q_HB_donor_asa_atom.out
# Usage  1A3Q 10.0 -> last agrument is the ASA cut off for the HB donor atom

    # Authors : Sucharita Dey, Arumay PAL                                          |
    # Contact : sdey@iitj.ac.in; contact.arumay@gmail.com
#*************************************************************************************************************************************************

$path = "./";

$file1 = $path . $ARGV[0] . "c1.int";
open(FH1,"$file1");

$file_out = $path . $ARGV[0] . "_HB_donor_asa_atom.hbd";
open(OUT,">$file_out");

$cut_off = $ARGV[1];

@fh1=<FH1>;
close FH1;

@data_all = ();
foreach $line (@fh1)
{
	###### For R,K,H residues	
	
	if($line =~ /^ATOM.........(N\w\d).(ARG).(.)(....)/)
	{	
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}		
	}
	if($line =~ /^ATOM.........(NE.).(ARG).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(N\w\d).(HIS).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(NZ.).(LYS).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}

	###### For other residues	

	if($line =~ /^ATOM.........(NE2).(GLN).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	
	if($line =~ /^ATOM.........(ND2).(ASN).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(OG ).(SER).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(OG1).(THR).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(OH.).(TYR).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}
	if($line =~ /^ATOM.........(NE1).(TRP).(.)(....)/)
	{
		$asa=substr($line,67,5);
		$atom=$1; $res=$2; $ch_int=$3; $r_no=$4;
		$data = $atom. " " . $res . " " . $ch_int . " " . $r_no . " " . $asa;
		if($asa > $cut_off)
		{
		push(@data_all,$data);
		}	
	}

}
	$donor_count = @data_all;
	printf(OUT "Number of Donor:%4d\n",$donor_count);
	for($i=0; $i<@data_all; $i++)
	{
		print OUT "$data_all[$i]\n";
	}
			


