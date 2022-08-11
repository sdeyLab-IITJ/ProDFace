#*************************************************************************************************************************************************
#perl program to find out the total No. of hydrogen bonds formed by a single interface (either A or B at a time) through water.
#Apart from HBs, it also calculates no of 1 side  Hbonded interface water molecules, number of interface water that do not hbond to protein, number that do not hbond to pr but form wate-water hbond and number that do not hbond at all

#****************************************************************************************************
#input files   :  1A3Q.hb2, 1A3Q.pdb

#                 "list_file" containing only 4 letter pdb code also all the .pdb.hb2 ,.pdb.water and  .int files
#		          "int_water_jan5_info" o/p of "int_water_jan5.pl"

#output files  :  1A3Q.water, "int_water_info", Water.dat

#                 " No_wa_1side_hb_all.dat", "No_1_side_WHBs_all.dat","No_wa_no_hb_pr.dat" ,"No_wa_no_hb_pr_but_wa_wa.dat","No_wa_no_hb_atall.dat"

#NOTE :  IT TAKES INTO ACCOUNT ALL HBONDS CONTRIBUTED BY INTERFACE WATERS ,NOT NECESSARILY INVOLVING INTERFACE ATOMS/RESIDUES ONLY.

# Usage 1A3Q
    # Authors : Sucharita Dey, Arumay PAL                                          |
    # Contact : sdey@iitj.ac.in; contact.arumay@gmail.com
#*************************************************************************************************************************************************

#!/usr/bin/perl -w

$line=$ARGV[0];  

#foreach $line(@da)
#{   
  	    chomp($line);
    	print $line,"\n";
	
  	#$hb2=$line.".pdb.hb2";
	$hb2 = "./" . $ARGV[0] . ".hb2";

	#$fwat=$line.".pdb.water";
	$fwat = "./" . $ARGV[0] . ".pdb";

	$output = "./" . $ARGV[0] . ".water";
	
	$fo_intWaterInfo = "./" . "int_water_info";

	$fo_water_dat = "./" . "Water.dat";

	
  	@whb=();
	
  	@cw=();$c3=0;$wa3=0;
	  	 	
	$wa1side=0;$wa_wa=0;$wa_no_hbatall=0;
  		
	$flag=0;
  
	unless (open(HB,$hb2))
	{
		print "file is not open.";
  	}
 		@filehb=<HB>;
  		close HB;
  
  	open(OUT2,">>$output");
#print OUT2 "****************************************************************\n";
print OUT2 "Interaction of rest of the interface waters that hydrogen bond with either protein or DNA\n"; 	
print OUT2 "-----------------------------------------------------------------------------------------\n"; 	
print OUT2 "With Protein:\n"; 	
 
 # Reading "int_water_jan5_info" (o/p of int_water.pl) and collecting interface water number and ch type.
  
  	
	open(FH4,$fo_intWaterInfo);
	while(<FH4>)
	{
		#if($_=~/^$fwat.*HOH.(....)(.).*/)
		if($_=~/^$fwat.*(.):.*HOH.(....)(.).*(.):.*/)
		{
       			push(@ch11,$1);push(@ch22,$4);
       			$cwtype=$3;
       			$wnumber=$2;
       			$wnumber=~s/ /0/g;

       			if($cwtype eq " ")
       			{
       				$cwtype=~s/ /-/g;
       			}	
       			$cw[$c3]=$cwtype.$wnumber."\-"."HOH"." "."O  ";
       			#print "$cw[$c3]\n";
			$c3++;
    		}
  	}
	close FH4;
  #print "Total int water $line -- $c3\n";
 
# Making a uniq set of chain IDs

	@ch11_sort= sort(@ch11); @ch22_sort= sort(@ch22);
while($i<@ch11_sort)
{
	while($ch11_sort[$i] eq $ch11_sort[$i+1])
	{
		$i++;
	}
	push(@ch11_fi,$ch11_sort[$i]);$i++;
}
$i=0;
while($i<@ch22_sort)
{
	while($ch22_sort[$i] eq $ch22_sort[$i+1])
	{
		$i++;
	}
	push(@ch22_fi,$ch22_sort[$i]);$i++;
}


 
 #Taking each interface water molecule, its hbonding partners are searched from the .hb2 file,  
  	for($k=0;$k<$c3;$k++)
  	{
  
          #print "$cw[$k]\n";
	  	
          $flag=0;
	  $temp_res1="";
	  $temp_res2="";
		 foreach $l(@filehb)
	  	{
	    #A0010-ARG N   A0006-ILE O   3.17 MH  -2 -1.00  -1.0 -1.00  -1.0 155.7     3  && $f ne '-1'
	    #A0008-LYS N   A0004-ALA O   2.87 MM  -1  6.00 168.7  1.89 169.0 165.6     1
	    #A0004-LYS NZ  Z0878-HOH O   3.35 SH  -2 -1.00 140.1  2.49  -1.0  -1.0     6
	


	   		if($l=~/^(.............).(.............)......(..)..........(..)\s*.*/)
	   		{
	     			$at1=$1;
	     			$at2=$2;
				$chch=$3;
				if(($at1=~/HOH/) && ($at2=~/HOH/))
				{
					push(@wa_only,$at1);
					push(@wa_only,$at2);
					goto LABEL1;
				}	
				
				if($cw[$k] eq $at1)
				{
					$temp_res1=$at2;
					if($at2 ne $temp_res2)
					{
						push(@whb,$at2);
						push(@chchtype,$chch);
					}
				}
				if($cw[$k] eq $at2)
				{
					$temp_res2=$at1;
					if($at1 ne $temp_res1)
					{
						push(@whb,$at1);
						push(@chchtype,$chch);
					}
				}	
				
			}
			LABEL1:
	    	}
		$size_whb=@whb;
		#print "$size_whb\n";
		if($size_whb eq 0)
		{
			push(@nobp,$cw[$k]);
			$wa3++;
			$flag =1;
		}	
	   #print "$whb[0]\n";
	   $test1="XXXX";
	   $test2="XXXX";
	    	for($ii=0;$ii<$size_whb;$ii++)
	    	{
	    		for($jj=$ii+1;$jj<$size_whb;$jj++)
			{
				if($whb[$ii] ne $whb[$jj]) 
				{
					for($ch1i=0;$ch1i<@ch11_fi;$ch1i++)
					{
						for($ch1j=0;$ch1j<@ch22_fi;$ch1j++)
						{
							$first_i=substr($whb[$ii],0,1);
							$first_j=substr($whb[$jj],0,1);

					
							#if( (($whb[$ii]=~/^A/) && ($whb[$jj]=~/^B/)) || (($whb[$ii]=~/^B/) && ($whb[$jj]=~/^A/)) )
							if( (($first_i eq $ch11_fi[$ch1i]) && ($first_j eq $ch22_fi[$ch1j])) || (($first_i eq $ch22_fi[$ch1j]) && ($first_j eq $ch11_fi[$ch1i])) )
							{
								push(@res1,$whb[$ii]);
								push(@res2,$whb[$jj]);
								push(@chch1,$chchtype[$ii]);
								push(@chch2,$chchtype[$jj]);
							}
						}
					}
				}
			        
			}
	    	}
	    	

		$size_res1=@res1;
# Calculating 1 sided Whbs in interface, also no. of interface waters with 1 side hbond		
		if($size_res1 eq 0 && $flag ne 1)
		{
			push(@wa_1sided,$cw[$k]);
			$wa1side++;
			
	# Separating the Hbond partners (protein and DNA) of each interface HOH 
			$prot_flag=0;
			foreach(@whb)
			{
				if($_=~/.*-\w\w\w .*/)
				{
					push(@whbProt,$_);
					$prot_flag=1;
					
					
				}
				else
				{
					push(@whbDna,$_);
					
				}
			}
			$hb1prot+=@whbProt;
			if($prot_flag eq 1)
			{
				$pp=join(' ',@whbProt);
				$temp=$cw[$k]."--".$pp;
				push(@protpair,$temp);

				$wa1prot++; 
			}
			else
			{
				$pp=join(' ',@whbDna);
				$temp=$cw[$k]."--".$pp;
				push(@dnapair,$temp);
			}

			#print OUT2 "1side pairs : $cw[$k]--@whb\n";
			$hb1side+=@whb;
			
			
		}	
		

		@whbDna=();
		@whbProt=();
		@whb=();
		@res1=();
		@res2=();
		@chch1=();
		@chch2=();
		@chchtype=();
		
		
			
	}
	@filehb=();
	@fint=();
	
	$ww2_size=@wa_only;

	print " no. of WATERS with 1 side hbonded  -- $wa1side\n";
	print " no. of  1 side HBONDS -- $hb1side\n";
	
	foreach $ww1(@nobp)
	{
		$sucha=0;
		for($ww2=0;$ww2<$ww2_size;$ww2++)
		{
			if($ww1 eq $wa_only[$ww2])
			{
				$sucha=1;
			}
		}
		if($sucha eq 0)
		{
			$wa_no_hbatall++;
		}
		if($sucha eq 1)
		{
				$wa_wa++;
		}
	}	


	open(OUT1,">>$fo_water_dat");

	printf(OUT1 "#Water mediated HBs bonded to only one subunit:%5d (protein %3d, DNA %3d)\n",$hb1side,$hb1prot,$hb1side-$hb1prot);

	printf(OUT1 "#Waters H-bonded to only one subunit          :%5d (protein %3d, DNA %3d)\n",$wa1side,$wa1prot,$wa1side-$wa1prot);
		
	printf(OUT1 "#Waters only H-bonded to other waters         :%5d \n",$wa_wa);


	foreach(@protpair)
	{
		print OUT2 "1side pairs : $_ \n";
	}
	#print OUT2 "--------------------------------------------------\n";
	print OUT2 "With DNA:\n"; 	

	foreach(@dnapair)
	{
		print OUT2 "1side pairs : $_ \n";
	}


	#open(OUT1,">>No_wa_1side_hb_all.dat");
	#printf(OUT1 "%s      %5d \n", $line,$wa1side );
	
	#open(OUT2,">>No_1_side_WHBs_all.dat");
	#printf(OUT2 "%s      %5d \n", $line,$hb1side );
	
	
	#open(OUT3,">>No_wa_no_hb_pr.dat");
	#printf(OUT3 "%s      %5d \n",$line,$wa3);
	
	
	#open(OUT4,">>No_wa_no_hb_atall.dat");
	#printf(OUT4 "%s      %5d \n",$line,$wa_no_hbatall);
	
	#open(OUT5,">>No_wa_no_hb_pr_but_wa_wa.dat");
	#printf(OUT5 "%s      %5d \n",$line,$wa_wa);
	
	#print " $line   wa_no_hb -- @nobp \n"; 
	#print " no. of wa with no bp    -- $wa3 \n\n";

	$hb1side=0;
	$wa1side=0;
	$wa_no_hbatall=0;
	$wa_wa=0;
	@wa_1sided =();
	@nobp=();
	@wa_only=();
	$wa1prot=0;$hb1prot=0;
	@protpair=();@dnapair=();
	
#}		
close OUT1;
close OUT2;

