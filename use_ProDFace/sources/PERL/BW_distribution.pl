#*************************************************************************************************************************************************
# For server
#This program gives a distribution of Bridging water(BW).It gives the total number of BW, and a distribution of total water mediated Hbonds.
#
# 

#input file   :   1234_org.hb2 ,int_water_info (o/p of int_water.pl ,containing interface waters)
#output files :  >> "Water.dat", "1A3Q.water"
#                      

#NOTE: A particular water may be involved in both SC-SC and SC-MC or MC-MC, but its contribution is counted only for SC-SC. In case of no SC-SC, considered for SC-MC and so on 

#USAGE:    

    # Authors : Sucharita Dey, Arumay PAL                                          |
    # Contact : sdey@iitj.ac.in; contact.arumay@gmail.com

#NOTE : IT TAKES INTO ACCOUNT ALL HBONDS CONTRIBUTED BY INTERFACE WATERS TO BOTH SUBUNITS, NOT NECESSARILY INVOLVING INTERFACE ATOMS/RESIDUES ONLY..
#*************************************************************************************************************************************************

 
#!/usr/bin/perl -w

#$file="list_file";

#unless (open(FH1,$file))
#{
#	print "file is not open.";
#}
#@da=<FH1>;
#close FH1;

$name1 = "./" . "int_water_info";

unless (open(FH2,$name1))
{
	print "file is not open.";
}
@int_wa=<FH2>;
close FH2;

$name2 = "./" . "Water.dat";

open(OUT1,">>$name2");

$name3 = "./" . "1A3Q.water";
open(OUT2,">>$name3");
print OUT2 "****************************************************************\n";
print OUT2 "H-bond information of bridging water at protein-DNA interface\n"; 	


#foreach $line(@da)
#{
  	
   
  	#chomp($line);
    	#print $line,"\n";
	
  	#$hb2=$line.".pdb.hb2";

	$hb2 = "./" . "1A3Q.hb2";

	
  	@ch1=();
	
  	@whb=();
	
  	@cw=();
	
  	$c1=0;
  	$c2=0;
  	$c3=0;
  	
  	$count=0;
	$no_br_water=0;
	$tot_Whbs=0;
	
	
	
  	$MM=0;$MS=0;$SM=0;$SS=0;$SB=0;
	$MM_BW=0;$MS_BW=0;$SM_BW=0;$SS_BW=0;$total_int_SB=0;
	$bw_ss=0;$bw_ms=0;$bw_mm=0;$bw_sb=0;
  
  
	unless (open(HB,$hb2))
	{
		print "file is not open.";
  	}
 		@filehb=<HB>;
  		close HB; 	
  	  
 # Collecting water number and ch type.
  
  	foreach $went(@int_wa)
	{
    #HETATM 5187  O   HOH Z   3      31.668  53.740   4.451  1.00 24.05           O        18.64     2.64
    
    		if($went=~/^$line.*(.):.*HOH.(....)(.).*(.):.*/)
    		{
				push(@ch11,$1);push(@ch22,$4);

				$cwtype=$3;
       			$wnumber=$2;
       			$wnumber=~s/ /0/g;
       
       			if($cwtype eq " ")
       			{
       				$cwtype=~s/ /-/g;
       			}	
       			$cw[$c3]=$cwtype.$wnumber."-"."HOH"." O  ";
       			
			$c3++;
    		}
  	}
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

  
  
 #Taking each water molecule, its hbonding partners are searched from the .hb2 file,  
  	for($k=0;$k<$c3;$k++)
  	{
  
          #print "$cw[$k]\n";	
          
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
		
	   $test1="XXXX";
	   $test2="XXXX";
	   $test_cw="XXXXX";

#Separating the bonding partners from  subunits A and B 
	   
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
		

# Removing redundancy of the bonding partners (from ch A & B) of the water molecule and counting  chainwise (ss,sm/ms,mm)		
	    	
		for($kk=0;$kk<$size_res1;$kk++)
		{
				if(($res1[$kk] ne $test1) || ($res2[$kk] ne $test2))
				{
					
				print OUT2 $res1[$kk]," : ",$cw[$k]," : ",$res2[$kk],"\n";
				
					$count++;
					if( (($chch1[$kk] eq "HS") || ($chch1[$kk] eq "SH")) && (($chch2[$kk] eq "HS") || ($chch2[$kk] eq "SH")) ) 
					{ 
						$SS++;$SS_BW++;
						#print $line,"--",$res1[$kk]," ",$cw[$k]," ",$res2[$kk],"\n";
						if( ((($res1[$kk] =~/ASP/) ||($res1[$kk] =~/GLU/)) && (($res2[$kk] =~/LYS/) || ($res2[$kk] =~/ARG/) || ($res2[$kk] =~/HIS/))) || ((($res2[$kk] =~/ASP/) ||($res2[$kk] =~/GLU/)) && (($res1[$kk] =~/LYS/) ||($res1[$kk] =~/ARG/) || ($res1[$kk] =~/HIS/))) )
						{
							$total_int_SB++;$SB++;
							#print $line,"--",$res1[$kk]," ",$cw[$k]," ",$res2[$kk],"\n";
						}	
					}
					if( (($chch1[$kk] eq "HM") || ($chch1[$kk] eq "MH")) && (($chch2[$kk] eq "HS") || ($chch2[$kk] eq "SH")) ) { $MS++;$MS_BW++;}
					if( (($chch1[$kk] eq "SH") || ($chch1[$kk] eq "HS")) && (($chch2[$kk] eq "HM") || ($chch2[$kk] eq "MH")) ) { $MS++;$MS_BW++;}
					if( (($chch1[$kk] eq "MH") || ($chch1[$kk] eq "HM")) && (($chch2[$kk] eq "HM") || ($chch2[$kk] eq "MH")) ) { $MM++;$MM_BW++;}
				
					$test1=$res1[$kk];
					$test2=$res2[$kk];
					if($cw[$k] ne $test_cw) 
					{ 	
						$no_br_water++; 
						$test_cw=$cw[$k];
					}
				
				}
					
		
	     	}
		#print "no. of hbs of $cw[$k] -- $count\n";
		$tot_Whbs+=$count;
		
		$count=0;
      		@whb=();
		@res1=();
		@res2=();
		@chch1=();
		@chch2=();
		@chchtype=();
		if($SS_BW > 0)
		{
			$bw_ss++;
			if($total_int_SB > 0)
			{
			 	$bw_sb++;
			}	
		}
		elsif(($SS_BW ==0) && ($MS_BW > 0))
		{
			$bw_ms++;
		}
		elsif(($SS_BW ==0) && ($MS_BW == 0) && ($MM_BW >0))
		{
			$bw_mm++;
		}	
		
		
		$SS_BW="0";$MS_BW="0";$MM_BW="0";$total_int_SB=0;	
				
	}
	@filehb=();
	@fint=();
	@wt=();
	
	
	printf(OUT1 "#Bridging Waters   :%5d\n",$no_br_water);
		
	printf(OUT1 "#DNA-HOH-Protein interactions:%5d\n",$tot_Whbs);
	 
	$#ch1=-1;
	$#cw=-1; 
#} 
print OUT2 "****************************************************************\n";
    
close OUT1;
close OUT2;
		
	    				
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
			  
	  
