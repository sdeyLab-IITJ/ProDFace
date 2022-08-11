#*******************************************************************************************************
#For Server
# Program to calculate Number of Water molecules at the interface
# Distances of HOH(s) were calculated only from .int atoms(considering only those atoms which are present in the .int file)
# Distances has to be <=4.5 from both subunits

#  I/P files: 1a3Qc1.int, 1A3Qc2.int , 1A3Q.pdb (pdb files containing the water coordinates)
#  O/P files: Water.dat, int_water_info(intermediate file) 

#usage : 1A3Q

    # Authors : Sucharita Dey, Arumay PAL                                          |
    # Contact : sdey@iitj.ac.in; contact.arumay@gmail.com
#*********************************************************************************************************
#!/usr/bin/perl
		
			$name=$ARGV[0];
			$name1 = "./" .$name."c1.int";
			push(@list_Aint,$name1);
			$name2 = "./" .$name."c2.int";
			push(@list_Bint,$name2);
				
			$name3 = "./" .$name.".pdb";
			push(@list_water,$name3);

			$out_file_1 = "./" ."int_water_info";		
			$out_file_2 = "./" ."Water.dat";	

	  		$size_water=@list_water; 

# Work starting for each entry, one by one: 

	for($i=0;$i<$size_water;$i++)
	{
	 	#Reading a.int file to collect coordinates of interface atoms
		
		open(FH_Aint,"./$list_Aint[$i]");

		 while($a=<FH_Aint>)
	 	{
			if($a=~/^ATOM.........(...).(\w\w\w).(.)(....).*/)
			{
		 	push(@atom1,$1);
			push(@res1,$2);
			push(@ch1,$3);
			$res_temp=int($4);
			push(@resno1,$res_temp);
			push(@c1,$a);
						
		 	}
		 
	 	}
	 	close(FH_Aint);
		
		#Reading b.int file to collect coordinates of interface atoms
		
		open(FH_Bint,"./$list_Bint[$i]");
	 
		 while($a=<FH_Bint>)
	 	{
			if($a=~/^ATOM.........(...).(..\w).(.)(....).*/)
			{

		 	push(@atom2,$1);
			push(@res2,$2);
			push(@ch2,$3);
			$res_temp=int($4);
			push(@resno2,$res_temp);
			push(@c2,$a);
						
		 	}
		 
	 	}
	 	close(FH_Bint);
	
		$res_size1=@res1;
		$n1=0;
	
		#Collecting the co-ordinates from a.int file
	
		foreach(@c1)
		{
			$co1[$n1][0]=substr($_,30,8); 
			$co1[$n1][1]=substr($_,38,8);
			$co1[$n1][2]=substr($_,46,8);
			$n1++;
		}

		
		$n2=0;
	
		#Collecting the co-ordinates from b.int file
	
		foreach(@c2)
		{
			$co2[$n2][0]=substr($_,30,8);
			$co2[$n2][1]=substr($_,38,8);
			$co2[$n2][2]=substr($_,46,8);
			$n2++;
		}
		
		#Reading .water file to collect	coordinates of HOH	
		
		open(FH_water,"./$list_water[$i]");
		while($d=<FH_water>)
		{
		 	if($d=~/^HETATM.......O...HOH.(.)(....).*/)
			{
				push(@c_wa,$d);
				push(@wa_ch,$1);
				push(@waterno,$2);
			}
		}
		close(FH_water);
		
		
		#Collecting the co-ordinates (from .water files)
		$o=0;
		
		foreach(@c_wa)
		{
			$co_wa[$o][0]=substr($_,30,8);
			$co_wa[$o][1]=substr($_,38,8);
			$co_wa[$o][2]=substr($_,46,8);
			$o++;
		}
		
		# Calculating Distances...
		
		for($k=0;$k<$o;$k++)
		{
			$flag1=0; $flag2=0;
			
			for($j=0;$j<$n1;$j++)
			{
				$x=$co_wa[$k][0]-$co1[$j][0];
				$y=$co_wa[$k][1]-$co1[$j][1];
				$z=$co_wa[$k][2]-$co1[$j][2];
				$dist1=sqrt($x*$x+$y*$y+$z*$z);
				
				if(($dist1 == 4.5) || ($dist1 < 4.5))
				{
					$flag1=1;$j1=$j;
					last;
				}
			}
			for($j=0;$j<$n2;$j++)
			{		
				$x=$co_wa[$k][0]-$co2[$j][0];
				$y=$co_wa[$k][1]-$co2[$j][1];
				$z=$co_wa[$k][2]-$co2[$j][2];
				$dist2=sqrt($x*$x+$y*$y+$z*$z);
				
				if(($dist2 == 4.5) || ($dist2 < 4.5)) 
				{
					$flag2=1;$j2=$j;
					last;
					
				}
			
			}

		
			if($flag1 ==1 && $flag2==1)
			{
				$no_int_wa++;
				#print OUT2 $list_water[$i]," "," $resno1[$j1] $res1[$j1] $atom1[$j1] $ch1[$j1]:$dist1-- HOH $waterno[$k]$wa_ch[$k] --$resno2[$j2] $res2[$j2] $atom2[$j2] $ch2[$j2]:$dist2\n";
				
				$out_final= $list_water[$i] . " " . " $resno1[$j1] $res1[$j1] $atom1[$j1] $ch1[$j1]:$dist1-- HOH $waterno[$k]$wa_ch[$k] --$resno2[$j2] $res2[$j2] $atom2[$j2] $ch2[$j2]:$dist2\n";
				push(@array_final,$out_final);

			}		
		}
		
		open(OUT2,">$out_file_1");
		print OUT2 @array_final;
		
		
		# Writing in output file
		
		open(fp_out,">$out_file_2");
		printf(fp_out "#interface water   :%5d\n",$no_int_wa);
		
		
		close(fp_out);	
		
		# Refreshing Arrays	
		
		$#c1=-1;	
		$#co1=-1;
		$#atom1=-1;
		$#res1=-1;
		$#ch1=-1;
		$#resno1=-1;
		
		$#c2=-1;	
		$#co2=-1;
		$#atom2=-1;
		$#res2=-1;
		$#ch2=-1;
		$#resno2=-1;
		
		$#c_wa=-1;
		$#co_wa=-1;
		$#wa_ch=-1;
		$#waterno=-1;
		$no_int_wa=0;
		
		$dist1=0;
		$dist2=0;
		
	}			
						
				
