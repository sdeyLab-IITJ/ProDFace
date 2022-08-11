/* **********************************************************************************************************
// Conservation of residues across protein-protein interfaces

//List of input files: .pdb, .int, .hssp are the 3 files required as input
//Usage: Give 4-letter PDBcode as the second argument, name of the executable is of course the first argument, percentage Id being the third.
// Ensure that all file names are in lowercase if the PDBcode is provided in lowercase. Often, the pdb file name appears in uppercase.

// For Server.....
// I/P: pdb, int, hssp.....
   O/P: .ent, _convmod.int
   
    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
********************************************************************************************************** */
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<conio.h>

float entropycalc(char *);      /* Function to calculate sequence entropy */
char aacode(char *);	/* Generates the single-letter amino acid code */
       
int main(int argc, char *argv[]) // Usage here 1A3Q 0.90
{
      FILE *fpdb, *fint, *fhssp, *fent, *fmodint, *ftab;
      int i,j,k,knewint,kpdb,kint,khssp,flag,flag_core,core,rim;
      float mon_asa,dim_asa,sasa,coreent,riment,coresasa,rimsasa,seqent,core_ent[2500],rim_ent[2500];
      char *alignment,*aminoacidname;
      char pdb[100],intfile[100],hssp[300],ent[100],modint[100],linepdb[100],check[5],chain_pdb[10],chain_pdb_check,clus[20];
      char chain_int_check,chain_int[10],lineint[100],chain_hssp_check,chain_hssp[10],linehssp[200],check_hssp[15];
      char newchain_int[10],linemodint[100],intres_old[6],temp[5],resnum_int[5],intchain_old;
      char amino_int[4],monomer_asa[6],dimer_asa[6],intamino_old[4],chid_int,num_hssp[9],aacid,acid;
      char align[100000],core_resnum[2500][6],core_chain[2500],core_amino[2500][4],rim_resnum[2500][6];
      char rim_chain[2500],rim_amino[2500][4];

	char AlignmentStart[20],TempHssp_1[10],TempHssp_2[10],align_new[100000],align_without_null[100000];
	int  s,k1,NumHomologs,HomologID[100000];

      strcpy(pdb,"./");     strcat(pdb,argv[1]);     strcat(pdb,".pdb");
      strcpy(intfile,"./"); strcat(intfile,argv[1]); strcat(intfile,".int");     
      strcpy(hssp,"./");    strcat(hssp,argv[1]);    strcat(hssp,".hssp");
      
      strcpy(ent,"./");    strcat(ent,argv[1]);      strcat(ent,".ent");
      strcpy(modint,"./"); strcat(modint,argv[1]);   strcat(modint,"_convmod.modint");

      coreent=riment=coresasa=rimsasa=0.0f;core=rim=0;
      
//----------------- Reading pdb file-------------------
	// -------------- For reading chain IDs -------

      fpdb=fopen(pdb,"r");

      kpdb=0;chain_pdb_check='?';
      while(fgets(linepdb,100,fpdb)!=NULL)
      {
       for(i=0;i<4;i++)
              check[i]=linepdb[i];
       check[4]='\0';
       
       if(strcmp(check,"ATOM")==0)
       {
        if(chain_pdb_check!=linepdb[21])
        {
         chain_pdb_check=linepdb[21];
         chain_pdb[kpdb++]=linepdb[21];
        }
       }
      }
      fclose(fpdb);

	chain_pdb[kpdb]='\0'; // contains chain ids of pdb file
// -------------- pdb reading ends -------------------------------
     
// -------------------- writing in ent file ---------------------- 
// --------------------- chain IDs of pdb ------------------------
      fent=fopen(ent,"w");
      fprintf(fent,"Chain id(s) in pdbfile %s.pdb -",argv[1]);
      for(i=0;i<kpdb;i++)
            fprintf(fent,"%c ",chain_pdb[i]);
      fprintf(fent,"\n");
      fclose(fent);

// ------------------ writing ends -------------------------------

// ---------------------- reading int file -----------------------
      // --------------- for reading chain Ids --------------

      fint=fopen(intfile,"r");
      kint=0;chain_int_check='?';
      while(fgets(lineint,100,fint)!=NULL)
      {
       for(i=0;i<4;i++)
             check[i]=lineint[i];
       check[4]='\0';
       
       if(strcmp(check,"ATOM")==0)
       {
        if(chain_int_check!=lineint[21])
        {
         chain_int_check=lineint[21];
         chain_int[kint++]=lineint[21];
        }
       }
      }
      fclose(fint);

// ------------------- reading int ends --------------------------
     
// -------------------- writing in ent file in append mode ---------------------- 
// --------------------- chain IDs of int ------------------------ 
      fent=fopen(ent,"a");
      fprintf(fent,"Chain id(s) in intfile %s.int -",argv[1]);
      for(i=0;i<kint;i++)
            fprintf(fent,"%c ",chain_int[i]);
      fprintf(fent,"\n");
      fclose(fent);

// ------------------- writing ends ------------------------------
     
// -------------------- reading hssp file ---------------------- 
      // --------------- for reading chain Ids --------------
 
      fhssp=fopen(hssp,"r");
      
    if(fhssp == NULL)
    {
            fmodint=fopen(modint,"w"); // just creating the .modint file of size 0
            fclose(fmodint);
               
            exit(0);   
    }
      
      khssp=0;flag=0;chain_hssp_check='?';
      while(fgets(linehssp,200,fhssp)!=NULL && flag==0)
      {
       for(i=0;i<13;i++)
              check_hssp[i]=linehssp[i];
       check_hssp[13]='\0';
       
         if(strcmp(check_hssp,"## ALIGNMENTS")==0)
         {
           flag=1;

           fgets(linehssp,200,fhssp); // dummy to go to the next line
        
		do
        	{
         	fgets(linehssp,200,fhssp);

         	  if(chain_hssp_check!=linehssp[12] && linehssp[14]!='!' && linehssp[1]!='#')
		  {
          		chain_hssp_check=linehssp[12];
          		chain_hssp[khssp++]=linehssp[12];
         	  }
         	
		  for(i=0;i<2;i++)
                  check_hssp[i]=linehssp[i];
         	  check_hssp[2]='\0';

        	}while(strcmp(check_hssp,"##")!=0);

          }
      } // while ends

        fclose(fhssp);
	chain_hssp[khssp]='\0';

// ------------------------reading hssp ends -------------------------- 

// -------------------- writing in ent file in append mode ---------------------- 
// --------------------- chain IDs of hssp ------------------------ 
 
      fent=fopen(ent,"a");
      fprintf(fent,"Chain id(s) in hsspfile %s.hssp -",argv[1]);
      for(i=0;i<khssp;i++)
            fprintf(fent,"%c ",chain_hssp[i]);
      fprintf(fent,"\n");
      fclose(fent);

// ------------------- writing ends ------------------------------

// -------------------------matching chain IDs of hssp and pdb ------------      
      knewint=0;flag=0;
      for(i=0;i<kpdb;i++)
      {
        for(j=0;j<khssp;j++)
        {
          if(chain_pdb[i]==chain_hssp[j])
          {
            flag=1;
            newchain_int[knewint++]=chain_pdb[i]; // contains chain IDs would be written in medint
          }
        }
      }

      if(flag==0 && chain_pdb[0]=='A' && chain_hssp[0]==' ')
      {
       knewint=1;flag=2;
       newchain_int[0]='A';
      }

// ------------------------- matching ends ---------------------------------
      
// --------------- Reading hssp to get the ID number of  >90% or >80% homologs --------------
 
      fhssp=fopen(hssp,"r");
      flag=0; NumHomologs=0;
      while(fgets(linehssp,200,fhssp)!=NULL)
      {
	for(i=0;i<4;i++) check_hssp[i]=linehssp[28+i]; check_hssp[4]='\0';
        if( strcmp(check_hssp,"%IDE")==0)
	{
	  flag=1;
	  fgets(linehssp,200,fhssp);
	}

	for(i=0;i<13;i++) AlignmentStart[i]=linehssp[i]; AlignmentStart[i]='\0';
	if(strcmp(AlignmentStart,"## ALIGNMENTS")==0) // alignment starts at this position 
	{
	break;
	}
	
         if( flag==1 )
         {
		for(i=0;i<5;i++) TempHssp_1[i]=linehssp[i]; TempHssp_1[i]=='\0';			
		for(i=0;i<4;i++) TempHssp_2[i]=linehssp[i+28]; TempHssp_2[i]=='\0';

		if( atof(TempHssp_2) >= atof(argv[2]) )
		{	
			HomologID[NumHomologs]= atoi(TempHssp_1);
			NumHomologs++;					
		}
         }


      } // while ends

        fclose(fhssp);

	 // for(i=0;i<NumHomologs;i++) printf("%d\n",HomologID[i]);

// ------------------------reading hssp ends -------------------------- 

//---- Write out the modified .int file containing only those chains which are present in the hssp file ----

      fint=fopen(intfile,"r");
      fmodint=fopen(modint,"w");

      while(fgets(lineint,100,fint)!=NULL)
      {
       for(i=0;i<4;i++)
              temp[i]=lineint[i];
       temp[4]='\0';
       
       if(strcmp(temp,"ATOM")==0)
       {
              if(flag==2 && lineint[21]=='A')
              {
               //Copy .int to .modint but replace chain id 'A' with ' ' to match hssp file
               lineint[21]=' ';
               fputs(lineint,fmodint);
              }
              else
              {
               for(i=0;i<knewint;i++)
               {
                      if(lineint[21]==newchain_int[i])
                      fputs(lineint,fmodint);
               }
              }
       }

      } // while ends
      fclose(fint);
      fclose(fmodint);

// ------------------------------------- writing in .modint ends -----------------------------

// ------------------------ Calculate conservation scores using modint and hssp files --------
   
      strcpy(intres_old,"?????");intchain_old='?';strcpy(intamino_old,"???");

      fmodint=fopen(modint,"r");
      while(fgets(linemodint,100,fmodint)!=NULL) // while of modint
      {
       
	for(i=0;i<4;i++) // res no. 
	resnum_int[i]=(linemodint[22+i]==' ')?'0':linemodint[22+i];
	if(linemodint[26]==' ')resnum_int[4]='\0';
        else
        {
           resnum_int[4]=linemodint[26];
           resnum_int[5]='\0';
        }

        chid_int=linemodint[21]; // ch id
       
        for(i=0;i<3;i++) // res name
        amino_int[i]=linemodint[17+i];
        amino_int[3]='\0';
       
       while(strcmp(intres_old,resnum_int)!=0 || intchain_old!=chid_int)
       {
        // Begin to work on a new interface residue
        // Get the sequence alignment details from the hssp file and calculate the conservation score
        	k=0;
        
		// ---------------- hssp reading --------------------------------------------------
		fhssp=fopen(hssp,"r");
		while(fgets(linehssp,150,fhssp)!=NULL) // for hssp reading...........
        	{
         		for(i=0;i<4;i++)
			num_hssp[i]=(linehssp[7+i]==' ')?'0':linehssp[7+i];
			if(linehssp[11]==' ')num_hssp[4]='\0';
		        else
        		{
				num_hssp[4]=linehssp[11];
				num_hssp[5]='\0';
        		}

         		aacid=linehssp[14];
         
         		aminoacidname=amino_int;
         		acid=aacode(aminoacidname); // function call to convert res nm of .int files to 1 letter code
         
         		   if(strcmp(resnum_int,num_hssp)==0 && acid==aacid && chid_int==linehssp[12])
         		   {
          			j=0;
          
				while(linehssp[51+j]!='\0')
          			{
           				if(linehssp[49]=='\n')
					break;

           				align[k++]=linehssp[51+j];
					++j;
          			}
         		   }
        	}
        	align[k]='\0';//What if there are no aligned sequences at this position?
       
// ----------------- changed for homolog %IDE --------------------- ARP,SD 09.06.09 ----------
		//alignment=align;
        	//seqent=entropycalc(alignment); // function call to calculate entropy at that int res position
      
		//printf("align=%s\n",align);
 
		// ------- omitting '\n's from align[]
		k1=0;
		for(i=0;i<k;i++)
		{
			if( align[i]!='\n' )
			{
			align_without_null[k1]=align[i];
			k1++;
			}
		}
		align_without_null[k1]='\0';
		// ------------------------------------

		//printf("k=%d,k1=%d\n",k,k1);	
		s=0;
		for(i=0;i<k1;i++)
		{
		   //if(align_without_null[i]=='\n')printf("null\n");

		   for(j=0;j<NumHomologs;j++)
		   {
			if( HomologID[j] == (i+1) )
			{
				align_new[s]=align_without_null[i];
				s++;
			}
		   }
		}
		align_new[s]='\0';	
		
		alignment=align_new;
		seqent=entropycalc(alignment);

        	fclose(fhssp);
// ---------------------------------------------------------------------------------------------
	// ------------------------------ hssp reading ends ---------------------------------
               
        
		flag_core=0;sasa=0.0f;
        	do
        	{        	
		strcpy(intres_old,resnum_int);intchain_old=chid_int;strcpy(intamino_old,amino_int);
         
         	for(i=0;i<5;i++)
                  monomer_asa[i]=(linemodint[67+i]==' ')?'0':linemodint[67+i];   // was 55
         	  monomer_asa[5]='\0';
         	mon_asa=atof(monomer_asa);
         
         	for(i=0;i<5;i++)
                  dimer_asa[i]=(linemodint[73+i]==' ')?'0':linemodint[73+i]; // was 61
         	  dimer_asa[5]='\0';
         	dim_asa=atof(dimer_asa);
         
         	//if(strcmp(dimer_asa,"00.00")==0)
         	if(dim_asa == 0.0)
                  flag_core=1;
                  
         	sasa+=(mon_asa-dim_asa);
         
         	if(fgets(linemodint,100,fmodint)==NULL)
                break;
                  
         	  for(i=0;i<4;i++)
                  resnum_int[i]=(linemodint[22+i]==' ')?'0':linemodint[22+i];
         
			if(linemodint[26]==' ')resnum_int[4]='\0';
         		else
         		{
             			resnum_int[4]=linemodint[26];
             			resnum_int[5]='\0';
         		}
         
         	  chid_int=linemodint[21];
         
         	   for(i=0;i<3;i++)
              	   amino_int[i]=linemodint[17+i];
         	   amino_int[3]='\0';
         
        	}while(strcmp(intres_old,resnum_int)==0 && intchain_old==chid_int && strcmp(intamino_old,amino_int)==0);
        
        //Save details of interface residue in array to print later to ent file
        //One set of arrays for the core residues, another for the rim residues
        
        if(flag_core==1)
        {
         strcpy(core_resnum[core],intres_old);core_chain[core]=intchain_old;strcpy(core_amino[core],intamino_old);

	if(k1>=1)//k<1 signifies that there are no aligned sequences for that residue position. Hence it is assigned a neagtive value and also not used in the average entropy calculation
		core_ent[core++]=seqent;
	else
		core_ent[core++]=-0.0; // If there are no aligned sequences for a particular interface residue

		if(k1>=1)
		{
                coreent+=(seqent*sasa);
                coresasa+=sasa;
		}
        }
        if(flag_core==0)
        {
         strcpy(rim_resnum[rim],intres_old);rim_chain[rim]=intchain_old;strcpy(rim_amino[rim],intamino_old);
	
		if(k1>=1)
		rim_ent[rim++]=seqent;
		else
		rim_ent[rim++]=-0.0;
	
		if(k1>=1)
		{
                riment+=(seqent*sasa);
                rimsasa+=sasa;
		}
        }

       } // while for each residue of modint ends

      } // while for modint ends
      fclose(fmodint);
      
      fent=fopen(ent,"a");
      fprintf(fent,"\n************************************\nNumber of homologs used=%6d\n",NumHomologs);
      fprintf(fent,"Mean entropy of core residues=%.4f\n",coreent/coresasa);
      fprintf(fent,"Mean entropy of rim residues=%.4f\n",riment/rimsasa);
      fprintf(fent,"Entropy ratio (Core/Rim): %.2f\n",(coreent/coresasa)/(riment/rimsasa));
      fprintf(fent,"************************************\n\n");
      fprintf(fent,"Res.# Chain Res.name Seq.entropy\n");
      fprintf(fent,"CORE residues:\n");
      for(i=0;i<core;i++)fprintf(fent,"%s %c %s %.4f\n",core_resnum[i],core_chain[i],core_amino[i],core_ent[i]);
      fprintf(fent,"\nRIM residues:\n");
      for(i=0;i<rim;i++)fprintf(fent,"%s %c %s %.4f\n",rim_resnum[i],rim_chain[i],rim_amino[i],rim_ent[i]);
      fclose(fent);

      //getch();

      return 0;

} // main ends



float entropycalc(char *al)
{
	int b=0;
	float class1,class2,class3,class4,class5,class6,class7;
	float entropy=0.0f;
	char c;

	class1=class2=class3=class4=class5=class6=class7=0.0f;

		//printf("last::%s\n",al); 
	while(*al != '\0')
	{
		c=*al;
		if(c=='A'||c=='V'||c=='L'||c=='I'||c=='M'||c=='C'||c=='a'||c=='v'||c=='l'||c=='i'||c=='m'||c=='c')
			++class1;
		else if(c=='G'||c=='S'||c=='T'||c=='s'||c=='t'||c=='g')
			++class2;
		else if(c=='D'||c=='E'||c=='d'||c=='e')
			++class3;
		else if(c=='N'||c=='Q'||c=='n'||c=='q')
			++class4;
		else if(c=='R'||c=='K'||c=='r'||c=='k')
			++class5;
		else if(c=='P'||c=='F'||c=='Y'||c=='W'||c=='p'||c=='f'||c=='y'||c=='w')
			++class6;
		else if(c=='H'||c=='h')
			++class7;
		else if(c==' '||c=='.'||c=='\n')
			--b;
/*		++b;
		if(c=='\n')
			--b;
		if(c==' '||c=='.')
			--b;*/
		++b;
		++al;
	}
	entropy=0.0f;
	if(class1 != 0)
		entropy=entropy+((class1/b)*log(class1/b));
	if(class2 != 0)
		entropy=entropy+((class2/b)*log(class2/b));
	if(class3 != 0)
		entropy=entropy+((class3/b)*log(class3/b));
	if(class4 != 0)
		entropy=entropy+((class4/b)*log(class4/b));
	if(class5 != 0)
		entropy=entropy+((class5/b)*log(class5/b));
	if(class6 != 0)
		entropy=entropy+((class6/b)*log(class6/b));
	if(class7 != 0)
		entropy=entropy+((class7/b)*log(class7/b));
	if(entropy != 0)
		entropy=-entropy;       /* adding a negative sign */
	return(entropy);
}

char aacode(char *aa)
{
	int i=0;
	char amino[4];

	while(*aa!='\0')
		amino[i++]=*aa++;
	amino[i]='\0';

	if(strcmp(amino,"ALA")==0)
		return('A');
	else if(strcmp(amino,"VAL")==0)
		return('V');
	else if(strcmp(amino,"LEU")==0)
		return('L');
	else if(strcmp(amino,"ILE")==0)
		return('I');
	else if(strcmp(amino,"GLY")==0)
		return('G');
	else if(strcmp(amino,"PRO")==0)
		return('P');
	else if(strcmp(amino,"TYR")==0)
		return('Y');
	else if(strcmp(amino,"TRP")==0)
		return('W');
	else if(strcmp(amino,"PHE")==0)
		return('F');
	else if(strcmp(amino,"ARG")==0)
		return('R');
	else if(strcmp(amino,"LYS")==0)
		return('K');
	else if(strcmp(amino,"ASP")==0)
		return('D');
	else if(strcmp(amino,"GLU")==0)
		return('E');
	else if(strcmp(amino,"ASN")==0)
		return('N');
	else if(strcmp(amino,"GLN")==0)
		return('Q');
	else if(strcmp(amino,"CYS")==0)
		return('C');
	else if(strcmp(amino,"HIS")==0)
		return('H');
	else if(strcmp(amino,"MET")==0)
		return('M');
	else if(strcmp(amino,"SER")==0)
		return('S');
	else if(strcmp(amino,"THR")==0)
		return('T');
}
