/* **********************************************************************************************************
// Gives average entropy <S> of a subunit

	I/P :: 1A3QA.pdb, 1A3Q.hssp
    O/P :: 1A3Q_c1_avg_subunit_entropy

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
********************************************************************************************************** */

#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>

float entropycalc(char *);      /* Function to calculate sequence entropy */
char aacode(char *);	/* Generates the single-letter amino acid code */
       
int main(int argc, char *argv[]) // Usage here 12as .7
{
      FILE *fpdb, *fhssp, *ftab, *fent;
      int i,j,k,knewint,kpdb,kint,khssp,flag,flag_core,core,rim;
      double mon_asa,dim_asa,sasa,coreent,riment,coresasa,rimsasa,seqent,core_ent[2500],rim_ent[2500];
      char *alignment,*aminoacidname;
      char pdb[90],hssp[300],out[90],linepdb[100],check[10],chain_pdb[10],chain_pdb_check,clus[20];
      char chain_hssp_check,chain_hssp[10],linehssp[200],check_hssp[15],newchain_int[10];
      char num_hssp[5],aacid,acid;
      char align[100000],core_resnum[2500][6],core_chain[2500],core_amino[2500][4],rim_resnum[2500][6];
      char rim_chain[2500],rim_amino[2500][4];

	  char AlignmentStart[20],AlignmentStop[35],TempHssp_1[10],TempHssp_2[10],align_new[100000],align_without_null[100000];
	  int  n,n1,s,k1,NumHomologs,HomologID[100000],TotResPdb,TotResSu;
	  char PdbResNo[100000][6],PdbChId[100000],TResNm[4],check_1[6],TResNo[6],hsspResNo[9];
	  double TotEnt;

      strcpy(pdb,"./");     strcat(pdb,argv[1]);   strcat(pdb,".pdb");
      strcpy(hssp,"./");    strcat(hssp,argv[1]);  strcat(hssp,".hssp");
      
      strcpy(out,"./");     strcat(out,argv[1]);   strcat(out,"_c1_avg_subunit_entropy");
      
		ftab=fopen(out,"w");
		fprintf(ftab,"%s ",argv[1]);

      
//----------------- Reading pdb file-------------------
	// -------------- For reading chain IDs -------

      fpdb=fopen(pdb,"r");

      kpdb=0;chain_pdb_check='?';
      while(fgets(linepdb,100,fpdb)!=NULL)
      {
          for(i=0;i<4;i++) check[i]=linepdb[i];
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

      fpdb=fopen(pdb,"r");

	// -----------------------------------------------

	// Reading again for res no. and ch ids

      strcpy(TResNm,"XXX"); strcpy(TResNo,"99999");
      TotResPdb=0;
      while(fgets(linepdb,100,fpdb)!=NULL)
      {
          for(i=0;i<4;i++) check[i]=linepdb[i];  
          check[4]='\0';

          if(strcmp(check,"ATOM")==0)
          {
              for(i=0;i<3;i++) check[i]=linepdb[17+i]; check[i]='\0';
	          for(i=0;i<5;i++) check_1[i]=linepdb[22+i]; check[i]='\0';

	          if( (strcmp(check,TResNm)!=0) || (strcmp(check_1,TResNo)!=0) )	
	          {
		          strcpy(PdbResNo[TotResPdb],check_1);
		          PdbChId[TotResPdb] = linepdb[21];
		          TotResPdb++;

		          strcpy(TResNm,check);
		          strcpy(TResNo,check_1);
	           }	

           }
      }
      fclose(fpdb);

	// -----------------------------------------------

// -------------- pdb reading ends -------------------------------

// -------------------- reading hssp file ---------------------- 
      // --------------- for reading chain Ids --------------
 
      fhssp=fopen(hssp,"r");
      
      if(fhssp == NULL)
      {               
            exit(0);   
      }
      
      khssp=0;flag=0;chain_hssp_check='?';
      while(fgets(linehssp,200,fhssp)!=NULL && flag==0)
      {
           for(i=0;i<13;i++) check_hssp[i]=linehssp[i];
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

               }  while(strcmp(check_hssp,"##")!=0);
           }
      } // while ends

        fclose(fhssp);
	    chain_hssp[khssp]='\0';

// ------------------------reading hssp ends -------------------------- 

// -------------------------matching chain IDs of hssp and pdb ------------      
      knewint=0;flag=0;
      for(i=0;i<kpdb;i++)
      {
        for(j=0;j<khssp;j++)
        {
          if(chain_pdb[i]==chain_hssp[j])
          {
            flag=1;
            newchain_int[knewint++]=chain_pdb[i]; // contains chain IDs existing in both pdb & hssp
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

// --------- Calculating conservation scores using pdb res no.+ ch id and hssp files --------

        TotEnt=0.0;   TotResSu=0;
        for(n1=0; n1<knewint; n1++) //newchain_int[knewint++
        {

           //TotEnt=0.0;   TotResSu=0;
           for(n=0; n<TotResPdb; n++)
           {
              if(newchain_int[n1]==PdbChId[n])
              {

	       // --------- reading hssp to collect the alignment for a particular res   ---------

                  fhssp=fopen(hssp,"r");
	              flag=0; k=0;
	              while(fgets(linehssp,200,fhssp)!=NULL)
                  {
		             for(i=0;i<5;i++) hsspResNo[i] = linehssp[7+i];
		             hsspResNo[i]='\0';

        	         for(i=0;i<13;i++) AlignmentStart[i]=linehssp[i]; AlignmentStart[i]='\0';
        	         
                     if(strcmp(AlignmentStart,"## ALIGNMENTS")==0) // alignment starts at this position 
        	         {
        	            flag=1;
		                fgets(linehssp,200,fhssp);
   	                 }
        
        	         for(i=0;i<30;i++) AlignmentStop[i]=linehssp[i]; AlignmentStop[i]='\0';
        	         
                     if(strcmp(AlignmentStart,"## SEQUENCE PROFILE AND ENTROPY")==0) // alignment stops at this position 
        	         {
        	            flag=0;
      	             }

         	         if( flag==1 )
         	         {
			            if(strcmp(PdbResNo[n],hsspResNo)==0 && PdbChId[n]==linehssp[12])
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
  	                  } // if ends

                  } // while ends
                  
                  fclose(fhssp);
        	
                  align[k]='\0';//What if there are no aligned sequences at this position?

	// -------------------- hssp reading ends ----------------------

	// ----------------- homolog %IDE ------------------

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
                
		        TotEnt += seqent;

                // ---------- writing in ent -----------
		           //fent=fopen("ent","a");
		           //fprintf(fent,"%s %c %.4f\n",PdbResNo[n],PdbChId[n],seqent);
		           //fclose(fent);
                // -------------------------------------

		        TotResSu++;

              } // if ends

           } // inner for ends

		// ---------- writing in tabulation file --------------

		//fprintf(ftab,"%c %.4lf ",newchain_int[n1],(TotEnt/(TotResSu*1.0)));

		// ----------------------------------------------------
		
        } // outer for ends
		
        fprintf(ftab,"%4.2lf ", (TotEnt/(TotResSu*1.0)));
        //fprintf(ftab,"\n");
		fclose(ftab);	

// ---------------------------------------------------------------------------------------------
        
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
