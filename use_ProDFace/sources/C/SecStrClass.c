/* **********************************************************************************************************
   THIS PROGRAM DIVIDES THE PROTEINS IN OUR DATASETS INTO SECONDARY STRUCTURE CLASSES DEPENDING UPON 
   THE SECONDARY STRUCTURE COMPOSITION OF THE INTERFACE (BOTH SUBUNITS TAKEN TOGETHER). 
   WE CONSIDER INTERFACES AS ALL-ALPHA, ALL-BETA, MIXED ALPHA-BETA AND NR INTERFACES 
   DEPENDING UPON THE FOLLOWING CRITERIA:-
	1) ALL-ALPHA INTERFACES - % INTERFACE RESIDUES IN ALPHA(>=60%), % IN BETA(<5%);OR,
				  % INTERFACE RESIDUES IN ALPHA(>40%), % IN BETA(=0%)
	2) ALL-BETA INTERFACES - % INTERFACE RESIDUES IN BETA(>=60%), % IN ALPHA(<5%);OR,
				 % INTERFACE RESIDUES IN BETA(>40%), % IN ALPHA(=0%)
	3) MIXED ALPHA-BETA INTERFACES - % INTERFACE RESIDUES IN ALPHA(>=20%), % IN BETA(>=20%)
	4) NR INTERFACES - % INTERFACE RESIDUES IN NR REGIONS(>=60%)

   FURTHER, THE SECONDARY STRUCTURE COMPOSITION AT THE PROTEIN-PROTEIN INTERFACES IS ALSO COMPARED TO THE 
   OVERALL SECONDARY STRUCTURE COMPOSITION OF THE PROTEIN. THIS WILL GIVE US AN IDEA OF HOW CORRELATED THESE 
   TWO QUANTITIES ARE- SEC.STR. IN THE OVERALL PROTEIN VS. SEC.STR. IN THE INTERFACE FORMED BY THE SAME PROTEIN.
   

	I/P :: 1A3Q.int, 1A3Q.dssp
    O/P :: 1A3Q.SecStrClass

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
   ********************************************************************************************************** */

/* This program is now being used for the protein-dna complexes. Date-12-04-07 */
#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdlib.h>
#include<stddef.h>

main(int argc, char *argv[])
{
	FILE *flist, *fmodint, *fdssp, *fout;
	int i, k, prot;
	float alpha_num,beta_num,mixed_num,nr_num,alpha_int,beta_int,nr_int,int_res,tot_res,alpha_tot,beta_tot;
    float nr_tot,alpha_frac_int,beta_frac_int,nr_frac_int,alpha_frac_tot,beta_frac_tot,nr_frac_tot;
	char pdb[100], out[100], pdbid[5],inputfile_int[100],oldnum[6],lineint[90],num[6],str_int,inputfile_dssp[100];
    char linedssp[150],str_dssp, chain_given[3];

	strcpy(out,argv[1]); strcat(out,".SecStrClass");
    fout=fopen(out,"w");

	prot=0;
	alpha_num=beta_num=mixed_num=nr_num=0.0f;
			
		strcpy(pdb,argv[1]);

		strcpy(inputfile_dssp,pdb);strcat(inputfile_dssp,".dssp");

// ---------------------------------------------------------------------------------
     alpha_int=beta_int=nr_int=0.0f;int_res=0.0f;
     for(k=0; k<(argc-2); k++)
     {
	//	strcpy(inputfile_int,pdb); strcat(inputfile_int,"mod_sssB.int");
        strcpy(inputfile_int,pdb); strcat(inputfile_int,"mod_sss"); strcat(inputfile_int,argv[k+2]); // Here is argv[2]
        strcat(inputfile_int,".int");
        
		fmodint=fopen(inputfile_int,"r");

		strcpy(oldnum,"?????");
		
		while(fgets(lineint,90,fmodint)!=NULL)
		{
			for(i=0;i<4;i++)
			{
				if(lineint[22+i]==' ')
					num[i]='0';
				else
					num[i]=lineint[22+i];
			}
			if(lineint[26]!=' ')
			{
				num[4]=lineint[26];num[5]='\0';
			}
			else
				num[4]='\0';

			if(strcmp(oldnum,num)!=0)
			{
				++int_res;

				str_int=lineint[79];

				if(str_int=='H')      alpha_int++;
                else if(str_int=='S') beta_int++;
                else if(str_int=='T') nr_int++;

				strcpy(oldnum,num);
			}
		}
		fclose(fmodint);

      } // end of for 'k'	
// --------------------------------------------------------------------------------		

		alpha_frac_int=(alpha_int/int_res)*100;
        beta_frac_int=(beta_int/int_res)*100;
        nr_frac_int=(nr_int/int_res)*100;

/*		if((alpha_frac_int>=60.0f && beta_frac_int<5.0f) || (alpha_frac_int>40.0f && beta_frac_int==0.0f))
			++#alpha;
		else if((beta_frac_int>=60.0f && alpha_frac_int<5.0f) || (beta_frac_int>40.0f && alpha_frac_int==0.0f))
			++#beta;
		else if(alpha_frac_int>=20.0f && beta_frac_int>=20.0f)
			++#mixed;
		else if(nr_frac_int>=60.0f)
			++#nr;*/

		fdssp=fopen(inputfile_dssp,"r");

		alpha_tot=beta_tot=nr_tot=0.0f;tot_res=0.0f;

		while(fgets(linedssp,150,fdssp)!=NULL)
		{
			if(linedssp[2]=='#')
			{
				while(fgets(linedssp,150,fdssp)!=NULL)
				{
					if(linedssp[13]!='!')
					{
						++tot_res;

						str_dssp=linedssp[16];

						if(str_dssp=='H'||str_dssp=='G')++alpha_tot;else if(str_dssp=='E'||str_dssp=='B')++beta_tot;else if(str_dssp=='T'||str_dssp=='S'||str_dssp==' ')++nr_tot;
					}
				}
			}
		}
		fclose(fdssp);

		alpha_frac_tot=(alpha_tot/tot_res)*100;beta_frac_tot=(beta_tot/tot_res)*100;nr_frac_tot=(nr_tot/tot_res)*100;

		if(alpha_frac_int>=40.0f && beta_frac_int<10.0f)
/*		if((alpha_frac_int>=60.0f && beta_frac_int<5.0f) || (alpha_frac_int>40.0f && beta_frac_int==0.0f))*/
		{
			++alpha_num;++prot;
			fprintf(fout,"ALPHA  %6.2f  %6.2f  %6.2f\n",alpha_frac_int, beta_frac_int, nr_frac_int);
		}
		else if(beta_frac_int>=40.0f && alpha_frac_int<10.0f)
/*		else if((beta_frac_int>=60.0f && alpha_frac_int<5.0f) || (beta_frac_int>40.0f && alpha_frac_int==0.0f))*/
		{
			++beta_num;++prot;
			fprintf(fout,"BETA   %6.2f  %6.2f  %6.2f\n",alpha_frac_int,beta_frac_int,nr_frac_int);
		}
/*		else if(alpha_frac_int>=20.0f && beta_frac_int>=20.0f)
		{
			++mixed_num;
			fprintf(fout,"%s\tMIXED\t%f\t%f\t%f\t%f\t%f\n",pdb,alpha_frac_int,beta_frac_int,nr_frac_int,alpha_frac_tot,beta_frac_tot);
		}*/
		else if(nr_frac_int>=60.0f)
		{
			++nr_num;++prot;
			fprintf(fout,"NR     %6.2f  %6.2f  %6.2f\n",alpha_frac_int,beta_frac_int,nr_frac_int);
		}
		else
		{
			++mixed_num;++prot;
			fprintf(fout,"MIXED  %6.2f  %6.2f  %6.2f\n",alpha_frac_int,beta_frac_int,nr_frac_int);
		}

	fclose(fout);
	//getch();
}
