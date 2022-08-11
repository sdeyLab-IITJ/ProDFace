/* **********************************************************************************************************
 THIS PROGRAM GENERATES THE LIST OF INTERFACE SECONDARY STRUCTURE SEGMENTS THAT ARE PRESENT IN THE DATASETS.
   THE GENERATED FILE CONTAINS THE SSS NAME, INTERFACE SEGMENT AND THE ENTIRE SEGMENT IN THE PROTEIN
   
  	I/P :: 1A3Qmod_sssA.int (similar to 1A3Q.int)
    O/P :: 1A3Q.SSSegment

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com 

 This program is now being made ready for working on the protein-dna complexes.
********************************************************************************************************** */

#include<stdio.h>
#include<string.h>
#include<ctype.h>
#include<stdlib.h>
#include<stddef.h>

main(int agrc, char *argv[])
{
	FILE *flist, *fmodint, *fdssp, *fout;
	int i,j,flag,helix,strand,nr,start,end,length;

	char linelist[10],lineint[90],pdb[100],inputfile_int[50],oldstr[5],str[5],resnum[5];
    char end_res[5],start_res[5],chain, outfile[100];

	strcpy(outfile,argv[1]);
    strcat(outfile,".SSSegment"); //strcat(outfile,argv[2]); 
    fout=fopen(outfile,"a");
    
        strcpy(pdb,argv[1]);
        fprintf(fout,"%s%s\n",argv[1],argv[2]);

		strcpy(inputfile_int,pdb); strcat(inputfile_int,"mod_sss"); strcat(inputfile_int,argv[2]); 
        strcat(inputfile_int,".int");

		fmodint=fopen(inputfile_int,"r");

		strcpy(oldstr,"????");flag=0;helix=strand=nr=0;

		while(fgets(lineint,90,fmodint)!=NULL)
		{
			j=0;
			while(lineint[79+j]!='\n')
			{
				str[j]=lineint[79+j];
				++j;
			}
			str[j]='\0';

			if(strcmp(oldstr,str)!=0 )
			{
				if(strcmp(oldstr,"????")!=0 )
				{
					flag=0;
					strcpy(end_res,resnum);end=atoi(end_res);//Trial modification
					length=(end-start)+1;
					fprintf(fout,"%5s  %3d\n",end_res,length);
				}

				strcpy(oldstr,str);

				for(i=0;i<4;i++)
				{
					if(lineint[22+i]==' ')
						resnum[i]='0';
					else
						resnum[i]=lineint[22+i];
				}
				resnum[4]='\0';

				chain=lineint[21];

				if(flag==0 && lineint[79]=='H')
				{
					strcpy(start_res,resnum);fprintf(fout,"Helix  #%3d  %c  %5s  ",++helix,chain,start_res);

					start=atoi(start_res);

					flag=1;
				}
				else if(flag==0 && lineint[79]=='S')
				{
					strcpy(start_res,resnum);fprintf(fout,"Strand #%3d  %c  %5s  ",++strand,chain,start_res);
					start=atoi(start_res);//Trial modification
					flag=1;
				}
				else if(flag==0 && lineint[79]=='T')
				{
					strcpy(start_res,resnum);fprintf(fout,"NR     #%3d  %c  %5s  ",++nr,chain,start_res);
					start=atoi(start_res);//Trial modification
					flag=1;
				}
			}

			for(i=0;i<4;i++)
				resnum[i]=(lineint[22+i]==' ')?'0':lineint[22+i];
			resnum[4]='\0';
		}
		fclose(fmodint);

		flag=0;
		strcpy(end_res,resnum);end=atoi(end_res);//Trial modification
		length=(end-start)+1;
		fprintf(fout,"%5s  %3d\n",end_res,length);
	
	fclose(fout);
}

