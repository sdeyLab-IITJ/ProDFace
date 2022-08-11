/*  **********************************************************************************************************
	1.	This program reads an .int file (interface file made by intf.f)
            and writes the res name,res no and chain id in a '.out' file,10 res per line.

	I/P::	1a4ia.int
	O/P::	1a4ia.out

eg. of 1a4ia.out
------------------------------------------------------------------------------
1a4ia int    38
GLU 112A ASN 116A SER 129A ILE 130A ALA 132A GLY 133A ARG 134A LEU 135A ALA 136A ARG 137A
ASP 139A ASP 142A PRO 146A ALA 164A GLY 165A ARG 166A HIS 167A ARG 173A LYS 175A ALA 179A
PRO 180A HIS 182A ASP 183A LEU 186A TRP 187A ASN 189A ALA 190A THR 191A VAL 192A THR 193A
THR 194A CYS 195A HIS 196A LYS 198A THR 199A ALA 200A GLU 205A LYS 208A
------------------------------------------------------------------------------

	2.	Then it reads this .out file and calculates the gap between consicutive
	   residues.If the gap is within the given gap length,it considers it in the same
	   fragment,otherwise in a different fragment.It writes the informations in
	   .segment file.

	O/P::	1a4ia.segment

eg. of 1a4ia.segment
------------------------------------------------------------------------------
GLU 112A ASN 116A
SER 129A ILE 130A ALA 132A GLY 133A ARG 134A LEU 135A ALA 136A ARG 137A ASP 139A ASP 142A PRO 146A
ALA 164A GLY 165A ARG 166A HIS 167A
ARG 173A LYS 175A ALA 179A PRO 180A HIS 182A ASP 183A LEU 186A TRP 187A ASN 189A ALA 190A THR 191A
VAL 192A THR 193A THR 194A CYS 195A HIS 196A LYS 198A THR 199A ALA 200A GLU 205A LYS 208A
4:::2 11 4 21
-----------------------------------------------------------------------------

	3.	Also it writes all information in 'segment_master'.

	O/P::	segment_master

eg. of segment_master
-----------------------------------------------------------------------------
12asa.int  52 1889.34   8:::   7 10 10  4  9  2  9  1
12asb.int  27  852.98   3:::   2  2 23
-----------------------------------------------------------------------------
Here,the 3rd column is for Buried asa i.e "asa[A]=asa[AB]".........

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
************************************************************************************************************ */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>

#define MAX 200
#define MAX1 20

main(int argc, char *argv[])
{
    FILE *fp_list_int,*fp_list_out;
    FILE *fp_name_int,*fp_name_out;
    FILE *fp_name_out_read;
    FILE *fp_name_segment;
    FILE *fp_segment_master;
    
    char line_list_int[MAX1],line_int[MAX],line_out[MAX];
    char test1[10],RES_NAME[10],RES_NO[9],CHAIN_ID,INT_NAME[9],ipfile1[100],ipfile2[100],ipfile3[100]; 
    int i,c,c1,res_no1,test2;

    char asa_1[8],asa_2[8];
    float d_asa;

    char test[5],res_nm[4],res_no[5],ch_id[1000],a_res_nm[1000][4];
    int c2,c3,k,z,gap,reqd_gap,length,a_res_no[1000],chain_len[1000];

	reqd_gap=4;
	
	strcpy(ipfile1,argv[1]); strcat(ipfile1,argv[2]); strcat(ipfile1,".int");
    strcpy(ipfile2,argv[1]); strcat(ipfile2,argv[2]); strcat(ipfile2,"pps.out");
    strcpy(ipfile3,argv[1]); strcat(ipfile3,argv[2]); strcat(ipfile3,".segment");

/* ***************************************************************************
        		WORKING WITH .int AND .out FILE
*************************************************************************** */ 
  
			fp_name_int=fopen(ipfile1,"r");
			fp_name_out=fopen(ipfile2,"w");			

	/* --- while loop to calculate no. of res in .int file --- */
      strcpy(test1,"XXX");
      test2=9999;
      c1=0;d_asa=0.0;
       while(fgets(line_int,MAX,fp_name_int)!=NULL)                
       {   
          for(i=0;i<3;i++)RES_NAME[i]=line_int[17+i];
          RES_NAME[3]='\0';

          for(i=0;i<4;i++)RES_NO[i]=line_int[22+i];
          RES_NO[i]='\0';

          res_no1=atoi(RES_NO);
          if(strcmp(test1,RES_NAME)!=0 || test2!=res_no1)
            {
              c1++;
              strcpy(test1,RES_NAME);
              test2=res_no1;
            }

          for(i=0;i<6;i++)asa_1[i]=line_int[54+i];             /* ****************** */
          asa_1[6]='\0';
          for(i=0;i<6;i++)asa_2[i]=line_int[60+i];             /* ASA calculation    */
          asa_2[6]='\0';
          d_asa=d_asa+( (atof(asa_1)) - (atof(asa_2)) );   /* ****************** */


       }                                             
          fprintf(fp_name_out,"%s    %d\n",ipfile1,c1); /* 1st line of .out file */ 

      			fclose(fp_name_int);
	/* -------------------------------------------------------- */    

	/* --- while loop to read res name,no.,and chain id from .int file
	       and print them in .out file,  10 res per line............. --- */
      
			fp_name_int=fopen(ipfile1,"r");

	  strcpy(test1,"XXX");
          test2=9999;
          c=0;
       while(fgets(line_int,MAX,fp_name_int)!=NULL)            
       {
          for(i=0;i<3;i++)RES_NAME[i]=line_int[17+i];
          RES_NAME[3]='\0';

          for(i=0;i<4;i++)RES_NO[i]=line_int[22+i];
          RES_NO[i]='\0';

          res_no1=atoi(RES_NO);
        
          CHAIN_ID=line_int[21];

          if(strcmp(test1,RES_NAME)!=0 || test2!=res_no1)
            {
               fprintf(fp_name_out,"%s%4d%c ",RES_NAME,res_no1,CHAIN_ID);
               c++; 
            }
          if(c==10){fprintf(fp_name_out,"\n");c=0;}
          
          strcpy(test1,RES_NAME);
          test2=res_no1;
        }
			fclose(fp_name_int);  
			fclose(fp_name_out);  


/*********************************************************************
	WORKING WITH .out AND.segment FILE
******************************************************************* */

	fp_name_segment=fopen(ipfile3,"w");
	fp_name_out_read=fopen(ipfile2,"r");

	/* ---- To open the .out file nd read res names & no.s in arrays --- */

c2=0;
while(fgets(line_out,MAX,fp_name_out_read)!=NULL)
{
    for(i=0;i<4;i++)test[i]=line_out[25+i];
    test[4]='\0';

    if(strcmp(test,".int")!=0)
     {
      length=strlen(line_out);

       for(k=0;k<(length-1);)
       {
        for(i=0;i<3;i++)res_nm[i]=line_out[i+k];
        res_nm[3]='\0';
        strcpy(a_res_nm[c2],res_nm);

        for(i=0;i<4;i++)res_no[i]=line_out[3+i+k];
        res_no[4]='\0';
        a_res_no[c2]=atoi(res_no);

        ch_id[c2]=line_out[7+k];

        k=k+9;
        c2++;
       }

     }
}
	/* --- To calculate the 'gap' between two successive residue
	       and if the gap is within given 'gap length',print the
	       res name,no. and chain id in .segment file.

	Also to print total no. of fragment and no. of each fragment
	       at the bottom of .segment file. ------------------- */

	c3=0;z=0;
    for(i=0;i<c2;i++)
    {
        if(i<c2-1)
        {
         gap=(a_res_no[i+1]-a_res_no[i]);
         if(gap<=reqd_gap)
          {
            fprintf(fp_name_segment,"%s%4d%c ",a_res_nm[i],a_res_no[i],ch_id[i]);
            z++;
           }

         if(gap>reqd_gap)
          {
            fprintf(fp_name_segment,"%s%4d%c ",a_res_nm[i],a_res_no[i],ch_id[i]);
            fprintf(fp_name_segment,"\n");
           chain_len[c3]=z+1;
           c3++;
           z=0;
          }
        }
         if(i==(c2-1))
         {
          gap=(a_res_no[i]-a_res_no[i-1]);
          if(gap<=reqd_gap)
          {fprintf(fp_name_segment,"%s%4d%c ",a_res_nm[i],a_res_no[i],ch_id[i]);
           chain_len[c3]=z+1;
           }
          if(gap>reqd_gap)
          {
            fprintf(fp_name_segment,"%s%4d%c ",a_res_nm[i],a_res_no[i],ch_id[i]);
           chain_len[c3]=z+1;
          }
         }

       }
	fprintf(fp_name_segment,"\n%3d:::",(c3+1));

	for(i=0;i<=c3;i++)
	fprintf(fp_name_segment,"%d ",chain_len[i]);
	
	fclose(fp_name_out_read);
	fclose(fp_name_segment);
 
}
