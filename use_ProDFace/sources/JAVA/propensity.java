/* ********************************************************************
   PROG TO FIND OUT THE DENSITY OF INTERFACE ATOMS
   ALSO CALCULATES "Residue Propensity Score"

    I/P :: 1A3Q_protein.int
    O/P :: stat1A3Q.dat

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
*********************************************************************** */

import java.text.*;
import java.io.*;

class propensity
{
   public static void main(String args[])
   {
   RandomAccessFile in, out;
   String path, pdb, pdb_fil1, pdb_fil2;
   String s, s1 = "", test1, res_nm;
   int test2, res_no;
   float prop_ALA, prop_ARG, prop_ASN, prop_ASP, prop_CYS;
   float prop_GLN, prop_GLU, prop_GLY, prop_HIS, prop_ILE;
   float prop_LEU, prop_LYS, prop_MET, prop_PHE, prop_PRO;
   float prop_SER, prop_THR, prop_TRP, prop_TYR, prop_VAL;
   float prop_TOT;
   DecimalFormat rnd = new DecimalFormat("0.00");

   if (args.length != 1)
   {
      System.out.println("Please pass the .int filename as argument. Halting Program.");
      return;
   }

   path = args[0].substring(0,7);
   pdb = args[0].substring(7,11);
   pdb_fil1 = args[0].substring(7,13) + ".int";
   pdb_fil2 = "stat" + pdb + ".dat";

   test1 = "XXX";
   test2 = 9999;
   prop_TOT = 0;

   try
   {

   in = new RandomAccessFile (path + pdb_fil1,"r");
   out = new RandomAccessFile (path + pdb_fil2,"rw");
   //out.setLength(0);
   out.seek(out.length());

   while ((s = in.readLine()) != null)   // start of while
   {
   if (s.startsWith("ATOM"))   // start of if
   {
      res_nm = s.substring(17,20);
      res_no = Integer.parseInt(s.substring(22,26).trim());
      if ((res_nm.equals(test1)==false) || (res_no != test2))
      {
         test1 = res_nm;
         test2 = res_no;
         if (res_nm.equals("ALA")==true)
         {
            prop_ALA = (float) -0.17;
            prop_TOT = prop_TOT + prop_ALA;
         }
         if (res_nm.equals("ARG")==true)
         {
            prop_ARG = (float)  0.59;
            prop_TOT = prop_TOT + prop_ARG;
         }
         if (res_nm.equals("ASN")==true)
         {
            prop_ASN = (float)  0.16;
            prop_TOT = prop_TOT + prop_ASN;
         }
         if (res_nm.equals("ASP")==true)
         {
           prop_ASP = (float) -1.21;
           prop_TOT = prop_TOT + prop_ASP;
         }
         if (res_nm.equals("CYS")==true)
         {
            prop_CYS = (float) 0.10;
            prop_TOT = prop_TOT + prop_CYS;
         }
         if (res_nm.equals("GLN")==true)
         {
            prop_GLN = (float) -0.12;
            prop_TOT = prop_TOT + prop_GLN;
         }
         if (res_nm.equals("GLU")==true)
         {
            prop_GLU = (float) -1.46;
            prop_TOT = prop_TOT + prop_GLU;
         }
         if (res_nm.equals("GLY")==true)
         {
            prop_GLY = (float)  0.23;
            prop_TOT = prop_TOT + prop_GLY;
         }
         if (res_nm.equals("HIS")==true)
         {
            prop_HIS = (float) 0.17;
            prop_TOT = prop_TOT + prop_HIS;
         }
         if (res_nm.equals("ILE")==true)
         {
            prop_ILE = (float) -0.27;
            prop_TOT = prop_TOT + prop_ILE;
         }
         if (res_nm.equals("LEU")==true)
         {
            prop_LEU = (float) -0.71;
            prop_TOT = prop_TOT + prop_LEU;
         }
         if (res_nm.equals("LYS")==true)
         {
            prop_LYS = (float) 0.06;
            prop_TOT = prop_TOT + prop_LYS;
         }
         if (res_nm.equals("MET")==true)
         {
            prop_MET = (float) -0.06;
            prop_TOT = prop_TOT + prop_MET;
         }
         if (res_nm.equals("PHE")==true)
         {
            prop_PHE = (float) 0.36;
            prop_TOT = prop_TOT + prop_PHE;
         }
         if (res_nm.equals("PRO")==true)
         {
            prop_PRO = (float) -0.59;
            prop_TOT = prop_TOT + prop_PRO;
         }
        if (res_nm.equals("SER")==true)
        {
           prop_SER = (float)  0.26;
           prop_TOT = prop_TOT + prop_SER;
        }
        if (res_nm.equals("THR")==true)
        {
           prop_THR = (float)  0.27;
           prop_TOT = prop_TOT + prop_THR;
        }
        if (res_nm.equals("TRP")==true)
        {
           prop_TRP = (float) 0.43;
           prop_TOT = prop_TOT + prop_TRP;
        }
        if (res_nm.equals("TYR")==true)
        {
           prop_TYR = (float) 0.47;
           prop_TOT = prop_TOT + prop_TYR;
        }
        if (res_nm.equals("VAL")==true)
        {
           prop_VAL = (float) -0.17;
           prop_TOT = prop_TOT + prop_VAL;
        }
     }
   }   // end of if
   }   // end of while

   in.close();

   if (("" + rnd.format(prop_TOT)).length() == 4)
      s1 = "      " + rnd.format(prop_TOT);
   else if (("" + rnd.format(prop_TOT)).length() == 5)
      s1 = "     " + rnd.format(prop_TOT);

   System.out.println("Residue Propensity Score = " + rnd.format(prop_TOT));

   out.writeBytes(s1 + "\n");

   out.close();

   } // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   } // main end

} // class propensity end