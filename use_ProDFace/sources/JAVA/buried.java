/* ********************************************************************************
    PROG. TO FIND THE FRACTION OF TOTAL BURIED ATOMS, FRACTION OF POLAR
    BURIED ATOMS AND FRACTION OF NON-POLAR BURIED ATOMS.

    FRAC_TOT=(TOTAL # OF FULLY BURIED ATOMS/TOTAL # OF ATOMS)IN a.int FILES.
    FRAC_TOT=(ATM_BU/TOTAL#ATOM) (IN THE O/P FILE)

    FRAC_POL=(TOTAL # OF FULLY BURIED POLAR ATOMS/TOTAL # OF POLAR ATOMS)IN
             a.int FILES.
    FRAC_POL=(POL_BU/POL_TOT) (IN THE O/P FILE)

    FRAC_APO=(TOTAL # OF FULLY BURIED APOLAR ATOMS/TOTAL # OF APOLAR ATOMS)IN
              a.int FILES.
    FRAC_APO=(APO_BU/APO_TOT) (IN THE O/P FILE)

    IT TAKES THE I/P FILE (AS 12asa.int) FROM THE FILE list_int (WHICH CONTAINS
    THE LIST OF ALL a.int FILES) & GENERATES THE O/P FILES INT_FRAC_BU

    I/P :: *a.int
    O/P :: INT_FRAC_BU


    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com

******************************************************************************** */

import java.text.*;
import java.io.*;

class buried
{
   public static void main(String args[])
   {
   RandomAccessFile in, out;
   String path, pdb, pdb_fil1, pdb_fil2;
   String s, s1 = "";
   String asa_comp, atom_nm;
   int b1, t1, pol, apol, polB, ApolB;
   float buried;
   DecimalFormat rnd = new DecimalFormat("0.00");

   if (args.length != 1)
   {
      System.out.println("Please pass the .int filename as argument. Halting Program.");
      return;
   }

   path = args[0].substring(0,20);
   pdb = args[0].substring(20,24);
   pdb_fil1 = args[0].substring(20,26) + ".int";
   pdb_fil2 = "stat" + pdb + ".dat";

   b1 = 0;
   t1 = 0;
   pol = 0;
   apol = 0;
   polB = 0;
   ApolB = 0;

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
      asa_comp = s.substring(74,78);

      if (asa_comp.equals("0.00"))
      {
      b1++;
      }

      atom_nm = s.substring(13,14);

      if ((atom_nm.equals("N")==true) || (atom_nm.equals("O")==true) || (atom_nm.equals("S")==true) || (atom_nm.equals("P")==true))
      {
         pol++;
         asa_comp = s.substring(74,78);
         if (asa_comp.equals("0.00"))
            polB++;
      }
      if ((atom_nm.equals("C")==true))
      {
         apol++;
         asa_comp = s.substring(74,78);
         if (asa_comp.equals("0.00"))
            ApolB++;
      }

      t1++;

   }   // end of if
   }   // end of while

   in.close();

   buried = (float) (b1/(t1*1.0));

   s1 = "      " + rnd.format(buried);

   System.out.println("Fraction of fully buried atoms in the interface = " + rnd.format(buried));

   out.writeBytes(s1 + "\n");

   out.close();

   } // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   } // main end

} // class buried end