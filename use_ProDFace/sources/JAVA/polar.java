/* ***************************************************************************
      PROG. TO FIND OUT THE POLAR AND APOLAR FRACTION OF AN INTERFACE IN
      PROTEIN-PROTEIN COMPLEX

      CALCULATES THE POLAR AND APOLAR FRACTION FOR SIDE CHAIN AND BACKBONE
      ALSO CALCULATES NUMBER OF POLAR AND APOLAR ATOMS PRESENT IN THE .int FILES

      MODIFIED FOR FIND OUT THE % OF POLAR/APOLAR FOR MC/MC MC/SC SC/MC SC/SC

      I/P FILE  .int file(1a3ca.int)
      O/P FILE  .out

      # Authors : Arumay PAL                                          |
      # Contact : contact.arumay@gmail.com
**************************************************************************** */

import java.text.*;
import java.io.*;

class polar
{
   public static void main(String args[])
   {
   RandomAccessFile in, out;
   String path, pdb, pdb_fil1, pdb_fil2;
   String s, s1 = "", s2 = "", s3 = "";
   String atom_nm, res_nm;
   int c1, res_no, polarcount, apolarcount, pbbcount, psccount, apobbcount, aposccount;
   float polar_bb, polar_sc, apolar_bb, apolar_sc;
   float asa_su, asa_cmp, bsa, polar, apolar, tot_int_area, percent_polar, percent_apolar, percent_polar_atoms, percent_apolar_atoms;
   float percent_polar_bb, percent_polar_sc, percent_apolar_bb, percent_apolar_sc;
   float percent_polar_bb_no, percent_polar_sc_no, percent_apolar_bb_no, percent_apolar_sc_no;
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

   c1 = 0;
   bsa = 0;
   polar_bb = 0;
   pbbcount = 0;
   polar_sc = 0;
   psccount = 0;
   apolar_bb = 0;
   apobbcount = 0;
   apolar_sc = 0;
   aposccount = 0;
   polar = 0;
   apolar = 0;

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
      atom_nm = s.substring(13,16);
      res_nm = s.substring(17,20);
      res_no = Integer.parseInt(s.substring(22,26).trim());
      asa_su = Float.parseFloat(s.substring(66,72).trim());
      asa_cmp = Float.parseFloat(s.substring(72,78).trim());
      bsa = asa_su - asa_cmp;

      if ((atom_nm.equals("O  ")==true)  || (atom_nm.equals("N  ")==true))
      {
         polar_bb = polar_bb + bsa;
         pbbcount++;
      }
      // For DNA, added on 16.02.10
      if ((atom_nm.equals("O5'")==true) || (atom_nm.equals("O4'")==true) || (atom_nm.equals("O3'")==true) ||
          (atom_nm.equals("N9 ")==true) || (atom_nm.equals("N7 ")==true) || (atom_nm.equals("N6 ")==true) ||
          (atom_nm.equals("N4 ")==true) || (atom_nm.equals("N3 ")==true) || (atom_nm.equals("N2 ")==true) ||
          (atom_nm.equals("N1 ")==true) || (atom_nm.equals("P  ")==true) || (atom_nm.equals("OP1")==true) ||
          (atom_nm.equals("OP2")==true) || (atom_nm.equals("O2 ")==true) || (atom_nm.equals("O4 ")==true) ||
          (atom_nm.equals("O6 ")==true))
      {
         polar_bb = polar_bb + bsa;
         pbbcount++;
      }

      // --------------------------
      if ((atom_nm.equals("OD1")==true) || (atom_nm.equals("OD2")==true) || (atom_nm.equals("OE1")==true) ||
          (atom_nm.equals("OE2")==true) || (atom_nm.equals("OG ")==true) || (atom_nm.equals("OG1")==true) ||
          (atom_nm.equals("OH ")==true) || (atom_nm.equals("NE ")==true) || (atom_nm.equals("NE1")==true) ||
          (atom_nm.equals("NE2")==true) || (atom_nm.equals("ND1")==true) || (atom_nm.equals("ND2")==true) ||
          (atom_nm.equals("NZ ")==true) || (atom_nm.equals("NH1")==true) || (atom_nm.equals("NH2")==true) ||
          (atom_nm.equals("SD ")==true) || (atom_nm.equals("SG ")==true))
      {
         polar_sc = polar_sc + bsa;
         psccount++;
      }

      // --------------------------------------------------------------------------

      if ((atom_nm.equals("C  ")==true) || (atom_nm.equals("CA ")==true))
      {
         apolar_bb = apolar_bb + bsa;
         apobbcount++;
      }
      // For DNA, added on 16.02.10
      if ((atom_nm.equals("C1'")==true) || (atom_nm.equals("C2'")==true) || (atom_nm.equals("C3'")==true) ||
          (atom_nm.equals("C4'")==true) || (atom_nm.equals("C5'")==true) || (atom_nm.equals("C2 ")==true) ||
          (atom_nm.equals("C4 ")==true) || (atom_nm.equals("C5 ")==true) || (atom_nm.equals("C6 ")==true) ||
          (atom_nm.equals("C7 ")==true) || (atom_nm.equals("C8 ")==true))
      {
         apolar_bb = apolar_bb + bsa;
         apobbcount++;
      }
      // --------------------------
      if ((atom_nm.equals("CB ")==true) || (atom_nm.equals("CD ")==true) || (atom_nm.equals("CD1")==true) ||
          (atom_nm.equals("CD2")==true) || (atom_nm.equals("CG ")==true) || (atom_nm.equals("CG1")==true) ||
          (atom_nm.equals("CG2")==true) || (atom_nm.equals("CE ")==true) || (atom_nm.equals("CE1")==true) ||
          (atom_nm.equals("CE2")==true) || (atom_nm.equals("CE3")==true) || (atom_nm.equals("CZ ")==true) ||
          (atom_nm.equals("CZ2")==true) || (atom_nm.equals("CZ3")==true) || (atom_nm.equals("CH3")==true))
      {
         apolar_sc = apolar_sc + bsa;
         aposccount++;
      }
   c1++;
   }   // end of if
   }   // end of while

   in.close();

   polar = polar_bb + polar_sc;
   apolar = apolar_bb + apolar_sc;

   polarcount = pbbcount + psccount;
   apolarcount = apobbcount + aposccount;

   percent_polar_bb = polar_bb/polar;
   percent_polar_sc = polar_sc/polar;
   percent_apolar_bb = apolar_bb/apolar;
   percent_apolar_sc = apolar_sc/apolar;

   percent_polar_bb_no = (float) (pbbcount/((pbbcount+psccount)*1.0));
   percent_polar_sc_no = (float) (psccount/((pbbcount+psccount)*1.0));
   percent_apolar_bb_no = (float) (apobbcount/((apobbcount+aposccount)*1.0));
   percent_apolar_sc_no = (float) (aposccount/((apobbcount+aposccount)*1.0));

   tot_int_area = polar + apolar;
   percent_polar = polar/tot_int_area;
   percent_apolar = apolar/tot_int_area;
   percent_polar_atoms = (float) (polarcount/(c1*1.0));
   percent_apolar_atoms = (float) (apolarcount/(c1*1.0));

   if (("" + rnd.format(tot_int_area)).length() == 4)
      s1 = "      " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 5)
      s1 = "     " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 6)
      s1 = "    " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 7)
      s1 = "   " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 8)
      s1 = "  " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 9)
      s1 = " " + rnd.format(tot_int_area);
   else if (("" + rnd.format(tot_int_area)).length() == 10)
      s1 = "" + rnd.format(tot_int_area);

   s2 = "      " + rnd.format(percent_apolar_atoms);

   if (("" + rnd.format(apolar)).length() == 4)
      s3 = "      " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 5)
      s3 = "     " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 6)
      s3 = "    " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 7)
      s3 = "   " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 8)
      s3 = "  " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 9)
      s3 = " " + rnd.format(apolar);
   else if (("" + rnd.format(apolar)).length() == 10)
      s3 = "" + rnd.format(apolar);

   System.out.println("Total area in the interface = " + rnd.format(tot_int_area));
   System.out.println("Fraction of non-polar atoms in the interface = " + rnd.format(percent_apolar_atoms));
   System.out.println("Non-polar area in the interface = " + rnd.format(apolar));

   out.writeBytes(s1 + s2 + s3 + "\n");

   out.close();

   } // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   } // main end

} // class polar end