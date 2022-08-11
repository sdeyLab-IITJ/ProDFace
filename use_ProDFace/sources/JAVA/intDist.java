/* ********************************************************************
   PROG TO CALCULATE DISTANCES BETWEEN ALL POSSIBLE ATOM PAIRS

    I/P :: 1A3QA.int
    O/P :: 1A3QA.d, dist1A3Q.dat

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
*********************************************************************** */

import java.text.*;
import java.io.*;

class intDist
{
   public static void main(String args[])
   {
   RandomAccessFile in;
   RandomAccessFile out;
   RandomAccessFile out2;
   String path, pdb, su;
   String s, s2 = "", s3 = "", s4 = "";
   String pdb_fil, pdb_fil2, pdb_fil3;
   int i, j;
   int index = 0;
   int natoms = 0;
   int n_sl[] = new int[1500];   //
   int n_res[] = new int[1500];   //
   double atm[] = new double[1500];   //
   double atm2[] = new double[1500];   //
   double atm3[] = new double[1500];   //
   double dist, dmax = -999, dmin = 999;
   String dist2;
   DecimalFormat rnd = new DecimalFormat("0.0");

   if (args.length != 2)
   {
      System.out.println("Please pass the pdb filename and sub-unit id as argument. Halting Program.");
      return;
   }

   path = args[0].substring(0,20);
   pdb = args[0].substring(20,24);
   su = args[1];
   pdb_fil = pdb + su + ".int";
   pdb_fil2 = pdb + su + ".d";
   pdb_fil3 = "dist" + pdb + ".dat";

   try
   {

   in = new RandomAccessFile (path + pdb_fil,"r");
   out = new RandomAccessFile(path + pdb_fil2,"rw");
   out.setLength(0);
   out2 = new RandomAccessFile(path + pdb_fil3,"rw");
   out2.seek(out2.length());

   while ((s = in.readLine()) != null)    // start of while
   {
   if ((s.startsWith("ATOM") == true) || (s.startsWith("HETATM") == true))   // start of if
   {
      n_sl[natoms] = Integer.parseInt(s.substring(7,11).trim());
      n_res[natoms] = Integer.parseInt(s.substring(22,26).trim());
      atm[natoms] = Double.parseDouble(s.substring(30,38).trim());
      atm2[natoms] = Double.parseDouble(s.substring(38,46).trim());
      atm3[natoms] = Double.parseDouble(s.substring(46,54).trim());
      natoms = natoms + 1;
   }   // end of if
   }   // end of while
   in.close();
   //System.out.println("Total no. of atoms in the input file = " + natoms);

   for (i=0;i<natoms;i++)   // start of for i
   {
   for (j=i+1;j<natoms;j++)   // start of for j
	{
   s2 = "";
	dist = Math.pow((atm[i]-atm[j]),2.0) + Math.pow((atm2[i]-atm2[j]),2.0) + Math.pow((atm3[i]-atm3[j]),2.0);
   dist = Math.pow((dist),0.5);
   if (dist > dmax)
      dmax = dist;
   if (dist < dmin)
      dmin = dist;
   switch((""+n_sl[i]).length())
   {
   case 1:
      s2 += "   " + n_sl[i] + " ";
      break;
   case 2:
      s2 += "  " + n_sl[i] + " ";
      break;
   case 3:
      s2 += " " + n_sl[i] + " ";
      break;
   case 4:
      s2 += n_sl[i] + " ";
      break;
   }

   switch((""+n_sl[j]).length())
   {
   case 1:
      s2 += "   " + n_sl[j] + " ";
      break;
   case 2:
      s2 += "  " + n_sl[j] + " ";
      break;
   case 3:
      s2 += " " + n_sl[j] + " ";
      break;
   case 4:
      s2 += n_sl[j] + " ";
      break;
   }

   dist2 = rnd.format(dist);

   switch(dist2.length())
   {
   case 1:
      s2 += "   " + dist2 + " ";
      break;
   case 2:
      s2 += "  " + dist2 + " ";
      break;
   case 3:
      s2 += " " + dist2 + " ";
      break;
   case 4:
      s2 += dist2 + " ";
      break;
   }

   switch((""+n_res[i]).length())
   {
   case 1:
      s2 += "  " + n_res[i] + " ";
      break;
   case 2:
      s2 += " " + n_res[i] + " ";
      break;
   case 3:
      s2 += n_res[i] + " ";
      break;
   }

   switch((""+n_res[j]).length())
   {
   case 1:
      s2 += "  " + n_res[j] + " ";
      break;
   case 2:
      s2 += " " + n_res[j] + " ";
      break;
   case 3:
      s2 += n_res[j] + " ";
      break;
   }
   index++;
   //System.out.println(s2);
   out.writeBytes(s2+"\n");
   }   // end of for i
   }   // end of for j

   if (("" + rnd.format(dmax)).length() == 3)
      s3 = "       " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 4)
      s3 = "      " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 5)
      s3 = "     " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 6)
      s3 = "    " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 7)
      s3 = "   " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 8)
      s3 = "  " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 9)
      s3 = " " + rnd.format(dmax);
   else if (("" + rnd.format(dmax)).length() == 10)
      s3 = "" + rnd.format(dmax);

   if (("" + rnd.format(dmin)).length() == 3)
      s4 = "       " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 4)
      s4 = "      " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 5)
      s4 = "     " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 6)
      s4 = "    " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 7)
      s4 = "   " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 8)
      s4 = "  " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 9)
      s4 = " " + rnd.format(dmin);
   else if (("" + rnd.format(dmin)).length() == 10)
      s4 = "" + rnd.format(dmin);

   //System.out.println(s3 + s4);
   out2.writeBytes(s3 + s4 +"\n");

   out.close();
   out2.close();
   } // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   } // main end

} // class intdist end