/* ********************************************************************
   PROG TO CALCULATE 'LOCAL DENSITY' OF INTERFACE ATOMS
   [algo- No. of atoms within 12A of each atom]

    I/P :: 1A3Q_protein.int
    O/P :: stat1A3Q.dat

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
*********************************************************************** */

import java.text.*;
import java.math.*;
import java.io.*;

class localDensity
{

   static int STORE = 8000;

   public static void main(String args[])
   {
   RandomAccessFile in, out;
   String path, pdb, pdb_fil1, pdb_fil2;
   String s, s1 = "";
   String atom_name[] = new String[STORE], res_name[] = new String[STORE], ch_id[] = new String[STORE], x_char[] = new String[STORE];
   int c1, c2, i, j, atom_num[] = new int[STORE], res_no[] = new int[STORE], store_partner[] = new int[STORE], temp;
   float fcoord_x[] = new float[STORE], fcoord_y[] = new float[STORE], fcoord_z[] = new float[STORE];
   float occupancy[] = new float[STORE], temp_fact[] = new float[STORE], d, avg, dev, std, sum;
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
       atom_num[c1] = Integer.parseInt(s.substring(7,11).trim());
       atom_name[c1] = s.substring(13,16);
       res_name[c1] = s.substring(17,20);
       ch_id[c1] = s.substring(21,22);
       res_no[c1] = Integer.parseInt(s.substring(22,26).trim());
       x_char[c1] = s.substring(26,27);
       fcoord_x[c1] = Float.parseFloat(s.substring(30,38).trim());
       fcoord_y[c1] = Float.parseFloat(s.substring(38,46).trim());
       fcoord_z[c1] = Float.parseFloat(s.substring(46,54).trim());
       occupancy[c1] = Float.parseFloat(s.substring(54,60).trim());
       temp_fact[c1] = Float.parseFloat(s.substring(60,66).trim());
       c1++;
   }   // end of if
   }   // end of while

   in.close();

   for (i=0;i<c1;i++)
   {
      c2 = 0;
      for (j=0;j<c1;j++)
      {
         d = (float) Math.sqrt((double)((fcoord_x[j]-fcoord_x[i])*(fcoord_x[j]-fcoord_x[i])+(fcoord_y[j]-fcoord_y[i])*(fcoord_y[j]-fcoord_y[i])+(fcoord_z[j]-fcoord_z[i])*(fcoord_z[j]-fcoord_z[i])));
         if (d<=12.0 && d!=0)
            c2++;
      }
      store_partner[i] = c2;
   }

   // sort the array store_partner

   for(i=0;i<c1-1;i++)
   {
     for(j=i+1;j<c1;j++)
     {
        if (store_partner[i]>store_partner[j])
        {
           temp = store_partner[i];
           store_partner[i] = store_partner[j];
           store_partner[j] = temp;
        }
     }
   }

   // CALCULATION OF AVERAGES PARTNER NUMBERS AND STANDARD DEVIATIONS

   sum = 0;
   avg = 0;
   std = 0;
   dev = 0;

   for (i=0;i<c1;i++)
      sum = sum + store_partner[i];
   avg = (float) (sum/c1);

   for (i=0;i<c1;i++)
      dev = dev + ((store_partner[i]-avg)*(store_partner[i]-avg));
   std = (float) Math.sqrt((double)(dev/c1));

   if (("" + rnd.format(avg)).length() == 4)
      s1 = "      " + rnd.format(avg);
   else if (("" + rnd.format(avg)).length() == 5)
      s1 = "     " + rnd.format(avg);
   else if (("" + rnd.format(avg)).length() == 6)
      s1 = "    " + rnd.format(avg);

   System.out.println("Local Density = " + rnd.format(avg));

   out.writeBytes(s1 + "\n");

   out.close();

   } // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   } // main end

} // class local_density end