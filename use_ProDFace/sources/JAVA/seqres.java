/* ********************************************************************
   CREATES THE FILE "seqres1A3Q.txt" CONTAINING RESIDUE AND SECONDARY
   STRUCTURE INFO

    I/P :: 1A3Q.dssp
    O/P :: seqres1A3Q.txt

    # Authors : Arumay PAL                                          |
    # Contact : contact.arumay@gmail.com
*********************************************************************** */

/*
 *
 * Creates the file "seqrestemp.txt" containing residue and secondary structure
 * information for displaying via the various plots later on.
 *
 * Take care that the only regular console output statements in seqres.java
 * are for the error trapping conditions handled in upload.php
 *
 * @version J2SDK 1.4.2_04
 * @author  Saptarshi Mandal
 * @date    16.4.2004
 *
 */

import java.io.*;

class seqres
{
   public static void main(String args[])
   {
   RandomAccessFile in;
   RandomAccessFile out;
   String s = "",s2 = "",s3 = "",s4 = "",s5 = "",s6 = "";
   String path, pdb;
   String chainsfile;
   String dsspfile;
   String seqresfile;
   String name[] = new String[20];
   int rescnt[] = new int[20];
   String totres[] = new String[20];
   int i,j;
   int rstart[] = new int[20];
   int gap[] = new int[20];
   int hstart = 0,hstop = 0;
   int sstart = 0,sstop = 0;
   int cnt = 0;
   int cnt2 = 0;
   int cnt3 = 0;
   int cnt4 = 0;
   int flag1 = 0;
   int flag2 = 0;
   int flag3 = 0;
   int chainflag = 0;
   int chcnt = 0;
   int store[] = new int[2];
   int diff = 1;

   if (args.length != 1)
   {
      System.out.println("Please pass the dssp filename as argument. Halting Program.");
      return;
   }

   path = args[0].substring(0,20);
   pdb = args[0].substring(20,24);

   dsspfile = pdb + ".dssp";
   chainsfile = "chains" + pdb + ".txt";
   seqresfile = "seqres" + pdb + ".txt";

   try
   {

   in = new RandomAccessFile(path + chainsfile,"r");

   while ((s = in.readLine()) != null)
   {
      chainflag = 1;

      if (cnt == 1)   // 2nd line contains no. of chains
      {
         chcnt =  Integer.parseInt(s);
      }
      else if (cnt > 1)   // 3rd line onwards contain sub-unit ids
      {
         name[cnt2++] = s.toUpperCase();
         // System.out.println(name[cnt2-1]);
      }
      cnt++;
   }

   if (chainflag == 0)   // useful for detecting whether no valid information is present in dssp file in cases like 1bdx.dssp
   {
      System.out.println("No valid information present in dssp file. Halting Program.");
      return;
   }

   in.close();
   // System.out.println(chcnt);

   for (i=0;i<chcnt;i++)
   {
      in = new RandomAccessFile(path + dsspfile,"r");

      do
      {
      s = in.readLine();
      } while (s.startsWith("  #  RESIDUE AA STRUCTURE") == false);   // positions the cursor for reading secondary structure information

      while ((s = in.readLine()) != null)
      {

      if (name[i].equals(" "))   // useful for detecting whether no chains are present in dssp file in cases like 1bds.dssp
      {
         System.out.println("No chains present in dssp file. Halting Program.");
         return;
      }

      if (s.substring(11,12).toUpperCase().equals(name[i]))
      {
         rescnt[i] = Integer.parseInt(s.substring(6,10).trim());
         totres[i] = s.substring(6,10);
      }
      }
      in.close();
   }

   out = new RandomAccessFile(path + seqresfile,"rw");
   out.setLength(0);

   for (i=0;i<chcnt;i++)
   {
      in = new RandomAccessFile(path + dsspfile,"r");
      flag3 = 0;
      cnt4 = 0;

      do
      {
      s = in.readLine();
      } while (s.startsWith("  #  RESIDUE AA STRUCTURE") == false);   // positions the cursor for reading secondary structure information

      while ((s = in.readLine()) != null)
      {
      if (s.substring(11,12).toUpperCase().equals(name[i]))
      {
         if (flag3 == 0)
         {
            rstart[i] = Integer.parseInt(s.substring(6,10).trim());
            gap[i] = rstart[i] - 1;
            if (gap[i] < 0)   // useful for detecting non-standard residue line numbering in cases like 1bdm.dssp
            {
               System.out.println("Non-standard line numbering of residues in dssp file. Halting Program.");
               return;
            }
            flag3 = 1;
         }

         store[cnt4%2] = Integer.parseInt(s.substring(6,10).trim());

         if (cnt4%2 == 0)
            diff = store[0] - store[1];   // store[0] contains larger among the two consecutive numbers for even cases
         else
            diff = store[1] - store[0];   // store[1] contains larger among the two consecutive numbers for odd cases

         if (cnt4 == 0)
         {
            s2 += s.substring(13,14).toUpperCase();
            cnt3++;
         }
         else
         {

         if (cnt4%2 == 0)   // useful for detecting non-standard residue line numbering in cases like 1bd7.dssp
         {
            if (store[0] <= store[1])   // store[0] should be > store[1] for even cases
            {
               System.out.println("Non-standard line numbering of residues in dssp file. Halting Program.");
               return;
            }
         }
         else
         {
            if (store[1] <= store[0])   // store[1] should be > store[0] for odd cases
            {
               System.out.println("Non-standard line numbering of residues in dssp file. Halting Program.");
               return;
            }
         }

         if (diff == 1)
         {
            s2 += s.substring(13,14).toUpperCase();
            cnt3++;
         }
         else
         {
            for (j = 0; j < diff-1; j++)
            {
               s2 += " ";
               cnt3++;
               if (cnt3%13 == 0)   // 13 residues per line
               {
                  // "gap[i]" denotes the difference between residue numbering startings in dssp and pdb files

                  if (("" + gap[i]).length() == 1)
                     s5 = "   " + gap[i];
                  else if (("" + gap[i]).length() == 2)
                     s5 = "  " + gap[i];
                  else if (("" + gap[i]).length() == 3)
                     s5 = " " + gap[i];
                  else if (("" + gap[i]).length() == 4)
                     s5 = "" + gap[i];

                  // System.out.println("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s5 + " ");
                  out.writeBytes("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s5 + " " + "\n");
                  s2 = "";
                  cnt3 = 0;
               }   // end of if
            }   // end of for
            s2 += s.substring(13,14).toUpperCase();
            cnt3++;
         }

         }

         cnt4++;

         if (cnt3%13 == 0)   // 13 residues per line
         {
            // "gap[i]" denotes the difference between residue numbering startings in dssp and pdb files

            if (("" + gap[i]).length() == 1)
               s5 = "   " + gap[i];
            else if (("" + gap[i]).length() == 2)
               s5 = "  " + gap[i];
            else if (("" + gap[i]).length() == 3)
               s5 = " " + gap[i];
            else if (("" + gap[i]).length() == 4)
               s5 = "" + gap[i];

            // System.out.println("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s5 + " ");
            out.writeBytes("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s5 + " " + "\n");
            s2 = "";
            cnt3 = 0;
         }
      }
      }

      if (s2.equals("") == false)   // last line contains less than 13 residues
      {
         int left = 13 - s2.length();
         for (j=0; j<left; j++)   // fill up remaining space of last line with spaces
         {
            s2 += " ";
         }

         // "gap[i]" denotes the difference between residue numbering startings in dssp and pdb files

         if (("" + gap[i]).length() == 1)
            s6 = "   " + gap[i];
         else if (("" + gap[i]).length() == 2)
            s6 = "  " + gap[i];
         else if (("" + gap[i]).length() == 3)
            s6 = " " + gap[i];
         else if (("" + gap[i]).length() == 4)
            s6 = "" + gap[i];

         // System.out.println("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s6 + " ");
         out.writeBytes("SEQRES " + totres[i] + " " + name[i]+ "  " + s2 + " " + s6 + " " + "\n");
         s2 = "";
         cnt3 = 0;
      }

      in.close();
   }

   for (i=0;i<chcnt;i++)
   {
      in = new RandomAccessFile(path + dsspfile,"r");

      do
      {
      s = in.readLine();
      } while (s.startsWith("  #  RESIDUE AA STRUCTURE") == false);   // positions the cursor for reading secondary structure information

      while ((s = in.readLine()) != null)
      {
      if ((s.substring(11,12).toUpperCase().equals(name[i])) && (s.substring(16,17).toUpperCase().equals("H")))
      {
         if (flag1 == 0)
         {
         hstart = Integer.parseInt(s.substring(6,10).trim()) - gap[i];   // helix start as appears on screen
         flag1 = 1;
         }
         else
         {
         hstop = Integer.parseInt(s.substring(6,10).trim()) - gap[i];   // helix stop as appears on screen
         }
      }
      else if (flag1 == 1)
      {
         if (("" + hstart).length() == 1)
            s3 = "   " + hstart;
         else if (("" + hstart).length() == 2)
            s3 = "  " + hstart;
         else if (("" + hstart).length() == 3)
            s3 = " " + hstart;
         else if (("" + hstart).length() == 4)
            s3 = "" + hstart;

         if (("" + hstop).length() == 1)
            s4 = "   " + hstop;
         else if (("" + hstop).length() == 2)
            s4 = "  " + hstop;
         else if (("" + hstop).length() == 3)
            s4 = " " + hstop;
         else if (("" + hstop).length() == 4)
            s4 = "" + hstop;

         // System.out.println("HELIX " + name[i] + " " + s3 + " " + s4);
         out.writeBytes("HELIX " + name[i] + " " + s3 + " " + s4 + "\n");
         flag1 = 0;
      }
      }
      in.close();
   }

   for (i=0;i<chcnt;i++)
   {
      in = new RandomAccessFile(path + dsspfile,"r");

      do
      {
      s = in.readLine();
      } while (s.startsWith("  #  RESIDUE AA STRUCTURE") == false);   // positions the cursor for reading secondary structure information

      while ((s = in.readLine()) != null)
      {
      if ((s.substring(11,12).toUpperCase().equals(name[i])) && (s.substring(16,17).toUpperCase().equals("E")))
      {
         if (flag2 == 0)
         {
         sstart = Integer.parseInt(s.substring(6,10).trim()) - gap[i];   // sheet start as appears on screen
         flag2 = 1;
         }
         else
         {
         sstop = Integer.parseInt(s.substring(6,10).trim()) - gap[i];   // sheet stop as appears on screen
         }
      }
      else if (flag2 == 1)
      {
         if (("" + sstart).length() == 1)
            s3 = "   " + sstart;
         else if (("" + sstart).length() == 2)
            s3 = "  " + sstart;
         else if (("" + sstart).length() == 3)
            s3 = " " + sstart;
         else if (("" + sstart).length() == 4)
            s3 = "" + sstart;

         if (("" + sstop).length() == 1)
            s4 = "   " + sstop;
         else if (("" + sstop).length() == 2)
            s4 = "  " + sstop;
         else if (("" + sstop).length() == 3)
            s4 = " " + sstop;
         else if (("" + sstop).length() == 4)
            s4 = "" + sstop;

         // System.out.println("SHEET " + name[i] + " " + s3 + " " + s4);
         out.writeBytes("SHEET " + name[i] + " " + s3 + " " + s4 + "\n");
         flag2 = 0;
      }
      }
      in.close();
   }

   out.close();

   }   // end of try

   catch (Exception e)
   {
      System.out.println(e.toString());
   }

   }   // main end

} // class seqres end