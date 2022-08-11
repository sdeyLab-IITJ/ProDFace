********************************************************************
********************************************************************
|                                                                  |
|                           ProDFace                               |
|                                                                  |
|         Web-tool for the analysis of Protein-DNA structures      |
|         (experimental structures, models, MD snapshots,          |
|         docked structures etc.)                                  |
|                                                                  |
|                             by                                   |
|        Arumay PAL, Pinak Chakrabarti, Sucharita Dey              |
|                                                                  |
********************************************************************									
******************************************************************** 


1. INTRODUCTION
----------------
   ProDFace is a web tool that characterizes the binding region of a protein-DNA complex
   based on amino acid propensity, hydrogen bond (HB) donor capacity (number of solvent accessible
   HB donor groups), sequence conservation at the interface core and rim region, and geometry.
   The program takes as input the structure of a protein-DNA complex in PDB (Protein Data Bank) format,
   and outputs various physicochemical and geometric parameters of the interface, as well as conservation
   of the interface residues in the protein component.

2. SYSTEM REQUIREMENT
---------------------
   ProDFace is a web tool that characterizes the binding region of a protein-DNA complex based on
   amino acid propensity, hydrogen bond (HB) donor capacity (number of solvent accessible HB donor groups),
   sequence conservation at the interface core and rim region, and geometry. The program takes as input the
   structure of a protein-DNA complex in PDB (Protein Data Bank) format, and outputs various physicochemical
   and geometric parameters of the interface, as well as conservation of the interface residues in the protein
   component.

   To run different ProDFace components that calculate important parameters, following softwares are required-
      i) JDK (Java Development Kit)
     ii) C compiler
    iii) Perl Compiler (Interpreter)

3. EXTERNAL PROGRAMS
---------------------
   ProDFace uses NACCESS for accessible surface area calculation, DSSP for secondary structure assignments,
   HBPLUS for hydrogen bond calculation, and PRINCIP and TRANSFORM module from SURFNET to transform the
   3D coordinates of the molecules along their principal axes (2D).

4. STEPS FOR RUNNING ProDFace web tool
---------------------------------------   
   Step 1
   Provide the protein-DNA complex structure by
   (a) providing the PDB code (for experimental structures)
                              OR
   (b) uploading your own complex structure (for modelled structures, docking solutions, or MD snapshots) in PDB format.

   Step 2
   Provide chain ids of the protein and the DNA components that form the interface to be analyzed.

   Step 3
   Provide the threshold distance for clustering interface atoms into spatial patches (defalut 20 Ang.)

   Step 4
   Provide Multiple Sequence Alignment (MSA) for the protein component to calculate residue conservation scores by
   (a) Using HSSP database
              OR
   (b) Providing user-defined MSA in FASTA format

5. IMPORTANT ProDFace components
---------------------------------
   C programs
   ----------
     i) InterfaceConsvSeqIdenArea.c
        Calculates conservation score (area normalized shannon entropy) of interface residues
    ii) SubunitAvgEntropy_c1.c
        Calculates 'average conservation score' of a protein subunit
   iii) SecStrClass.c
        Calculates percentage of SSE components (helix, strand, Non-regular) of a protein subunit
    iv) PPSegment.c
        Calculates 'number of interface segments' of a subunit and length of each segment
     v) SSSegment.c
        Determines SSE of each interface segment

    Java programs
    -------------
     i) intDist.java
        Calculates distances between all possible interface atom pairs
    ii) polar.java
        Finds the polar and apolar fraction of an interface
   iii) buried.java
        Finds fraction of total buried, polar buried, and non-polar buried atoms in an interface
    iv) seqres.java
        Finds DSSP assigned secondary structure of each interface residue
     v) propensity.java
        Calculates "Residue Propensity Score" of an interface
    vi) localDensity.java
        Calculates 'Local Density' of the interface

    Perl programs
    -------------
     i) int_water.pl
        Calculates Number of interface Waters
    ii) Whbs_one_sided.pl
        Calculates total No. of water mediated hydrogen bonds formed by an interface 
   iii) BW_distribution.pl
        Calculates total number of Bridging-water and a distribution of total water mediated Hbonds
    iv) int_HB_donor_atm_asa.pl
        Caclculates accesibility of Hbond donor groups at the interface


If you have any problem or suggestion, please contact:
contact.arumay@gmail.com OR sdey@iitj.ac.in
