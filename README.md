# ProDFace
### Web-tool for the analysis of Protein-DNA structures (experimental structures, models, MD snapshots, docked structures etc.)

#### 1. INTRODUCTION  
[ProDFace](http://structbioinfo.iitj.ac.in/resources/bioinfo/pd_interface/) is a web tool that characterizes the binding region of a protein-DNA complex based on amino acid propensity, hydrogen bond (HB) donor capacity (number of solvent accessible HB donor groups), sequence conservation at the interface core and rim region, and geometry. The program takes as input the structure of a protein-DNA complex in PDB (Protein Data Bank) format, and outputs various physicochemical and geometric parameters of the interface, as well as conservation of the interface residues in the protein component.

#### 2. SYSTEM REQUIREMENT  
The web version of ProDFace has been tested in Chrome, Firefox and Edge. However, Java plug-ins should be pre-installed to display Java applets properly. The stand alone version of the program will be released in the next update of ProDFace.

#### 3. EXTERNAL PROGRAMS
ProDFace uses [NACCESS](http://www.bioinf.manchester.ac.uk/naccess/) for accessible surface area calculation, [DSSP](https://swift.cmbi.umcn.nl/gv/dssp/) for secondary structure assignments, [HBPLUS](https://www.ebi.ac.uk/thornton-srv/software/HBPLUS/) for hydrogen bond calculation, and PRINCIP and TRANSFORM module from [SURFNET](https://www.ebi.ac.uk/thornton-srv/software/SURFNET/) to transform the 3D coordinates of the molecules along their principal axes (2D).

#### 4. STEPS FOR RUNNING ProDFace (web version)
##### Step 1 #####
a) Provide the protein-DNA complex structure by EITHER a) providing the PDB code (for experimental structures), OR b) uploading your own complex structure (for modelled structures, docking solutions, or MD snapshots) in PDB format.
##### Step 2 #####
Provide chain ids of the protein and the DNA components that form the interface to be analyzed.
##### Step 3 #####
Provide the threshold distance for clustering interface atoms into spatial patches (defalut 20 Ang.)
##### Step 4 #####
Provide Multiple Sequence Alignment (MSA) for the protein component to calculate residue consrvation scores.  
a) Using HSSP database  
           or   
b) Providing user-defined MSA in FASTA format  
