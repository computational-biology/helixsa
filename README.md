## HelixSA : *A standalone command line tool for alignment of secondary structures of RNA*



##Synopsis  


> HelixSA is a stand alone connand line based linux tool that accepts two structure files in mmCIF/PDB format and aligns their secondary sequence.


## Installation 
>Download the binary executable helixsa.linux from bin directory and bin directory. **DONE**.



## Compilation
> If the executable from downloaded from the release does not work for you, then download the source file given in the release and extract them. go to the directory where the Makefile resides. Run the make command from that directory. The 'helix.linux' will be created and will be stored in bin directory. place it to a suitable directory.
## Dependencies
> The compilation process needs a C and a FORTRAN compiler. We have assumed *gcc* for C and *gfortran* for FORTRAN. We have made the Makefile accordingly. If the user wants a different compiler, change the makefile accordingly. 

## A word about the Makefile
> 
We have given a ready made make file for easy compilation.  The Makefile assumes *gcc* as C compiler and *gfortran* as FORTRAN compiler. If, however, the user do not have access to these compilers and to something else instead, they need a small modification in the Makefile. For example, if you have *clang* as C compiler and *f77* as fortran compiler, then change the following in the Makefile. 
	change `CC := gcc -std=c99` to `CC := clang` and change `FF := gfortran` to `FF := f77`

## How to run
> To run the program use the following command.
`helixsa.linux [OPTIONS] [mmCIF/PDB]  [mmCIF/PDB]`<br>
options are given in detail in the command-line option section of this document.
**EXAMPLE:**<br>
`helixsa.linux 1ehz.cif  1n77.cif`<br>



## Command-line options

>###  Software Related command-line options
>>**-help:** Shows a small help on the terminal.<br> 
**-genhelp:** Generates this help.md file in the current directory. <br>
**-ver:** Prints the version of the program.<br>
**--pub:** Shows how to refer to the software if it is used by anyone.<br>
**--contact:** Prints the contact information for query or bug report.<br>


> ### Domain related command line options.
>> 
**–hbd=[value]:** HelixSA uses 3.8A as default distance for donor-acceptor atom distance for hydrogen bonding interactions. However, this can be changed with this option.  <br><br>
**EXAMPLE** <br>
    `helixsa.linux 1n32.cif -hbdist=4.1`<br><br><br>
**-nsg=[true/false]:** By default MetBP considers pentose sugar O2’ atom for sugar edge based pairing. This is called a sugar mediated base pair. If the user does not want that, then they can set it to ‘false’. <br><br>
**EXAMPLE**<br>
`metbp.linux 1n32.cif -sugmed=false` <br><br><br>
**-nch=[true/false]:** By default MetBP considers C-H...O/N mediated hydrogen bonds for base pair formations. If the user does not want it, the C-H…O/N type hydrogen bonds can be excluded. <br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -chmed=false`<br><br><br>
**-hetatm=[true/false]:** By default, the MetBP software considers all modified nucleic acids. The modified nucleic acid residues can be excluded from computations with the following command line option. <br><br>
**EXAMPLE**<br>
`metbp.linux 1n32.cif -hetatm=false`<br><br><br>
**-chain=name:** The MetBP program can handle a specific chain also. If the user supplies the chain name, the program considers base pairs of only that specific chain and thus metal and base pair binding of that specific chain.<br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -chain=A`<br><br>
**CAUTION !!!**  *If a specific chain is selected then all analyses are restricted to that specific chain only. So, if a base of the given chain has a pair with the bases of the other chain, the program will not report them.* <br><br><br>
**-nmrmdl=[number]:** If an NMR structure is supplied, by default the MetBP software reads the first model and works on that. But it can take a specific model of an NMR structure and can compute the metal and base pair interaction on that model also. <br><br><br>
**EXAMPLE** <br>
`metbp.linux 2ll9.cif -nmrmdl=4`<br><br><br>
**-mode=[bp/nuc/all/dev]:** The MetBP program works in three different modes. The ‘bp’ mode, the ‘nuc’ mode and finally the ‘all’ mode.  In ‘bp’ mode, the details of only those metals are reported which are directly coordinating with an atom of nucleic acid (i.e. base part) that forms a base pair with some other bases. (Note: It may be noted that we consider O2’ of pentose sugar as an important atom for base pairs. So, any bond with O2’ will be treated as an atom of base.) In ‘nuc’ mode, it reports the details of all metals which bind with any atom of a nucleic acid residue. i.e. in this mode, if a metal binds with the oxygen of a phosphate group or any atom of the pentose sugar, the details of that metal is considered. In 'dev' mode, the program runs internally in 'bp' mode but it generates two extra json files. One contains all base pair information and the other contains all metals that has a contact with nucleic acids (The nucleic acis may not form a base pair). <br><br><br>
**-paramfile=[PATH]:** Default the program takes the metal parameters from the ‘metal.params’ file stored in the NUCLEIC_ACID_DIR path location. If the user wants to alter that metal file, she/he may do that. However, if the user does not want to modify the default metal.params file, she may take a copy of the same in any other location and can supply the same by this switch. <br><br>
**EXAMPLE** <br>
`metbp.linux 1n32.cif -paramfile=./metal.params`<br><br>
**Note:**  *If the user supplies the parameter file in this way, then all the parameters will be taken from this file. * <br><br>
**Caution!!!** *This is not an update of the existing values. It is the resetting all the values. For example, if the user supplies the values for Magnesium, then the system will consider only Magnesium and no other metals.* <br><br><br>
**-diff:** With this command if we run the program, files names will be different in different modes. In the MetBP program, there are four files whose contents are changed when modes are changed. These files are .met (metal-base pair info), .det (metal-base pair details), .pml (the pymol script) and _metbp.json (Newly added after the revision). All other files are invariant of mode, as for example, .sum (summary file). So, we have given this mode tag with the file names for only those four files.
However, the default mode will still produce the files names without the mode tag. To get this facility, we need to supply the "-diff" flag as a command-line argument. If we supply -diff, then the files will be different for different modes.
 <br>
**Example:**
	`metbp.linux   1n32.cif      -diff  -mode=bp`
The file will generate files like 1n32_diff_bp.met, 1n32_diff_bp.det etc.  The reason we have added the extra _diff_ in the file name is because the user can easily work with these files when the directory contains files without  -diff mode, like 1n32.met. So, user can easily copy or remove these files using    cp \*_diff_\* or rm \*_diff_\* and can bypass the other files.
<br><br><br>





























##  List of output files
**Let us assume that we have supplied 1n32.cif file to the MetBP program. The following will be the list of output files.**<br>
- **1n32.sum:** This is the summary file. This file reports the number of nucleic acids and protein residues, metals, water molecules found in the structure. Also it gives a summary of every metal binding site, i.e. its coordination, number of nucleic acids, proteins etc to which it binds. <br><br>
- **1n32.met:** This is the main output file. This file gives the metal and base pair interaction details. This file also gives the secondary structure of the base with which the metal binds.<br><br>
- **1n32.det:** This file gives the details of every metal binding site. It gives the metal to the coordinated atom distances, different atom-metal-atom distances.<br><br>
- **1n32.pml:** MetBP generates a [pymol](https://pymol.org/) script. The script, when run through pymol, shows the metal and its different binding residues. To give a better visual effect, the different types of bases are shown in different colors. For example, white bases indicate that they form a base pair. Orange colored base indicates that the  metal binds to any atom of the pentose sugar. Green colored base indicates that the metal binds to the oxygen of the phosphate group.<br><br>
- **1n32_metbp.json:** This file gives the metal-base pair interaction details in JSON format.<br><br>
 - **1n32_basepair.json:** This file gives all the base pair details in JSON format. The get this file the user has to run the program in developer mode . i.e. the program has to be run as
 `metbp.linux 1n32.cif  -mode=dev`<br><br>
 - **1n32_metnuc.json:** This file gives the information of all metals that interact with a nucleic acid resisue whether or not that residue forms a pair. The get this file the user has to run the program in developer mode . i.e. the program has to be run as
 `metbp.linux 1n32.cif  -mode=dev`<br><br>

## Output file format

**Note: We use _atom_site.auth_asym_id, _atom_site.auth_seq_id and _atom_site.pdbx_PDB_ins_code when consider mmCIF files. **
### .met file format
#### BP TAG
          |  Metal Detail       |   Base Pair Details                    |   Outcome                         |
          ----------------------------------------------------------------------------------------------------
           resid  chn  mtl  lnk  loc  res  atm   resid  chn   resid  chn       pair      E-val    dist  numbp
          ----------------------------------------------------------------------------------------------------
	BP      1563  A    MG     1  NUC    G  N7      858  A       828  A       G:A-S:HT    0.22     2.450    1
	BP      1568  A    MG     1  NUC    G  O6     1370  A      1352  A       G:C-W:WC    0.65     2.296    1
	BP      1576  A    MG     2  NUC    G  O6      299  A       566  A       G:G-W:HC    0.16     2.191    1
	BP      1577  A    MG     1  NUC    G  N7      324  A       109  A       G:A-S:HT    0.48     2.360    1

The output has three parts. first the metal detail which gives the residue ID, chain, name of the metal and a link. This link says how may bases does a metal bind.
#### MOTF TAG
          |  Metal Detail  |   Base Pair Details                                                       |
          ----------------------------------------------------------------------------------------------
           resid  chn  mtl sec  bp    type   atm         resid 
          ----------------------------------------------------------------------------------------------
	MOTF                    H   C-G   W:WC 
	MOTF    1563  A   MG->  H   G-A   S:HT   N7           A     858
	MOTF                    H   A-U   W:SC 


	MOTF                    W   C-G   W:WC 
	MOTF    1563  A   MG->  C   G-           N7           A     869
	MOTF                    C   U-

This tag gives another view for the metal- base pair interaction. Here the *sec* tag is indicating the secondary positions. This tag also gives the base pairs above and below the target base pair where the ion binds.

### .det file format
If a metal is bind to a base that forms a base pair, then that metal's other coordination details come under this file. 
#### BIND tag
         |  Metal Detail  |      Coordination atom detail    |
         -----------------------------------------------------
           resid  chn  mtl   loc  res  atm   resid  chn   dist
         -----------------------------------------------------
	BIND    8004  0    MG    NUC  G    O6      456  0    2.077
	BIND    8004  0    MG    PHP  A    OP1     459  0    1.929
	BIND    8004  0    MG    H2O  HOH  O      3305  0    1.930
	BIND    8004  0    MG    H2O  HOH  O      7267  0    2.058
	BIND    8004  0    MG    H2O  HOH  O      8567  0    1.958
	BIND    8004  0    MG    H2O  HOH  O      9300  0    2.050
In this tag, the other binding details of the Magnasium has been shown. 

#### ANGLE tag
          |  Metal Detail  |  C-1 atom detail      |   C-2 atom detail      |      outcome            |
          ---------------------------------------------------------------------------------------------
           RESID  CHN  MTL   resid  chn  res  atm      resid  chn  res  atm     angle (C1-MTL-C2)
          ---------------------------------------------------------------------------------------------
	ANGL    8004  0    MG      456  0    G    O6   :     459  0    A    OP1     87.24
	ANGL    8004  0    MG      456  0    G    O6   :    3305  0    HOH  O       87.66
	ANGL    8004  0    MG      456  0    G    O6   :    7267  0    HOH  O       94.68
	ANGL    8004  0    MG      456  0    G    O6   :    8567  0    HOH  O       85.21
	ANGL    8004  0    MG      456  0    G    O6   :    9300  0    HOH  O      174.00
	ANGL    8004  0    MG      459  0    A    OP1  :    3305  0    HOH  O       90.74
	ANGL    8004  0    MG      459  0    A    OP1  :    7267  0    HOH  O      170.30
	ANGL    8004  0    MG      459  0    A    OP1  :    8567  0    HOH  O       99.02
	ANGL    8004  0    MG      459  0    A    OP1  :    9300  0    HOH  O       90.78
	ANGL    8004  0    MG     3305  0    HOH  O    :    7267  0    HOH  O       79.86
	ANGL    8004  0    MG     3305  0    HOH  O    :    8567  0    HOH  O      167.60
	ANGL    8004  0    MG     3305  0    HOH  O    :    9300  0    HOH  O       98.03
	ANGL    8004  0    MG     7267  0    HOH  O    :    8567  0    HOH  O       90.61
	ANGL    8004  0    MG     7267  0    HOH  O    :    9300  0    HOH  O       88.21
	ANGL    8004  0    MG     8567  0    HOH  O    :    9300  0    HOH  O       89.51
Under this tag, all the coordination atom-metal-atoms are shown. Angles are shown in degree.

### .json files
The software generates a json file for metal-base pair interactions. The software also generates two other json files in developer mode (-mode=dev). Assuming the mmCIF file or PDB fle accn number is 1n32, then the file names will be 1n32_metbp.json for the metal-base pair interaction file. The other two files that are generated in the developer mode shows the following. The 1n32_basepair.json file stores all the base pair information and the 1n32_metnuc.json file stores the information of all the metal that binds to a base (whether or not the base forms a pair).

#### 1n32_metbp.json
A sample record of metbp.json file will be as follows. Here acc. is the accn number, mode is the mode in which the file has been generated. Then the metal sites as an array. where a site id is given. then the metal detail, then the *bases* tag which gives the details of the bases that forms pair, then the attrib tag gives the base pair type, loc of the bind, atom name, E-value, distance of the atom from the metal and an extra tartiary information to indicate whether the base to which the metal binds has made more than one base pairs or not. If tartiary=true, then it forms more than one base pairs. It can be noted that the siteid may not be unique. If the metal binds with more than one nucleic acids, then for all those sites the siteid will be same. In fact this is a handy tool to check whether a metal has binds to more than one nucleic acids or not.

		{
		  "accn":"1n32",
		  "mode":"bp",
		  "metbp_sites":[
				  {
				  "siteid":"1n32_MG_1563_A",
				  "metal":
				  {
					"name":"MG", 
					"resid":1563, 
					"chain":"A"
				  },
				  "bases": 
				  {
					"name1":"G", 
					"resid1":858, 
					"chain1":"A",
					"ins1":null,
					"name2":"A", 
					"resid2":828, 
					"chain2":"A",
					"ins2":null
				  },
				  "attrib": 
				  {
					"bptype":"S:HT",
					"loc": "NUC",
					"atom":"N7",
					"eval":0.22  ,
					"dist":2.45  
				  },
				  "tartiary":false 
			    }
			    ]
		}
#### 1n32_metnuc.json
The metnuc.json file gives the information about metal-nucleic acid interactions whether or not the nucleic acid forms a pair. This file is generated only in developer mode (-mode=dev). siteid may not be unique as stated abobe.

	{
	  "accn":"1n32",
	  "mode":"dev",
	  "metal_sites":[
			{
				"siteid":"1n32_MG_1546_A",
				"metal_name":"MG", 
				"metal_resid":1546, 
				"metal_chain":"A", 
				"base_name":"G", 
				"base_resid":944, 
				"base_chain":"A",
				"base_ins":null,
				"loc": "PHP",
				"atom":"OP1",
				"dist":2.02  ,
				"paired":true 
			}
			]
	}
#### 1n32_basepair.json.
This file presents the base pair details. Here all the information are self explanatory. The ins code is set to null if the ins code is missing.

	{
	    "accn":"1n32",
	    "basepairs":[
		{
		    "resnum1":9,
		    "chain1":"A",
		    "ins1":null,
		    "resname1":"G",
		    "resnum2":25,
		    "chain2":"A",
		    "ins2":null,
		    "resname2":"C",
		    "basepair":"W:WC",
		    "eval":0.20
		}
		]
	}

## bug-report:

		Parthajit Roy, 
			`roy.parthajit@gmail.com`

		Dhananjay Bhattacharyya,
			`dhananjay.bhattacharyya@saha.ac.in` 
