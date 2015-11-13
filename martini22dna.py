################################
## 6 # FORCE FIELD PARAMETERS ##  -> @FF <-
################################

# New martini 2.2 parameters.
# Changed: 
#   Unstructured Pro backbone bead
#   Proline side chains
#   Phe sidechain
#   Trp sidechain
#   Helix BB-bonds to constraint      
import MAP

class martini22dna:
    def __init__(self):
        import SS,FUNC,IO

        # parameters are defined here for the following (protein) forcefields:
        self.name = 'martini22dna'
        
        # Charged types:
        self.charges = {"Qd":1, "Qa":-1, "SQd":1, "SQa":-1, "RQd":1, "AQa":-1}                                                           #@#
        self.bbcharges = {"BB1":-1}                                                                                                      #@#
        
        
        #----+---------------------+
        ## A | BACKBONE PARAMETERS |
        #----+---------------------+
        #
        # bbss  lists the one letter secondary structure code
        # bbdef lists the corresponding default backbone beads
        # bbtyp lists the corresponding residue specific backbone beads
        #
        # bbd   lists the structure specific backbone bond lengths
        # bbkb  lists the corresponding bond force constants
        #
        # bba   lists the structure specific angles
        # bbka  lists the corresponding angle force constants
        #
        # bbd   lists the structure specific dihedral angles
        # bbkd  lists the corresponding force constants
        #
        # -=NOTE=- 
        #  if the secondary structure types differ between bonded atoms
        #  the bond is assigned the lowest corresponding force constant 
        #
        # -=NOTE=-
        # if proline is anywhere in the helix, the BBB angle changes for 
        # all residues
        #
        
        ###############################################################################################
        ## BEADS ##                                                                         #                 
        #                              F     E     H     1     2     3     T     S     C    # SS one letter   
        self.bbdef    =    FUNC.spl(" N0   Nda    N0    Nd    Na   Nda   Nda    P5    P5")  # Default beads   #@#
        self.bbtyp    = {                                                                   #                 #@#
                    "ALA": FUNC.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4"), # ALA specific    #@#
                    "PRO": FUNC.spl(" C5    N0    C5    N0    Na    N0    N0    P4    P4"), # PRO specific    #@#
                    "HYP": FUNC.spl(" C5    N0    C5    N0    N0    N0    N0    P4    P4")  # HYP specific    #@#
        }                                                                                   #                 #@#
        ## BONDS ##                                                                         #                 
        self.bbldef   =             (.365, .350, .310, .310, .310, .310, .350, .350, .350)  # BB bond lengths #@#
        self.bbkb     =             (1250, 1250, None, None, None, None, 1250, 1250, 1250)  # BB bond kB      #@#
        self.bbltyp   = {}                                                                  #                 #@#
        self.bbkbtyp  = {}                                                                  #                 #@#
        ## ANGLES ##                                                                        #                 
        self.bbadef   =             ( 119.2,134,   96,   96,   96,   96,  100,  130,  127)  # BBB angles      #@#
        self.bbka     =             ( 150,   25,  700,  700,  700,  700,   20,   20,   20)  # BBB angle kB    #@#
        self.bbatyp   = {                                                                   #                 #@#
               "PRO":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127), # PRO specific    #@#
               "HYP":               ( 119.2,134,   98,   98,   98,   98,  100,  130,  127)  # PRO specific    #@#
        }                                                                                   #                 #@#
        self.bbkatyp  = {                                                                   #                 #@#
               "PRO":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25), # PRO specific    #@#
               "HYP":               ( 150,   25,  100,  100,  100,  100,   25,   25,   25)  # PRO specific    #@#
        }                                                                                   #                 #@#
        ## DIHEDRALS ##                                                                     #                 
        self.bbddef   =             ( 90.7,   0, -120, -120, -120, -120)                    # BBBB dihedrals  #@#
        self.bbkd     =             ( 100,   10,  400,  400,  400,  400)                    # BBBB kB         #@#
        self.bbdmul   =             (   1,    1,    1,    1,    1,    1)                    # BBBB mltplcty   #@#
        self.bbdtyp   = {}                                                                  #                 #@#
        self.bbkdtyp  = {}                                                                  #                 #@#
                                                                                            #                 
        ###############################################################################################               
        
        # Some Forcefields use the Ca position to position the BB-bead (me like!)
        # martini 2.1 doesn't
        self.ca2bb = False 
        
        # BBS angle, equal for all ss types                                                         
        # Connects BB(i-1),BB(i),SC(i), except for first residue: BB(i+1),BB(i),SC(i)               
        #                 ANGLE   Ka                                                                
        self.bbsangle =      [   100,  25]                                                               #@#
        
        # Bonds for extended structures (more stable than using dihedrals)                          
        #               LENGTH FORCE                                                                
        self.ebonds   = {                                                                                #@#
               'short': [ .640, 2500],                                                              #@#
               'long' : [ .970, 2500]                                                               #@#
        }                                                                                           #@#
        
        
        #----+-----------------------+
        ## B | SIDE CHAIN PARAMETERS |
        #----+-----------------------+
        
        # To be compatible with Elnedyn, all parameters are explicitly defined, even if they are double.
        self.sidechains = {
            #RES#   BEADS                       BONDS                                                   ANGLES                      DIHEDRALS
            #                                   BB-SC          SC-SC                                    BB-SC-SC  SC-SC-SC
            "TRP": [FUNC.spl("SC4 SNd SC5 SC5"),[(0.300,5000)]+[(0.270,None) for i in range(5)],        [(210,50),(90,50),(90,50)], [(0,50),(0,200)]],
            "TYR": [FUNC.spl("SC4 SC4 SP1"),    [(0.320,5000), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "PHE": [FUNC.spl("SC5 SC5 SC5"),    [(0.310,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIS": [FUNC.spl("SC4 SP1 SP1"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "HIH": [FUNC.spl("SC4 SP1 SQd"),    [(0.320,7500), (0.270,None), (0.270,None),(0.270,None)],[(150,50),(150,50)],        [(0,50)]],
            "ARG": [FUNC.spl("N0 Qd"),          [(0.330,5000), (0.340,5000)],                           [(180,25)]],
            "LYS": [FUNC.spl("C3 Qd"),          [(0.330,5000), (0.280,5000)],                           [(180,25)]],
            "CYS": [FUNC.spl("C5"),             [(0.310,7500)]],
            "ASP": [FUNC.spl("Qa"),             [(0.320,7500)]],
            "GLU": [FUNC.spl("Qa"),             [(0.400,5000)]],
            "ILE": [FUNC.spl("AC1"),            [(0.310,None)]],
            "LEU": [FUNC.spl("AC1"),            [(0.330,7500)]],
            "MET": [FUNC.spl("C5"),             [(0.400,2500)]],
            "ASN": [FUNC.spl("P5"),             [(0.320,5000)]],
            "PRO": [FUNC.spl("C3"),             [(0.300,7500)]],
            "HYP": [FUNC.spl("P1"),             [(0.300,7500)]],
            "GLN": [FUNC.spl("P4"),             [(0.400,5000)]],
            "SER": [FUNC.spl("P1"),             [(0.250,7500)]],
            "THR": [FUNC.spl("P1"),             [(0.260,None)]],
            "VAL": [FUNC.spl("AC2"),            [(0.265,None)]],
            "ALA": [],
            "GLY": [],
            }
        
        # Not all (eg Elnedyn) forcefields use backbone-backbone-sidechain angles and BBBB-dihedrals.
        self.UseBBSAngles          = True 
        self.UseBBBBDihedrals      = True

        # Martini 2.2p has polar and charged residues with seperate charges.
        self.polar   = []
        self.charged = []

        # If masses or charged diverge from standard (45/72 and -/+1) they are defined here.
        self.mass_charge = {
        #RES   MASS               CHARGE
        }

        # Defines the connectivity between between beads
        self.aa_connectivity = {
        #RES       BONDS                                   ANGLES             DIHEDRALS              V-SITE
        "TRP":     [[(0,1),(1,2),(1,3),(2,3),(2,4),(3,4)], [(0,1,2),(0,1,3)], [(0,2,3,1),(1,2,4,3)]],  
        "TYR":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]], 
        "PHE":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIS":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "HIH":     [[(0,1),(1,2),(1,3),(2,3)],             [(0,1,2),(0,1,3)], [(0,2,3,1)]],
        "GLN":     [[(0,1)]],
        "ASN":     [[(0,1)]],
        "SER":     [[(0,1)]],
        "THR":     [[(0,1)]],
        "ARG":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "LYS":     [[(0,1),(1,2)],                         [(0,1,2)]],
        "ASP":     [[(0,1)]],
        "GLU":     [[(0,1)]],
        "CYS":     [[(0,1)]],
        "ILE":     [[(0,1)]],
        "LEU":     [[(0,1)]],
        "MET":     [[(0,1)]],
        "PRO":     [[(0,1)]],
        "HYP":     [[(0,1)]],
        "VAL":     [[(0,1)]],
        "ALA":     [],
        "GLY":     [],
        }
        
        #----+----------------+
        ## C | DNA/RNA bases  |
        #----+----------------+

        # DNA BACKBONE PARAMETERS
        self.dna_bb = {
            'atom'  : FUNC.spl("Q0 SN0 SC2"),
            'bond'  : [(1,  0.360, 20000),          
                       (1,  0.198, 80000),          
                       (1,  0.353, 10000)],         
            'angle' : [(2,  110.0, 200),            
                       (2,  102.0, 150),           
                       (2,  106.0,  75)],           
            'dih'   : [(2,   95.0,  25),
                       (1,  180.0,   2, 3),
                       (9,   85.0,   2, 2,  9,  160.0,  2, 3)],
            'excl'  : [(), (), ()],
            'pair'  : [],
        }
        # DNA BACKBONE CONNECTIVITY
        self.dna_con  = {
            'bond'  : [(0, 1),
                       (1, 2),
                       (2, 0)],
            'angle' : [(0, 1, 2),
                       (1, 2, 0),
                       (2, 0, 1)],
            'dih'   : [(0, 1, 2, 0),
                       (1, 2, 0, 1),
                       (2, 0, 1, 2)],
            'excl'  : [(0, 2), (1, 0),(2, 1)],
            'pair'  : [],
        }

## FOR PLOTTING ONLY
#        # DNA BACKBONE PARAMETERS
#        self.dna_bb = {
#            'atom'  : FUNC.spl("Q0 SN0 SC2"),
#            'bond'  : [(1,  0.360, 30000),          
#                       (1,  0.400, 10000),          
#                       (1,  0.200, 50000),          
#                       (1,  0.355, 10000)],         
#            'angle' : [(2,  115.0,  85),           
#                       (2,  102.0, 105),           
#                       (2,  110.0,  60)],           
#            'dih'   : [(2,  100.0,  1),           
#                       (2, -120.0,  5),           
#                       (2,  140.0,  5)],          
#            'excl'  : [(), (), ()],
#        }
#        # DNA BACKBONE CONNECTIVITY
#        self.dna_con  = {
#            'bond'  : [(0, 1),
#                       (0, 2),
#                       (1, 2),
#                       (2, 0)],
#            'angle' : [(0, 1, 2),
#                       (1, 2, 0),
#                       (2, 0, 1)],
#            'dih'   : [(0, 1, 2, 0),
#                       (1, 2, 0, 1),
#                       (2, 0, 1, 2)],
#            'excl'  : [(0, 2), (1, 0), (2, 1)],
#        }

        # RNA BACKBONE PARAMETERS
        self.rna_bb = {
            'atom'  : FUNC.spl("Q0 N0 C2"),
            'bond'  : [(0.120,5000),(0.220,5000),(0.320,5000)],
            'angle' : [(10.0, 100), (20.0, 100), (30.0, 100)],
            'dih'   : [(100, 10), (100, 10), (100, 10),],
            'excl'  : [],
        }
        # RNA BACKBONE CONNECTIVITY
        self.rna_con  = {
            'bond'  : [(0,1),(1,2),(2,0)],
            'angle' : [(0,1,2),(1,2,0),(2,0,1)],
            'dih'   : [(0,1,2,0),(1,2,0,1),(2,0,1,2)],
            'excl'  : [],
        }

        # For bonds, angles, and dihedrals the first parameter should always 
        # be the type. It is pretty annoying to check the connectivity from 
        # elsewhere so we update these one base at a time.

        # ADENINE
        self.bases = {
            "DA": [FUNC.spl("TN0 TA2 TA3 TNa"),                                      
            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.348, 20000), (1,  0.229,  None), (1,  0.266,  None),             # BONDS BB3-SC1 bond lengthened by 0.048 nm.
                   [(1,  0.300, 30000), (1,  0.229,  None), (1,  0.266,  None),             # BONDS BB3-SC1 bond lengthened by 0.048 nm.
                    (1,  0.326, 20000), (1,  0.288,  None), (1,  0.162,  None),],     
                   [(2,   94.0,   250), (2,  160.0,   200), (2,  140.0,   200),             # ANGLES
                    (1,   85.0,   200), (2,  158.0,   200), (1,  125.0,   200),
                    (1,   74.0,   200), (1,   98.0,   200)],                           
                   [(2,  -90.0,    20), (2, -116.0,   0.5), (2,   98.0,    15)],            # DIHEDRALS
                   [],                                                                      # IMPROPERS 
                   [],                                                                      # VSITES
                   [(), (), (), (), (), (), (), (), (), (), (), (), (), ()],                # EXCLUSIONS
                   []],                                                                     # PAIRS
            }
        self.base_connectivity = {
            "DA": [[(2, 3),             (3, 4),             (4, 5),                         # BONDS
                    (4, 6),             (5, 6),             (6, 3)],   
                   [(1, 2, 3),          (2, 3, 4),          (2, 3, 6),                      # ANGLES
                    (3, 4, 5),          (3, 2, 7),          (4, 3, 6),
                    (4, 5, 6),          (5, 6, 3)], 
                   [(0, 1, 2, 3),       (1, 2, 3, 4),       (1, 2, 3, 6),],                # DIHEDRALS        
                   [],                                                                      # IMPROPERS
                   [],                                                                      # VSITES
                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
                    (0, 6),             (1, 3),             (1, 4),
                    (1, 5),             (1, 6),             (2, 3),
                    (2, 4),             (2, 5),             (2, 6),
                    (3, 5),             (4, 6)],
                   []],                                                                     # PAIRS                     
            }

        # CYTOSINE
        self.bases.update({
            "DC": [FUNC.spl("TN0 TY2 TY3"),                                                     
            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.303, 20000), (1,  0.220,  None), (1,  0.285,  None),             # BONDS BB3-SC1 bond lenghtened by 0.033 nm.
                   [(1,  0.270, 30000), (1,  0.220,  None), (1,  0.285,  None),             # BONDS BB3-SC1 bond lenghtened by 0.033 nm.
                    (1,  0.268,  None),],
                   [(2,   95.0,   210), (2,   95.0,   300), (1,  150.0,   500),             # ANGLES
                    (1,  180.0,    30), (1,   61.0,   200), (1,   71.0,   200), 
                    (1,   47.0,   200)],
                   [(2,  -78.0,    25), (2,  -90.0,    20), (2, -142.0,    50)],            # DIHEDRALS
                   #[(2,  -78.0,    25), (2, -108.0,    10), (2,   40.0,    15)],            # DIHEDRALS
                   [],                                                                      # IMPROPERS
                   [],                                                                      # VSITES
                   [(), (), (), (), (), (), (), (), ()],                                    # EXCLUSIONS
                   []],                                                                     # PAIRS                     
        })
        self.base_connectivity.update({
            "DC": [[(2, 3),           (3, 4),             (4, 5),                         # BONDS
                    (5, 3)],
                   [(1, 2, 3),        (2, 3, 4),          (1, 3, 5),                      # ANGLES
                    (3, 2, 6),        (3, 4, 5),          (4, 3, 5),
                    (4, 5, 3)],
                   [(0, 1, 2, 3),     (1, 2, 3, 4),       (2, 1, 3, 5)],                  # DIHEDRALS
                   [],                                                                    # IMPROPERS
                   [],                                                                    # VSITES
                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
                    (1, 3),             (1, 4),             (1, 5),             
                    (2, 3),             (2, 4),             (2, 5)],                                           
                   []],                                                                     # PAIRS                     
        })

        # GUANINE
        self.bases.update({
            "DG": [FUNC.spl("TN0 TG2 TG3 TNa"),
            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.353, 20000), (1,  0.295,  None), (1,  0.295,  None),             # BONDS BB3-SC1 bond lengthened by 0.053 nm.
                   [(1,  0.300, 30000), (1,  0.295,  None), (1,  0.295,  None),             # BONDS BB3-SC1 bond lengthened by 0.053 nm.
                    (1,  0.389, 20000), (1,  0.285,  None), (1,  0.161,  None),],     
                   [(2,   94.5,   250), (2,  137.0,   300), (2,  130.0,   250),             # ANGLES
                    (1,   69.5,   200), (2,  157.0,   150), (1,  125.0,   200),
                    (1,   84.0,   200), (1,   94.0,   200)],                           
                   [(2,  -90.0,    20), (2, -117.0,     1), (2,   92.0,    15)],            # DIHEDRALS  
                   [],                                                                      # IMPROPERS 
                   [],                                                                      # VSITES
                   [(), (), (), (), (), (), (), (), (), (), (), (), (), ()],                # EXCLUSIONS
                   []],                                                                     # PAIRS                     
        })
        self.base_connectivity.update({
            "DG": [[(2, 3),             (3, 4),             (4, 5),                         # BONDS
                    (4, 6),             (5, 6),             (6, 3)],
                   [(1, 2, 3),          (2, 3, 4),          (2, 3, 6),                      # ANGLES
                    (3, 4, 5),          (3, 2, 7),          (4, 3, 6), 
                    (4, 5, 6),          (5, 6, 3)],
                   [(0, 1, 2, 3),       (1, 2, 3, 4),       (1, 2, 3, 6),],                 # DIHEDRALS        
                   [],                                                                      # IMPROPERS
                   [],                                                                      # VSITES
                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
                    (0, 6),             (1, 3),             (1, 4),
                    (1, 5),             (1, 6),             (2, 3),
                    (2, 4),             (2, 5),             (2, 6),
                    (3, 5),             (4, 6)],                                           
                   []],                                                                     # PAIRS                     
        })

        # THYMINE
        self.bases.update({
            "DT": [FUNC.spl("TN0 TT2 TT3"),
            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.326, 20000), (1,  0.217,  None), (1,  0.322,  None),             # BONDS BB3-SC1 bond lengthened by 0.056 nm.
                   [(1,  0.270, 30000), (1,  0.217,  None), (1,  0.322,  None),             # BONDS BB3-SC1 bond lengthened by 0.056 nm.
                    (1,  0.265,  None),],
                   [(2,   92.0,   220), (2,  107.0,   300), (1,  145.0,   400),             # ANGLES
                    (1,  180.0,    30), (1,   55.0,   100), (1,   83.0,   100),
                    (1,   42.0,   100)],
                   [(2,  -75.0,    40), (2, -110.0,    15), (2, -145.0,    65)],            # DIHEDRALS
                   [],                                                                      # IMPROPERS
                   [],                                                                      # VSITES
                   [(), (), (), (), (), (), (), (), ()],                                    # EXCLUSIONS
                   []],                                                                     # PAIRS                     
        })
        self.base_connectivity.update({
            "DT": [[(2, 3),           (3, 4),             (4, 5),                         # BONDS
                    (5, 3)],
                   [(1, 2, 3),        (2, 3, 4),          (1, 3, 5),                      # ANGLES
                    (3, 2, 6),        (3, 4, 5),          (4, 3, 5), 
                    (4, 5, 3)],
                   [(0, 1, 2, 3),     (1, 2, 3, 4),       (2, 1, 3, 5)],                  # DIHEDRALS
                   [],                                                                    # IMPROPERS
                   [],                                                                    # VSITES
                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
                    (1, 3),             (1, 4),             (1, 5),             
                    (2, 3),             (2, 4),             (2, 5)],                                           
                   []],                                                                     # PAIRS                     
        })

## FOR PLOTTING ONLY
#        # ADENINE
#        self.bases = {
#            "DA": [FUNC.spl("TN0 TA2 TA3 TNa"),                                      
#            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.330, 30000), (1,  0.229, 30000), (1,  0.266, 30000),             # BONDS BB3-SC1 bond lengthened by 0.030 nm.
#                    (1,  0.325, 30000), (1,  0.288, 30000), (1,  0.162, 30000),],     
#                   [(2,   93.0,   250), (2,  160.0,   200), (2,  140.0,   200),             # ANGLES
#                    (2,   85.0,   200), (2,  148.0,   350), (2,  125.0,   200),
#                    (2,   74.0,   200), (2,   98.0,   200)],                           
#                   [(2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1)],
#                   [(2,    0.0,   500)],                                                    # IMPROPERS 
#                   [],                                                                      # VSITES
#                   [(), (), (), (), (), (), (), (), (), (), (), (), (), ()]],               # EXCLUSIONS
#            }
#        self.base_connectivity = {
#            "DA": [[(2, 3),             (3, 4),             (4, 5),                         # BONDS
#                    (4, 6),             (5, 6),             (6, 3)],   
#                   [(1, 2, 3),          (2, 3, 4),          (2, 3, 6),                      # ANGLES
#                    (3, 4, 5),          (3, 2, 7),          (4, 3, 6),
#                    (4, 5, 6),          (5, 6, 3)], 
#                   [(0, 1, 2, 3),       (0, 2, 3, 4),       (0, 2, 3, 6),                  # DIHEDRALS        
#                    (1, 2, 3, 4),       (1, 2, 3, 6),       (2, 8, 9,10),
#                    (3, 2, 7, 8),       (3, 2, 7, 9),       (3, 2, 7,10),
#                    (3, 7, 8,10),       (3, 8, 9,10),       (4, 2, 7, 8),
#                    (7, 2, 3, 4),       (7, 2, 3, 6)],
#                   [(3, 4, 5, 6)],                                                          # IMPROPERS
#                   [],                                                                      # VSITES
#                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
#                    (0, 6),             (1, 3),             (1, 4),
#                    (1, 5),             (1, 6),             (2, 3),
#                    (2, 4),             (2, 5),             (2, 6),
#                    (3, 5),             (4, 6)]],                                           
#            }
#
#        # CYTOSINE
#        self.bases.update({
#            "DC": [FUNC.spl("TN0 TY2 TY3"),                                                     
#            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.290, 30000), (1,  0.220, 30000), (1,  0.285, 30000),             # BONDS BB3-SC1 bond lenghtened by 0.020 nm.
#                    (1,  0.268, 30000),],
#                   [(2,   93.0,   200), (2,  108.0,   250), (2,  170.0,   350),             # ANGLES
#                    (2,  180.0,     1), (2,   62.0,   200), (2,   71.0,   200), 
#                    (2,   47.0,   200), (2,  100.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1)],
#                   [(2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1)],
#                   [],                                                                      # IMPROPERS
#                   [],                                                                      # VSITES
#                   [(), (), (), (), (), (), (), (), ()]],                                   # EXCLUSIONS
#        })
#        self.base_connectivity.update({
#            "DC": [[(2, 3),           (3, 4),             (4, 5),                         # BONDS
#                    (5, 3)],
#                   [(1, 2, 3),        (2, 3, 4),          (2, 3, 5),                      # ANGLES
#                    (3, 2, 6),        (3, 4, 5),          (4, 3, 5),
#                    (4, 5, 3),        (1, 3, 5),          (1, 5, 3),
#                    (2, 3, 6),        (2, 1, 3),          (2, 1, 5)],
#                   [(0, 1, 2, 3),       (0, 2, 3, 4),       (0, 2, 3, 5),                  # DIHEDRALS        
#                    (1, 2, 3, 4),       (1, 2, 3, 5),       (2, 7, 8, 9),
#                    (3, 2, 6, 7),       (3, 2, 6, 8),       (3, 2, 6, 9),
#                    (3, 6, 7, 9),       (3, 7, 8, 9),       (4, 2, 6, 7),
#                    (6, 2, 3, 4),       (6, 2, 3, 5),       (2, 1, 3, 5),
#                    (2, 1, 5, 3)],
#                   [],                                                                    # IMPROPERS
#                   [],                                                                    # VSITES
#                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
#                    (1, 3),             (1, 4),             (1, 5),             
#                    (2, 3),             (2, 4),             (2, 5)]],                                           
#        })
#
#        # GUANINE
#        self.bases.update({
#            "DG": [FUNC.spl("TN0 TG2 TG3 TNa"),
#            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.300, 30000), (1,  0.295, 30000), (1,  0.295, 30000),             # BONDS BB3-SC1 bond stays the same.
#                    (1,  0.390, 30000), (1,  0.285, 30000), (1,  0.161, 30000),],     
#                   [(2,   95.0,   250), (2,  137.0,   300), (2,  128.0,   250),             # ANGLES
#                    (2,   69.0,   200), (2,  145.0,   350), (2,  125.0,   200),
#                    (2,   84.0,   200), (2,   94.0,   200)],                           
#                   [(2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1)],
#                   [(2,    0.0,   150)],                                                    # IMPROPERS 
#                   [],                                                                      # VSITES
#                   [(), (), (), (), (), (), (), (), (), (), (), (), (), ()]],               # EXCLUSIONS
#        })
#        self.base_connectivity.update({
#            "DG": [[(2, 3),             (3, 4),             (4, 5),                         # BONDS
#                    (4, 6),             (5, 6),             (6, 3)],
#                   [(1, 2, 3),          (2, 3, 4),          (2, 3, 6),                      # ANGLES
#                    (3, 4, 5),          (3, 2, 7),          (4, 3, 6), 
#                    (4, 5, 6),          (5, 6, 3)],
#                   [(0, 1, 2, 3),       (0, 2, 3, 4),       (0, 2, 3, 6),                  # DIHEDRALS        
#                    (1, 2, 3, 4),       (1, 2, 3, 6),       (2, 8, 9,10),
#                    (3, 2, 7, 8),       (3, 2, 7, 9),       (3, 2, 7,10),
#                    (3, 7, 8,10),       (3, 8, 9,10),       (4, 2, 7, 8),
#                    (7, 2, 3, 4),       (7, 2, 3, 6)],
#                   [(3, 4, 5, 6)],                                                          # IMPROPERS
#                   [],                                                                      # VSITES
#                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
#                    (0, 6),             (1, 3),             (1, 4),
#                    (1, 5),             (1, 6),             (2, 3),
#                    (2, 4),             (2, 5),             (2, 6),
#                    (3, 5),             (4, 6)]],                                           
#        })
#
#        # THYMINE
#        self.bases.update({
#            "DT": [FUNC.spl("TN0 TT2 TT3"),
#            #     TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS TYPE   EQUIL   OPTS
#                   [(1,  0.310, 30000), (1,  0.217, 30000), (1,  0.322, 30000),             # BONDS BB3-SC1 bond lengthened by 0.040 nm.
#                    (1,  0.265, 30000),],
#                   [(2,   93.0,   250), (2,  108.0,   350), (2,  165.0,   550),             # ANGLES
#                    (2,  165.0,   400), (2,   55.0,   200), (2,   83.0,   200),
#                    (2,   42.0,   200), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1)],
#                   [(2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1),
#                    (2,    0.0,     1), (2,    0.0,     1), (2,    0.0,     1)],
#                   [],                                                                      # IMPROPERS
#                   [],                                                                      # VSITES
#                   [(), (), (), (), (), (), (), (), ()]],                                   # EXCLUSIONS
#        })
#        self.base_connectivity.update({
#            "DT": [[(2, 3),           (3, 4),             (4, 5),                         # BONDS
#                    (5, 3)],
#                   [(1, 2, 3),        (2, 3, 4),          (2, 3, 5),                      # ANGLES
#                    (3, 2, 6),        (3, 4, 5),          (4, 3, 5), 
#                    (4, 5, 3),        (2, 3, 6),          (1, 3, 5),
#                    (2, 1, 3),        (2, 1, 5)],
#                   [(0, 1, 2, 3),       (0, 2, 3, 4),       (0, 2, 3, 5),                  # DIHEDRALS        
#                    (1, 2, 3, 4),       (1, 2, 3, 5),       (2, 7, 8, 9),
#                    (3, 2, 6, 7),       (3, 2, 6, 8),       (3, 2, 6, 9),
#                    (3, 6, 7, 9),       (3, 7, 8, 9),       (4, 2, 6, 7),
#                    (6, 2, 3, 4),       (6, 2, 3, 5),       (2, 1, 3, 5)],
#                   [],                                                                    # IMPROPERS
#                   [],                                                                    # VSITES
#                   [(0, 3),             (0, 4),             (0, 5),                         # EXCLUSIONS
#                    (1, 3),             (1, 4),             (1, 5),             
#                    (2, 3),             (2, 4),             (2, 5)]],                                           
#        })


        #----+----------------+
        ## D | SPECIAL BONDS  |
        #----+----------------+
        
        self.special = {
            # Used for sulfur bridges
            # ATOM 1         ATOM 2          BOND LENGTH   FORCE CONSTANT
            (("SC1","CYS"), ("SC1","CYS")):     (0.39,         5000),
            }
        
        # By default use an elastic network
        self.ElasticNetwork = False 

        # Elastic networks bond shouldn't lead to exclusions (type 6) 
        # But Elnedyn has been parametrized with type 1.
        self.EBondType = 6
        
        #----+----------------+
        ## D | INTERNAL STUFF |
        #----+----------------+
        
        
        ## BACKBONE BEAD TYPE ##                                                                    
        # Dictionary of default bead types (*D)                                                     
        self.bbBeadDictD  = FUNC.hash(SS.bbss,self.bbdef)                                                             
        # Dictionary of dictionaries of types for specific residues (*S)                            
        self.bbBeadDictS  = dict([(i,FUNC.hash(SS.bbss,self.bbtyp[i])) for i in self.bbtyp.keys()])                        

        # combine the connectivity records for different molecule types
        self.connectivity = dict(self.base_connectivity.items() + self.aa_connectivity.items())
        # XXX No need to do that, let's just use separate for DNA for now
        
        ## BB BOND TYPE ##                                                                          
        # Dictionary of default abond types (*D)                                                    
        self.bbBondDictD = FUNC.hash(SS.bbss,zip(self.bbldef,self.bbkb))                                                   
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbBondDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbltyp[i],self.bbkbtyp[i]))) for i in self.bbltyp.keys()])       
        # This is tricky to read, but it gives the right bondlength/force constant
        
        ## BBB ANGLE TYPE ##                                                                        
        # Dictionary of default angle types (*D)                                                    
        self.bbAngleDictD = FUNC.hash(SS.bbss,zip(self.bbadef,self.bbka))                                                  
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbAngleDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbatyp[i],self.bbkatyp[i]))) for i in self.bbatyp.keys()])      
                    
        ## BBBB DIHEDRAL TYPE ##                                                                    
        # Dictionary of default dihedral types (*D)                                                 
        self.bbDihedDictD = FUNC.hash(SS.bbss,zip(self.bbddef,self.bbkd,self.bbdmul))                                           
        # Dictionary of dictionaries for specific types (*S)                                        
        self.bbDihedDictS = dict([(i,FUNC.hash(SS.bbss,zip(self.bbdtyp[i],self.bbkdtyp[i]))) for i in self.bbdtyp.keys()])      

        ## DNA DICTIONARIES ##
        # Dictionary for the connectivities and parameters of bonds between DNA backbone beads
        self.dnaBbBondDictC = dict(zip(self.dna_con['bond'],self.dna_bb['bond']))
        # Dictionary for the connectivities and parameters of angles between DNA backbone beads
        self.dnaBbAngleDictC = dict(zip(self.dna_con['angle'],self.dna_bb['angle']))
        # Dictionary for the connectivities and parameters of dihedrals between DNA backbone beads
        self.dnaBbDihDictC = dict(zip(self.dna_con['dih'],self.dna_bb['dih']))
        # Dictionary for exclusions for DNA backbone beads
        self.dnaBbExclDictC = dict(zip(self.dna_con['excl'],self.dna_bb['excl']))
        # Dictionary for pairs for DNA backbone beads
        self.dnaBbPairDictC = dict(zip(self.dna_con['pair'],self.dna_bb['pair']))

        ## RNA DICTIONARIES ##
        # Dictionary for the connectivities and parameters of bonds between RNA backbone beads
        self.rnaBbBondDictC = dict(zip(self.rna_con['bond'],self.rna_bb['bond']))
        # Dictionary for the connectivities and parameters of angles between rna backbone beads
        self.rnaBbAngleDictC = dict(zip(self.rna_con['angle'],self.rna_bb['angle']))
        # Dictionary for the connectivities and parameters of dihedrals between rna backbone beads
        self.rnaBbDihDictC = dict(zip(self.rna_con['dih'],self.rna_bb['dih']))
        # Dictionary for exclusions for RNA backbone beads
        self.rnaBbExclDictC = dict(zip(self.rna_con['excl'],self.rna_bb['excl']))
        
        
    # The following function returns the backbone bead for a given residue and                   
    # secondary structure type.                                                                 
    # 1. Check if the residue is DNA/RNA and return the whole backbone for those
    # 2. Look up the proper dictionary for the residue                                          
    # 3. Get the proper type from it for the secondary structure                                
    # If the residue is not in the dictionary of specials, use the default                      
    # If the secondary structure is not listed (in the residue specific                         
    # dictionary) revert to the default.                                                        
    def bbGetBead(self,r1,ss="C"):                                                               
        if r1 in MAP.dnares3:
            return self.dna_bb['atom']
        elif r1 in MAP.rnares3:
            return self.rna_bb['atom']
        else:
            return self.bbBeadDictS.get(r1,self.bbBeadDictD).get(ss,self.bbBeadDictD.get(ss))                      
    
    def bbGetBond(self,r,ca,ss):
        # Retrieve parameters for each residue from tables defined above
        # Check is it DNA residue
        if r[0] in MAP.dnares3:
            return ca in self.dnaBbBondDictC.keys() and self.dnaBbBondDictC[ca] or None
        # RNA is not implemented properly yet
        elif r[0] in MAP.rnares3:
            return ca in self.rnaBbBondDictC.keys() and self.rnaBbBondDictC[ca] or None
        # If it's protein
        else:
            b1 = self.bbBondDictS.get(r[0],self.bbBondDictD).get(ss[0],self.bbBondDictD.get(ss[0]))
            b2 = self.bbBondDictS.get(r[1],self.bbBondDictD).get(ss[1],self.bbBondDictD.get(ss[1]))
            # Determine which parameters to use for the bond
            return ( (b1[0]+b2[0])/2, min(b1[1],b2[1]) )
    
    def bbGetAngle(self,r,ca,ss):
        # Check is it DNA residue
        if r[0] in MAP.dnares3:
            return ca in self.dnaBbAngleDictC.keys() and self.dnaBbAngleDictC[ca] or None
        # RNA is not implemented properly yet
        elif r[0] in MAP.rnares3:
            return ca in self.rnaBbAngleDictC.keys() and self.rnaBbAngleDictC[ca] or None
        # For protein
        else:
            # PRO in helices is dominant
            if r[1] == "PRO" and ss[1] in "H123":
                return self.bbAngleDictS["PRO"].get(ss[1])
            else:
                # Retrieve parameters for each residue from table defined above
                a = [ self.bbAngleDictS.get(r[0],self.bbAngleDictD).get(ss[0],self.bbAngleDictD.get(ss[0])),
                      self.bbAngleDictS.get(r[1],self.bbAngleDictD).get(ss[1],self.bbAngleDictD.get(ss[1])),
                      self.bbAngleDictS.get(r[2],self.bbAngleDictD).get(ss[2],self.bbAngleDictD.get(ss[2])) ]
                # Sort according to force constant
                a.sort(key=lambda i: (i[1],i[0]))
                # This selects the set with the smallest force constant and the smallest angle
                return a[0]

    def bbGetExclusion(self,r,ca,ss):
        if r[0] in MAP.dnares3:
            return ca in self.dnaBbExclDictC.keys() and ' ' or None
        # RNA is not implemented properly yet
        elif r[0] in MAP.rnares3:
            return ca in self.rnaBbExclDictC.keys() and ' ' or None
        else:
            return None

    def bbGetPair(self,r,ca,ss):
        if r[0] in MAP.dnares3:
            return ca in self.dnaBbPairDictC.keys() and ' ' or None
        # RNA is not implemented properly yet
        elif r[0] in MAP.rnares3:
            return ca in self.rnaBbPairDictC.keys() and ' ' or None
        else:
            return None

    def bbGetDihedral(self,r,ca,ss):
        # Retrieve parameters for each residue from table defined above
        # Check is it DNA residue
        if r[0] in MAP.dnares3:
            return ca in self.dnaBbDihDictC.keys() and self.dnaBbDihDictC[ca] or None
        # RNA is not implemented properly yet
        elif r[0] in MAP.rnares3:
            return ca in self.rnaBbDihDictC.keys() and self.rnaBbDihDictC[ca] or None
        # Apparently protein has none currently

    def getCharge(self,atype,aname):
        return self.charges.get(atype,self.bbcharges.get(aname,0))
        
    def messages(self):
        '''Prints any force-field specific logging messages.'''
        import logging
        logging.warning('#####################################################################################')
        logging.warning('This is a version of martinize for DNA and should NOT be used for proteins.')
        logging.warning('#####################################################################################')
        pass
