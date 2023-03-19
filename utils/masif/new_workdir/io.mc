##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Mon Sep  5 14:49:11 2022

##############################################################################
Hello world from PE 0
Vnm_tstart: starting timer 26 (APBS WALL CLOCK)..
NOsh_parseInput:  Starting file parsing...
NOsh: Parsing READ section
NOsh: Storing molecule 0 path temp1
NOsh: Done parsing READ section
NOsh: Done parsing READ section (nmol=1, ndiel=0, nkappa=0, ncharge=0, npot=0)
NOsh: Parsing ELEC section
NOsh_parseMG: Parsing parameters for MG calculation
NOsh_parseMG:  Parsing dime...
PBEparm_parseToken:  trying dime...
MGparm_parseToken:  trying dime...
NOsh_parseMG:  Parsing cglen...
PBEparm_parseToken:  trying cglen...
MGparm_parseToken:  trying cglen...
NOsh_parseMG:  Parsing fglen...
PBEparm_parseToken:  trying fglen...
MGparm_parseToken:  trying fglen...
NOsh_parseMG:  Parsing cgcent...
PBEparm_parseToken:  trying cgcent...
MGparm_parseToken:  trying cgcent...
NOsh_parseMG:  Parsing fgcent...
PBEparm_parseToken:  trying fgcent...
MGparm_parseToken:  trying fgcent...
NOsh_parseMG:  Parsing mol...
PBEparm_parseToken:  trying mol...
NOsh_parseMG:  Parsing lpbe...
PBEparm_parseToken:  trying lpbe...
NOsh: parsed lpbe
NOsh_parseMG:  Parsing bcfl...
PBEparm_parseToken:  trying bcfl...
NOsh_parseMG:  Parsing pdie...
PBEparm_parseToken:  trying pdie...
NOsh_parseMG:  Parsing sdie...
PBEparm_parseToken:  trying sdie...
NOsh_parseMG:  Parsing srfm...
PBEparm_parseToken:  trying srfm...
NOsh_parseMG:  Parsing chgm...
PBEparm_parseToken:  trying chgm...
MGparm_parseToken:  trying chgm...
NOsh_parseMG:  Parsing sdens...
PBEparm_parseToken:  trying sdens...
NOsh_parseMG:  Parsing srad...
PBEparm_parseToken:  trying srad...
NOsh_parseMG:  Parsing swin...
PBEparm_parseToken:  trying swin...
NOsh_parseMG:  Parsing temp...
PBEparm_parseToken:  trying temp...
NOsh_parseMG:  Parsing calcenergy...
PBEparm_parseToken:  trying calcenergy...
NOsh_parseMG:  Parsing calcforce...
PBEparm_parseToken:  trying calcforce...
NOsh_parseMG:  Parsing write...
PBEparm_parseToken:  trying write...
NOsh_parseMG:  Parsing end...
MGparm_check:  checking MGparm object of type 1.
NOsh:  nlev = 6, dime = (129, 129, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 2003 atoms
Valist_getStatistics:  Max atom coordinate:  (30.1, 49.261, 73.206)
Valist_getStatistics:  Min atom coordinate:  (-5.036, 6.725, 35.44)
Valist_getStatistics:  Molecule center:  (12.532, 27.993, 54.323)
NOsh_setupCalcMGAUTO(/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1855):  coarse grid center = 12.532 27.993 54.323
NOsh_setupCalcMGAUTO(/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1860):  fine grid center = 12.532 27.993 54.323
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.497184, 0.593087, 0.533282
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.448711, 0.505125, 0.469945
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.902506, 0.851687, 0.881232 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1972):  coarse mesh center = 12.532 27.993 54.323
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 44.3517 65.9506 88.453
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -19.2877 -9.9646 20.193
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 41.2495 60.321 84.3995
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -16.1855 -4.335 24.2465
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 41.2495 60.321 84.3995
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -16.1855 -4.335 24.2465
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 26.1387
Vpbe_ctor2:  solute dimensions = 37.435 x 44.656 x 40.153
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 74 x 75 x 75 table
Vclist_ctor2:  Using 74 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (46.212, 53.612, 48.842)
Vclist_setupGrid:  Grid lower corner = (-10.574, 1.187, 29.902)
Vclist_assignAtoms:  Have 2810337 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.295187
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 8.539470e-01
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.931100e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 4.266300e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.335036e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.247130e-01
Vprtstp: contraction number = 1.247130e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.636548e-02
Vprtstp: contraction number = 1.312252e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.375602e-03
Vprtstp: contraction number = 1.451593e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 4.440110e-04
Vprtstp: contraction number = 1.869046e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.079681e-04
Vprtstp: contraction number = 2.431653e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.204403e-05
Vprtstp: contraction number = 2.967917e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.055139e-05
Vprtstp: contraction number = 3.292779e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 3.883235e-06
Vprtstp: contraction number = 3.680306e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 1.332091e-06
Vprtstp: contraction number = 3.430364e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 5.111676e-07
Vprtstp: contraction number = 3.837332e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.526239e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.573565e+01
Vpmg_setPart:  lower corner = (-19.2877, -9.9646, 20.193)
Vpmg_setPart:  upper corner = (44.3517, 65.9506, 88.453)
Vpmg_setPart:  actual minima = (-19.2877, -9.9646, 20.193)
Vpmg_setPart:  actual maxima = (44.3517, 65.9506, 88.453)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 6.448376371190E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 5.869000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 2.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 26.1387
Vpbe_ctor2:  solute dimensions = 37.435 x 44.656 x 40.153
Vpbe_ctor2:  solute charge = 4
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 74 x 75 x 75 table
Vclist_ctor2:  Using 74 x 75 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (46.212, 53.612, 48.842)
Vclist_setupGrid:  Grid lower corner = (-10.574, 1.187, 29.902)
Vclist_assignAtoms:  Have 2810337 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -16.1855, -4.335, 24.2465
VPMG::focusFillBound -- New mesh maxs = 41.2495, 60.321, 84.3995
VPMG::focusFillBound -- Old mesh mins = -19.2877, -9.9646, 20.193
VPMG::focusFillBound -- Old mesh maxs = 44.3517, 65.9506, 88.453
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-16.1855, -4.335, 24.2465)
Vpmg_setPart:  upper corner = (41.2495, 60.321, 84.3995)
Vpmg_setPart:  actual minima = (-19.2877, -9.9646, 20.193)
Vpmg_setPart:  actual maxima = (44.3517, 65.9506, 88.453)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-16.1855, -4.335, 24.2465)
VPMG::extEnergy    Disj part upper corner = (41.2495, 60.321, 84.3995)
VPMG::extEnergy    Old lower corner = (-19.2877, -9.9646, 20.193)
VPMG::extEnergy    Old upper corner = (44.3517, 65.9506, 88.453)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.320314 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.302164
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.089505e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (129, 129, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 2.936600e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (065, 065, 065)
Vbuildops: Galer: (033, 033, 033)
Vbuildops: Galer: (017, 017, 017)
Vbuildops: Galer: (009, 009, 009)
Vbuildops: Galer: (005, 005, 005)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 3.889870e-01
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 1.837441e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.347496e-01
Vprtstp: contraction number = 1.347496e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.744615e-02
Vprtstp: contraction number = 1.294709e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.414742e-03
Vprtstp: contraction number = 1.384111e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 4.288411e-04
Vprtstp: contraction number = 1.775930e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.044530e-04
Vprtstp: contraction number = 2.435704e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 3.332397e-05
Vprtstp: contraction number = 3.190331e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.149209e-05
Vprtstp: contraction number = 3.448596e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 4.470234e-06
Vprtstp: contraction number = 3.889836e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 1.601515e-06
Vprtstp: contraction number = 3.582620e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 6.526784e-07
Vprtstp: contraction number = 4.075380e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 1.546854e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 1.598190e+01
Vpmg_setPart:  lower corner = (-16.1855, -4.335, 24.2465)
Vpmg_setPart:  upper corner = (41.2495, 60.321, 84.3995)
Vpmg_setPart:  actual minima = (-16.1855, -4.335, 24.2465)
Vpmg_setPart:  actual maxima = (41.2495, 60.321, 84.3995)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 8.067264099714E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.155000e-03
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 3.437922e+01
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Mon Sep  5 14:50:59 2022

##############################################################################
Vgrid_readDX:  Grid dimensions 129 x 129 x 129 grid
Vgrid_readDX:  Grid origin = (-16.1855, -4.335, 24.2465)
Vgrid_readDX:  Grid spacings = (0.448711, 0.505125, 0.469945)
Vgrid_readDX:  allocating 129 x 129 x 129 doubles for storage