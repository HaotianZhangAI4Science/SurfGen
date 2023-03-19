##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Sat Dec 10 14:53:57 2022

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
NOsh:  nlev = 4, dime = (97, 97, 129)
NOsh: Done parsing ELEC section (nelec = 1)
NOsh: Parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing PRINT section
NOsh: Done parsing file (got QUIT)
Valist_readPQR: Counted 1716 atoms
Valist_getStatistics:  Max atom coordinate:  (28.41, 44.706, 73.206)
Valist_getStatistics:  Min atom coordinate:  (-5.036, 10.782, 35.44)
Valist_getStatistics:  Molecule center:  (11.687, 27.744, 54.323)
NOsh_setupCalcMGAUTO(/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1855):  coarse grid center = 11.687 27.744 54.323
NOsh_setupCalcMGAUTO(/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1860):  fine grid center = 11.687 27.744 54.323
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1872):  Coarse grid spacing = 0.635765, 0.641343, 0.533282
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1874):  Fine grid spacing = 0.582313, 0.585594, 0.469945
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1876):  Displacement between fine and coarse grids = 0, 0, 0
NOsh:  2 levels of focusing with 0.915925, 0.913075, 0.881232 reductions
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1970):  starting mesh repositioning.
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1972):  coarse mesh center = 11.687 27.744 54.323
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1977):  coarse mesh upper corner = 42.2037 58.5284 88.453
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1982):  coarse mesh lower corner = -18.8297 -3.04045 20.193
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1987):  initial fine mesh upper corner = 39.638 55.8525 84.3995
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 1992):  initial fine mesh lower corner = -16.264 -0.3645 24.2465
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2053):  final fine mesh upper corner = 39.638 55.8525 84.3995
NOsh_setupCalcMGAUTO (/home/ubuntu/git/apbs-pdb2pqr/apbs/src/generic/nosh.c, 2058):  final fine mesh lower corner = -16.264 -0.3645 24.2465
NOsh_setupMGAUTO:  Resetting boundary flags
NOsh_setupCalc:  Mapping ELEC statement 0 (1) to calculation 1 (2)
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 23.0623
Vpbe_ctor2:  solute dimensions = 35.902 x 36.217 x 40.153
Vpbe_ctor2:  solute charge = 3
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 72 x 75 table
Vclist_ctor2:  Using 71 x 72 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.522, 45, 48.842)
Vclist_setupGrid:  Grid lower corner = (-10.574, 5.244, 29.902)
Vclist_assignAtoms:  Have 2681616 atom entries
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
Vacc_SASA: Time elapsed: 0.157646
Vpmg_fillco:  done filling coefficient arrays
Vpmg_fillco:  filling boundary arrays
Vpmg_fillco:  done filling boundary arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.156752e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.739700e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 065)
Vbuildops: Galer: (025, 025, 033)
Vbuildops: Galer: (013, 013, 017)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 1.997004e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 4.596003e+00
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.208960e-01
Vprtstp: contraction number = 1.208960e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.624377e-02
Vprtstp: contraction number = 1.343615e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.501024e-03
Vprtstp: contraction number = 1.539682e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 4.771339e-04
Vprtstp: contraction number = 1.907755e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.318591e-04
Vprtstp: contraction number = 2.763565e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.309219e-05
Vprtstp: contraction number = 3.268049e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.496593e-05
Vprtstp: contraction number = 3.473002e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 5.793251e-06
Vprtstp: contraction number = 3.870961e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 2.053823e-06
Vprtstp: contraction number = 3.545200e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 8.137709e-07
Vprtstp: contraction number = 3.962225e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 7.277871e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 7.621269e+01
Vpmg_setPart:  lower corner = (-18.8297, -3.04045, 20.193)
Vpmg_setPart:  upper corner = (42.2037, 58.5284, 88.453)
Vpmg_setPart:  actual minima = (-18.8297, -3.04045, 20.193)
Vpmg_setPart:  actual maxima = (42.2037, 58.5284, 88.453)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 4.512576520042E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 2.426100e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 1.000000e-06
Vnm_tstart: starting timer 27 (Setup timer)..
Setting up PBE object...
Vpbe_ctor2:  solute radius = 23.0623
Vpbe_ctor2:  solute dimensions = 35.902 x 36.217 x 40.153
Vpbe_ctor2:  solute charge = 3
Vpbe_ctor2:  bulk ionic strength = 0
Vpbe_ctor2:  xkappa = 0
Vpbe_ctor2:  Debye length = 0
Vpbe_ctor2:  zkappa2 = 0
Vpbe_ctor2:  zmagic = 7042.98
Vpbe_ctor2:  Constructing Vclist with 71 x 72 x 75 table
Vclist_ctor2:  Using 71 x 72 x 75 hash table
Vclist_ctor2:  automatic domain setup.
Vclist_ctor2:  Using 1.9 max radius
Vclist_setupGrid:  Grid lengths = (44.522, 45, 48.842)
Vclist_setupGrid:  Grid lower corner = (-10.574, 5.244, 29.902)
Vclist_assignAtoms:  Have 2681616 atom entries
Vacc_storeParms:  Surf. density = 10
Vacc_storeParms:  Max area = 191.134
Vacc_storeParms:  Using 1936-point reference sphere
Setting up PDE object...
Vpmp_ctor2:  Using meth = 2, mgsolv = 1
Setting PDE center to local center...
Vpmg_ctor2:  Filling boundary with old solution!
VPMG::focusFillBound -- New mesh mins = -16.264, -0.3645, 24.2465
VPMG::focusFillBound -- New mesh maxs = 39.638, 55.8525, 84.3995
VPMG::focusFillBound -- Old mesh mins = -18.8297, -3.04045, 20.193
VPMG::focusFillBound -- Old mesh maxs = 42.2037, 58.5284, 88.453
VPMG::extEnergy:  energy flag = 1
Vpmg_setPart:  lower corner = (-16.264, -0.3645, 24.2465)
Vpmg_setPart:  upper corner = (39.638, 55.8525, 84.3995)
Vpmg_setPart:  actual minima = (-18.8297, -3.04045, 20.193)
Vpmg_setPart:  actual maxima = (42.2037, 58.5284, 88.453)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
VPMG::extEnergy:   Finding extEnergy dimensions...
VPMG::extEnergy    Disj part lower corner = (-16.264, -0.3645, 24.2465)
VPMG::extEnergy    Disj part upper corner = (39.638, 55.8525, 84.3995)
VPMG::extEnergy    Old lower corner = (-18.8297, -3.04045, 20.193)
VPMG::extEnergy    Old upper corner = (42.2037, 58.5284, 88.453)
Vpmg_qmEnergy:  Zero energy for zero ionic strength!
VPMG::extEnergy: extQmEnergy = 0 kT
Vpmg_qfEnergyVolume:  Calculating energy
VPMG::extEnergy: extQfEnergy = 0 kT
VPMG::extEnergy: extDiEnergy = 0.179839 kT
Vpmg_fillco:  filling in source term.
fillcoCharge:  Calling fillcoChargeSpline2...
Vpmg_fillco:  filling in source term.
Vpmg_fillco:  marking ion and solvent accessibility.
fillcoCoef:  Calling fillcoCoefMol...
Vacc_SASA: Time elapsed: 0.276744
Vpmg_fillco:  done filling coefficient arrays
Vnm_tstop: stopping timer 27 (Setup timer).  CPU TIME = 1.086880e+00
Vnm_tstart: starting timer 28 (Solver timer)..
Vnm_tstart: starting timer 30 (Vmgdrv2: fine problem setup)..
Vbuildops: Fine: (097, 097, 129)
Vbuildops: Operator stencil (lev, numdia) = (1, 4)
Vnm_tstop: stopping timer 30 (Vmgdrv2: fine problem setup).  CPU TIME = 1.755500e-02
Vnm_tstart: starting timer 30 (Vmgdrv2: coarse problem setup)..
Vbuildops: Galer: (049, 049, 065)
Vbuildops: Galer: (025, 025, 033)
Vbuildops: Galer: (013, 013, 017)
Vnm_tstop: stopping timer 30 (Vmgdrv2: coarse problem setup).  CPU TIME = 2.171422e+00
Vnm_tstart: starting timer 30 (Vmgdrv2: solve)..
Vnm_tstop: stopping timer 40 (MG iteration).  CPU TIME = 8.244352e+01
Vprtstp: iteration = 0
Vprtstp: relative residual = 1.000000e+00
Vprtstp: contraction number = 1.000000e+00
Vprtstp: iteration = 1
Vprtstp: relative residual = 1.338085e-01
Vprtstp: contraction number = 1.338085e-01
Vprtstp: iteration = 2
Vprtstp: relative residual = 1.901154e-02
Vprtstp: contraction number = 1.420803e-01
Vprtstp: iteration = 3
Vprtstp: relative residual = 2.984374e-03
Vprtstp: contraction number = 1.569770e-01
Vprtstp: iteration = 4
Vprtstp: relative residual = 5.626097e-04
Vprtstp: contraction number = 1.885185e-01
Vprtstp: iteration = 5
Vprtstp: relative residual = 1.450069e-04
Vprtstp: contraction number = 2.577397e-01
Vprtstp: iteration = 6
Vprtstp: relative residual = 4.700901e-05
Vprtstp: contraction number = 3.241847e-01
Vprtstp: iteration = 7
Vprtstp: relative residual = 1.664390e-05
Vprtstp: contraction number = 3.540576e-01
Vprtstp: iteration = 8
Vprtstp: relative residual = 6.802247e-06
Vprtstp: contraction number = 4.086931e-01
Vprtstp: iteration = 9
Vprtstp: relative residual = 2.537231e-06
Vprtstp: contraction number = 3.729989e-01
Vprtstp: iteration = 10
Vprtstp: relative residual = 1.077934e-06
Vprtstp: contraction number = 4.248466e-01
Vprtstp: iteration = 11
Vprtstp: relative residual = 4.058141e-07
Vprtstp: contraction number = 3.764740e-01
Vnm_tstop: stopping timer 30 (Vmgdrv2: solve).  CPU TIME = 8.007201e+01
Vnm_tstop: stopping timer 28 (Solver timer).  CPU TIME = 8.395356e+01
Vpmg_setPart:  lower corner = (-16.264, -0.3645, 24.2465)
Vpmg_setPart:  upper corner = (39.638, 55.8525, 84.3995)
Vpmg_setPart:  actual minima = (-16.264, -0.3645, 24.2465)
Vpmg_setPart:  actual maxima = (39.638, 55.8525, 84.3995)
Vpmg_setPart:  bflag[FRONT] = 0
Vpmg_setPart:  bflag[BACK] = 0
Vpmg_setPart:  bflag[LEFT] = 0
Vpmg_setPart:  bflag[RIGHT] = 0
Vpmg_setPart:  bflag[UP] = 0
Vpmg_setPart:  bflag[DOWN] = 0
Vnm_tstart: starting timer 29 (Energy timer)..
Vpmg_energy:  calculating only q-phi energy
Vpmg_qfEnergyVolume:  Calculating energy
Vpmg_energy:  qfEnergy = 5.448772306771E+04 kT
Vnm_tstop: stopping timer 29 (Energy timer).  CPU TIME = 1.552700e-02
Vnm_tstart: starting timer 30 (Force timer)..
Vnm_tstop: stopping timer 30 (Force timer).  CPU TIME = 0.000000e+00
Vgrid_writeDX:  Opening virtual socket...
Vgrid_writeDX:  Writing to virtual socket...
Vgrid_writeDX:  Writing comments for ASC format.
printEnergy:  Performing global reduction (sum)
Vcom_reduce:  Not compiled with MPI, doing simple copy.
Vnm_tstop: stopping timer 26 (APBS WALL CLOCK).  CPU TIME = 1.629563e+02
##############################################################################
# MC-shell I/O capture file.
# Creation Date and Time:  Sat Dec 10 14:54:09 2022

##############################################################################
Vgrid_readDX:  Grid dimensions 97 x 97 x 129 grid
Vgrid_readDX:  Grid origin = (-16.264, -0.3645, 24.2465)
Vgrid_readDX:  Grid spacings = (0.582313, 0.585594, 0.469945)
Vgrid_readDX:  allocating 97 x 97 x 129 doubles for storage
