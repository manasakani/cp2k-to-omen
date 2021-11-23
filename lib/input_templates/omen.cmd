/*********************************************************************************************************************************
Definition of a wire structure:

1) the order of the parameters must not be changed
2) the points in the structure mat_coord(i,j) and ox_coord(i,j) must follow the order described below
3) comments must remain on one line or start again with //

Commands:

CB_Bandstructure	= conduction bandstructure of contacts
VB_Bandstructure	= valence bandstructure of contacts
CB_Transmission_UWF	= conduction band transmission with UMFPACK
VB_Transmission_UWF	= valence band transmission with UMFPACK
CB_Transmission_SWF	= conduction band transmission with SuperLU_DIST
VB_Transmission_SWF	= valence band transmission with SuperLU_DIST
CB_Transmission_MWF	= conduction band transmission with MUMPS
VB_Transmission_MWF	= valence band transmission with MUMPS
CB_Transmission_RGF	= conduction band transmission with recursive GF
VB_Transmission_RGF	= valence band transmission with recursive GF
EL_SC_UWF		= self-consistent electron simulation with UMFPACK
HO_SC_UWF		= self-consistent hole simulation with UMFPACK
EL_SC_SWF		= self-consistent electron simulation with SuperLU_DIST
HO_SC_SWF		= self-consistent hole simulation SuperLU_DIST
EL_SC_MWF		= self-consistent electron simulation with MUMPS
HO_SC_MWF		= self-consistent hole simulation MUMPS
EL_SC_RGF		= self-consistent electron simulation with recursive GF
HO_SC_RGF		= self-consistent hole simulation with recursive GF
Write_Layer_Matrix      = write file with atom positions + connections
Write_Grid_Matrix       = write file with grid points + file with index of atom position

*********************************************************************************************************************************/
/*Parameters*/

mat_name            	= mat_par;
lattice_type		= lattice_dat;
a0                      = 0.1;                //lattice constant
first_atom              = anion;               	//atom situated at [0 0 0]

poisson_solver		= 0;
//full_current		= 1;
transport_type		= 0;

//atomic_dos		= 1;

read_hamiltonian	= 1;
NDim			= 3;

injection_type		= [3 500];
//injection_parallel	= 0;

tb			= 10;                   	//tight-binding order
dsp3			= 0;                   	//passivation energy [eV]

Temp           		= 300;				//operation temperature

n_of_modes		= 16;              		//number of modes for bandstructure calculation
Nk			= 51;                           //number of k points in bandstructure calculation. High Nk important for good energy grid in pMOS
//bs_solver		= full;

max_bond_def       = 0.001;

last_first		= 0;                    	//last super cell is equal to first super cell if last_first=1

eta_res			= 0;				//imaginary part of the energy in the reservoir (if 0 less time to get BC), eta_res<=i*1e-6
eta			= 0;			//imaginary part of the energy in the device

Elimit			= 50e-3;                	//energy interval after a mode where the energy grid is finer (should not be changed)
Emin_tail		= 5.0e-3;			//energy interval below the lowest band (should not be changed)
EOffset                 = 10*UT;			//Emax = Emin + EOffset or Emax = max(Efl,Efr) + EOffset
dE_in			= 1.0e-3;			//smallest energy interval (should not be changed)
dE_f			= 1.0e-3;			//largest energy interval (should not be changed)
dE_sep			= 1.0e-4;			//distance between a channel turn-on and the following energy point (should not be changed)
NEmax			= 100;				//maximum number of calculated energy point

CPU_ppoint		= 2;
CPU_per_bc		= 8;
CPU_interleave		= 1;
spec_decomp		= 1;		

x                       = [1 0 0];			//transport direction
y                       = [0 1 0];			//direction of confinement
z     			= [0 0 1];              	//direction of confinement

NVG			= 1;
Vgmin			= 0.0;
Vgmax			= 0.0;

NVS			= 1;
Vsmin			= 0.0;
Vsmax			= 0.0;

NVD			= 1;
Vdmin 		 = 0.5;
Vdmax 		 = 0.5;

restart 		 = [2 0 0 0];
vact_file 		 = vact_dat;

update_energy		= 2;
energy_file		= E_dat;

update_cell_dim		= 1;
cell_file		= Smin_dat;

update_fermi		= 1;
fermi_level 		 = 2;

/*********************************************************************************************************************************/
/*Structure*/

sort_coordinates	= 0;
grid_accuracy		= 1;                    	//number of grid points added between two neighbor atoms (2 is a good value)

no_mat			= 1;				//number of pieces that form the nanowire (channel + oxide)
no_channel_mat          = 1;                    	//number of pieces that form the nanowire channel
no_oxide_mat            = 0;                    	//number of pieces that form the oxide around the wire  

Lc			=   2.288;
Ls			=  3.0;  			//source length
Ld			=  3.0;			//drain length
tc			=   1.642;
hc			=   1.639;
x0			=  -0.001;
z0			=  -0.017;

mat_type(1)		= square;			//type of material: square, circle, or sphere
mat_cs(1)		= yes;                          //does the material determine the nanowire cross section 
mat_coord(1,1)	        = [0.0 x0 z0];          	//[xmin ymin zmin]
mat_coord(1,2)		= [Ls+Lc+Ld x0 z0];		//[xmax ymin zmin]
mat_coord(1,3)		= [Ls+Lc+Ld tc z0];		//[xmax ymax zmin]
mat_coord(1,4)		= [0.0 tc z0];			//[xmin ymax zmin]
mat_coord(1,5)		= [0.0 x0 hc];			//[xmin ymin zmax]
mat_coord(1,6)		= [Ls+Lc+Ld x0 hc];		//[xmax ymin zmax]
mat_coord(1,7)		= [Ls+Lc+Ld tc hc];		//[xmax ymax zmax]
mat_coord(1,8)		= [0.0 tc hc];			//[xmin ymax zmax]

/*********************************************************************************************************************************/

command(1)		= Write_Layer_Matrix;
command(2)		= Write_Grid_Matrix;
//command(3)		= CB_Bandstructure;
//command(4)		= VB_Bandstructure;
command(5)		= EL_SC_SSWF;
