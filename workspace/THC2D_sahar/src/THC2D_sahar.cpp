//============================================================================
// Name        : THC2D_sahar.cpp
// Author      : B. Galvan
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================


#include "egsheader.h"
using namespace std;

int main() {

	omp_set_num_threads(12);
	char pause;

	///flags
	int g_flag = 1;
	int lte_flag = 0;

	///constants
	double g = g_flag * 9.806;										//gravitational acceleration
	double PI =  3.1415926559265;

	/// Domain dimensions = number of grid points in x and y direction
	int nx = 100;
	int ny = 100;
	double length = 600.0;
	double depth = 150.0;
	double dx = length/(nx-1);
	double dy = depth/(ny-1);
	double xmin=0;
	double xmax=length+xmin;
	double ymin=0;
	double ymax=length+ymin;

	///time
	double t = 0;
	double sec = 0;
	int hour =0;
	int tag=0;
	int hourcounter = 0;
	int coutcounter = 0;
	int savecounter = 0;
	int day = 3600*24;
	int mounth = 30*day;
	int tmax =20*12*mounth;
	double teq =mounth;									//time in which fracture and hydrostatic pressure become to a steady state
	double dt = double(7*day);

	/// Fluid parameters
	double viscw0 = 1.002e-3;								//viscosity of water
	double Kw = 2.37e9;									//Bulkmodulus of water
	double betaw = 1/Kw;									 //compressibility of water
	double cpw = 4200.0;									//heat capacity of water (J/kg °c)
	double rhow0 = 1000;								//density of water (kg/m^3)
	double tcw = 0.609;									//thermal conductivity of water (W/m °c)
	double tw0 = 200.0;									//injected water temperature °C
	double betaf = 10e-10;								//fluid compressibility from Nature paper of Steve

	double betaphi = 10e-8;								//crack compressibility from Nature paper of Steve
	double phiini = 0.10; 									//initial porosity of rock
	double phi0frac = 0.1;								//initial fracture porosity
	double phi_cap = 0.1; 									//initial porosity of rockn

	double phir= 0.01;
	double phia = -8e-9;									//pressure dependent porosity control variable

	double kappa0 = 4.845e-13;									//initial permeability
	double  perm_max = 1e-10;
	double c = 5;						  					//pressure dependent permeability control variable

	double s = 0.15;										//standard deviation for Kappa heterogenity
	double stc = 0.05;										//starndard deviation for TC Rock

	/// dt paremeters
	double kappatime1=0;									//max KappaP
	double kappatime2=0;									//max KappaT
	double kappatime3=0;									//save box
	double kappatime4=0;									//save box
	double kappatime5=0;									//save box
	double kappatime6=0;									//save box
	double kappatime7=0;									//save box
	double dffusvty_P = 0;									//pressure eq. diffusivity
	double dffusvty_T = 0;									// Temperature LTE eq. diffusivity
	double dffusvty_Tr = 0;									//Tr eq. diffusivity
	double dffusvty_Tf = 0;									//Tf eq. diffusivity
	double dffusvty_Cs = 0;									//Cs eq. diffusivity
	double dt_T =0;										//temparature dt

	///Temperature eq. parameters
	double cpr = 1170.0;								//heat capacity of rock (J/kg °c)
	double rhor0 = 2820.0;								//density of rock (kg/m^3)		//buntsandstein
	double tcr = 2.80;									//thermal conductivity of rock (W/m /°c) //3.20 thermal cond. of granite cho 2009
	double tr0 = 200.0;									//rock initial temperature  °C
	double htc = 1; 									//overall heat transfer coefficient (W/m^2 °c)
	double tgrad = 0.0;									//temperature gradient in rock

	///heat transfer area
	double rockr = 1e-3;									//radius of rock spheres
	int rockn = int((1-phiini)* dx *dy /(PI* rockr * rockr));			//number of rock spheres in one grid cell 2D	number of air bubbles in one gris cell
	double hta = 3*(1-phiini)/rockr *100;						//heat transfer area


	/// mass transfer eq. paremeters
	double csini = 0.00;
	double c_fresh =1e-20;								//freshwater salinity
	double rho_max = 1.2*rhow0;							//density of salt water in Henry's problem
	double dm = 3.565e-6;								//Simpsons or Acherer for Henry

	/// pressure eq. parameters
	double vwx0 = 6.6e-5;
	double pmax = 10e6;
	double pmin =0;//-1e6;
	double sigma_n = 120e6;

	///pumping hole
	double e_hole = 0.20*depth;
	double h_hole = 0.40*depth;
	double kappa_hole = 2e-15;
	double tc_hole = 10*tcr;								//????	//pressure dependent permeability control variable
	///CAP properties
	double k_cap = 1e-20;
	double phiini_cap = 0.01;
	double z_cap = ymin + 2*dy;
	///Distance from injection to calculate the Pw(depth )
	double d_pump= 200/dx;								//meters



	///save time parameters
	int tcout =6;				//cout/printing time in hour
	int tsave =2*24;				//saving in vtk format in hour
	int seccounter = 3600;		//how many seconds has an hour /OR change the TIME UNIT

	int 	timesave=10*3600;
	int 	timesave0=timesave;

	/// average variables
	double pwleft=0;
	double pwright=0;
	double pwmid=0;
	double value = 0;
	double trleft=0;
	double trright=0;
	double twleft=0;
	double twright=0;
	double ave_hta;

	/// min max
	double min_perm, max_perm;
	double min_phi, max_phi;
	double min_vx, max_vx;
	double min_vy, max_vy, ave_vy, ave_vx;
	double hta_min, hta_max;
	double min_rho, max_rho;
	double max_Re;
	double alpha = 0, gamma=1;

	createHostVariable(double, Y, nx*ny);
	createHostVariable(double, X, nx*ny);
	createHostVariable(double, Phi, nx*ny);
	createHostVariable(double, iPhi, nx*ny);
	createHostVariable(double, Pw, nx*ny);                              	               //fluid pressure of former time step
	createHostVariable(double, iPw, nx*ny);                                                 //last time step fluid pressure
	createHostVariable(double, Pw0, nx*ny);                                                 //relax fluid pressure
	createHostVariable(double, Kappa, nx*ny);
	createHostVariable(double, Kappa0, nx*ny);												///heterogenous permeability
	createHostVariable(double, KappaV, nx*ny);
	createHostVariable(double, Sxx, nx*ny);
	createHostVariable(double, Szz, nx*ny);
	createHostVariable(double, iTr, nx*ny);
	createHostVariable(double, Tr, nx*ny);													///rock temperature
	createHostVariable(double, TCm, nx*ny);													///matrix thermal conductivity
	createHostVariable(double, iTw, nx*ny);
	createHostVariable(double, Tw, nx*ny);
	createHostVariable(double, Q, nx*ny);													///heat transfer btw rock and fluid
	createHostVariable(double, Vwx, nx*ny);													///Darcy Velocity in x direction
	createHostVariable(double, Vwy, nx*ny);													///Darcy Velocity in y direction
	createHostVariable(double, Viscw, nx*ny);
	createHostVariable(double, Rhow, nx*ny);
	createHostVariable(double, iRhow, nx*ny);
	createHostVariable(double, Betaw, nx*ny);
	createHostVariable(double, Betaphi, nx*ny);
	createHostVariable(double, KappaTime, nx*ny);
	createHostVariable(double, Cs, nx*ny);													///salt concentration
	createHostVariable(double, iCs, nx*ny);
	createHostVariable(double, TCs, nx*ny);													///temporary concentration after diffusion before advection
	createHostVariable(double, Dm, nx*ny);
	createHostVariable(double, Hta, nx*ny);													///heat transfer area matrix
	createHostVariable(double, S, nx*ny);													///specific storativity of the porous medium [LT2 /M]
	createHostVariable(double, Cpw, nx*ny);
	createHostVariable(double, Frac, nx*ny);													///fracture network
	createHostVariable(double, Re, nx*ny);													///Reynold's number
	createHostVariable(double, Adv, nx*ny);													///depth dependent pressure at distance d from the source
	createHostVariable(double, Pd, 1*ny);
	createHostVariable(double, Rcp, nx*ny);													///rhowcp for lte
	min_perm=kappa0;
	max_perm=kappa0;
	min_phi=phiini;
	max_phi=phiini;
	double random1=0;


	///saving names
	string name;
	string Prx,Rhorx, Kapparx, Phirx,Cpwrx, Viscrx, Vwxrx, Vwyrx, Rerx, AdvCoeff,Tw2D, Pdepth, distave;
	name = "08-04-o-";
	distave ="08-04-o-";
	AdvCoeff = "Adv-08-04-o-";
	Tw2D = "Tw-08-04-o-";
	Pdepth= "Pd-08-04-o-";
	Prx ="PwRx-0804-o";
	Rhorx ="RhowRx-0804-o";
	Kapparx ="KappaRx-0804-o";
	Phirx ="PhiRx-0804-o";
	Cpwrx ="CpwRx-0804-o";
	Viscrx ="ViscRx-0804-o";
	Vwxrx ="VwxRx-0804-o";
	Vwyrx ="VwyRx-0804-o";
	Rerx ="ReRx-0804-o";

	///initializing the arrays
	readdata(Frac,nx, ny,"Soultz_Fractures_more_100x");
#pragma omp  for collapse(2)
	for (int j=0; j<ny; j++)
	{
		for (int i=0; i<nx; i++)
		{
			X[i+j*nx] = xmin+i*dx;
			Y[i+j*nx] = ymin+j*dy;
			Szz[i+j*nx] =  (rhor0*(1-Phi[i+j*nx])+iRhow[i+j*nx]*Phi[i+j*nx])*g*g_flag*Y[i+j*nx] ;                                     //density includes porosity
			Sxx[i+j*nx] =  0.7*Szz[i+j*nx];
			//			Frac[i+j*nx] =  0;
			// 			  	Phi[i+j*nx] = phir+(phiini- phir) * exp(phia * ( 0.5*(Sxx[i+j*nx]+Szz[i+j*nx]) ));
			Phi[i+j*nx] = phiini;//Frac[i+j*nx] * phi0frac + (1-Frac[i+j*nx] )*phiini;
			iPhi[i+j*nx] = phiini;//
			srand(time(NULL));
			Kappa0[i+j*nx] =kappa0;//+50*kappa0*Frac[i+j*nx];//* exp(-(random1*random1)/(2*s*s));;
			Kappa[i+j*nx] =Kappa0[i+j*nx];
			Tr[i+j*nx] = tr0+tgrad*Y[i+j*nx];
			iTr[i+j*nx] = Tr[i+j*nx];
			iTw[i+j*nx] = Tr[i+j*nx];
			Tw[i+j*nx] = Tr[i+j*nx];
			Cs[i+j*nx] = csini;
			iCs[i+j*nx] = csini;

			if(j==0 && X[i+j*nx] >=150 && X[i+j*nx]<=450){
				Cs[i+j*nx] = 1.0;
				iCs[i+j*nx] = 1.0;
			}

			Q[i+j*nx] = 0;
			Vwx[i+j*nx] = 1e-20;
			Vwy[i+j*nx] = 1e-20;
			//Viscw[i+j*nx] = 0.6612 * pow((Tw[i+j*nx]+ 273 - 229), -1.562);		//viscosity of water as a function of temperature
			Viscw[i+j*nx] = 10E-3;//Elders problem
			Rhow[i+j*nx] = 1E3;// Elders problem 9.992e2 + 9.539e-2 *Tw[i+j*nx]  + 7.999e-1 * Cs[i+j*nx]*1e3 -7.618e-3*Tw[i+j*nx]*Tw[i+j*nx] - 2.409e-3* Cs[i+j*nx]*1e3*Tw[i+j*nx];
			iRhow[i+j*nx] =  Rhow[i+j*nx];
			Pw[i+j*nx] = Rhow[i+j*nx]*g*g_flag*Y[i+j*nx]+10;
			iPw[i+j*nx] = Rhow[i+j*nx]*g*g_flag*Y[i+j*nx]+10;
			Betaw[i+j*nx] = betaw;
			Betaphi[i+j*nx] = betaphi;
			Dm[i+j*nx] = dm;
			Cpw[i+j*nx] = cpw;
			Hta[i+j*nx] =0;
			S[i+j*nx] = betaphi*(1-Phi[i+j*nx] ) + betaw*Phi[i+j*nx] ;
			TCm[i+j*nx] = (1-Frac[i+j*nx] )*((1-Phi[i+j*nx])*tcr + Phi[i+j*nx]*tcw)+Frac[i+j*nx]*tcw ;
			Rcp[i+j*nx] = (1-Frac[i+j*nx] )*((1-Phi[i+j*nx])*rhor0*cpr + Phi[i+j*nx]*tcw)+Frac[i+j*nx]*Rhow[i+j*nx]*Cpw[i+j*nx];
			Re[i+j*nx] = 1.0;//Rhow[i+j*nx]*(g_flag *Vwy[i+j*nx] +(1-g_flag) *Vwx[i+j*nx] )* 2*rockr/(Viscw[i+j*nx] *Phi [i+j*nx]);
			Adv[i+j*nx] = 1;
		}
	}
	//top and botoom cap
/*
#pragma omp  for collapse(2)
	for (int j=0; j<2 ; j++)
	{
		for (int i=0; i<nx; i++)
		{
			Kappa0[i+j*nx] = k_cap;
			Kappa0[i+(ny-1-j)*nx] = k_cap;
		}
	}
*/

	double rhocp=0,valxlt=0,valylt=0, valxnlt=0, valynlt=0;


	///determining dt

	StateEqTemperetureDensityTransport(Viscw , Rhow, Tw, Pw,Cs, Betaw,Cpw, rhow0,nx,ny);
	double isomass;

	for(int i=0; i< nx; i++)
	{										//hydrostatic pressure
		for(int j=0; j< ny; j++)
		{
			isomass = Rhow[i+j*nx]*Phi[i+j*nx];
			iPw[i+j*nx] += isomass*fabs(g)*dy;
		}
		isomass=0;
	}
	for(int i=0; i< nx; i++)
	{										//hydrostatic pressure
		for(int j=0; j< ny; j++)
		{
			Szz[i+j*nx] =  (rhor0*(1-Phi[i+j*nx])+Rhow[i+j*nx]*Phi[i+j*nx])*g*g_flag*Y[i+j*nx];
			Sxx[i+j*nx] =  0.7*Szz[i+j*nx];
		}
	}

	//PorosityFunction(Phi, iPhi, Pw, Sxx, Szz, phiini, phir, phia,nx, ny);
	//PermeabilitySinglePhase(Kappa, Phi, Kappa0,Pw, Sxx, Szz,  phiini, z_cap, k_cap,  c, h_hole, e_hole,  kappa_hole, dx,  dy,  nx, ny);

	TimeStepCalculation(dt, Phi, Kappa, Viscw, Betaphi, Betaw, Rhow, TCm, Dm,Cpw, tcr, tcw,  cpr, rhor0, dx, dy, nx,  ny);
	min_perm = *min_element(Kappa,Kappa+(nx*ny));
	max_perm = *max_element(Kappa,Kappa+(nx*ny));
	min_phi = *min_element(Phi,Phi+(nx*ny));
	max_phi = *max_element(Phi,Phi+(nx*ny));

	cout<<"min_perm "<<min_perm <<"       max_perm "<<max_perm<<"     min_phi "<<min_phi <<"       max_phi "<<max_phi<<endl;
	cout<<"dt: "<<dt<<"  nx: "<<nx<<"  ny: "<<ny<<"  dx: "<<dx<<"  dy: "<<dy<<setw(14)<<"  grainsize: "<<rockr<<endl;
	cout<< "t("<<seccounter<<"sec)"<<"\t"<<"dt(s)"<<"\t"<<setw(14)<<"Tr inj." <<"\t"<<setw(14)<<"Tw inj."<<"\t"<<setw(14)<<"pw left"<<"\t"<<setw(14)<<"pw right"<<"\t"<<setw(14)<<"min vx"<<"\t"<<setw(14)<<"min vy"<<endl<<endl;

	SaveVtk(Pw,Cs,Tw,Tr, Kappa, Rhow, Vwx,Vwy, name, dx, dy, hour/24, nx, ny);

	sec = 0;
	t = 0;
	//Pf relaxation loop
	cout<<"Begin Pf relaxation"<<endl;

	double maxpf=10e10;
	double rpfdiff=10e10;
	double lpfdiff=10e10;
	double maxvx=-10;
	double maxvy=-10;
	double meanvx=0;
	double meanvy=0;
	int loopcounter=0;

/*
	while( maxpf>0.5)
	{
		timesave=timesave-dt;

		//density includes porosity
		PressureRelaxFunction(Pw,Vwx, Vwy,iPw, Rhow, Phi, Kappa, Viscw ,Betaw, Betaphi, S,  vwx0 ,dx, dy, t, g_flag, dt, nx, ny);
		StateEqTemperetureDensityTransport(Viscw , Rhow, Tw, Pw,Cs, Betaw,Cpw, rhow0, nx,ny);
		PorosityFunction(Phi, iPhi, Pw, Sxx, Szz, phiini, phir, phia,nx, ny);
		PermeabilitySinglePhase(Kappa, Phi, Kappa0,Pw, Sxx, Szz,  phiini, z_cap, k_cap,  c, h_hole, e_hole,  PI, dx,  dy,  nx, ny);

		maxpf = -10;
		maxvx = -10;
		maxvy = -10;

		for(int j=1; j< ny-1; j++)
		{
			for(int i=1; i< nx-1; i++)
			{
				meanvx=+Vwx[i+j*nx];
				meanvy=+Vwy[i+j*nx];
				maxpf=max(maxpf,fabs(Pw[i+j*nx]-iPw[i+j*nx]));
			}
		}
		min_perm = *min_element(Kappa,Kappa+(nx*ny));
		max_perm = *max_element(Kappa,Kappa+(nx*ny));
		min_phi = *min_element(Phi,Phi+(nx*ny));
		max_phi = *max_element(Phi,Phi+(nx*ny));
		maxvx=*max_element(Vwx, Vwx+(nx*ny));
		maxvy=*max_element(Vwy, Vwy+(nx*ny));

		meanvx=meanvx/(nx*ny);
		meanvy=meanvy/(nx*ny);
		MemoryCopy(iPw,Pw,nx,ny);
		MemoryCopy(iTw, Tw, nx, ny);
		MemoryCopy(iTr, Tr, nx, ny);
		if(timesave<0){
			loopcounter++;
			cout<<"maxpf= "<<maxpf<<"   rightMid= "<<rpfdiff<<"   leftMid= "<<lpfdiff<<"  meanvx= "<<meanvx<<"   meanvy= "<<meanvy<<"  min_perm "<<min_perm <<"  max_perm "<<max_perm<<"   min_phi "<<min_phi <<"   max_phi "<<max_phi<<endl;
			timesave=timesave0;
		}
	}
*/

	//     	writedata(Pw, nx, ny,hour,Prx);
	// 	writedata(Rhow, nx, ny,hour,Rhorx);
	// 	writedata(Kappa, nx, ny,hour,Kapparx);
	// 	writedata(Phi, nx, ny,hour,Phirx);
	// 	writedata(Cpw, nx, ny,hour,Cpwrx);
	// 	writedata(Viscw, nx, ny,hour,Viscrx);
	// 	writedata(Vwx, nx, ny,hour,Vwxrx);
	// 	writedata(Vwy, nx, ny,hour,Vwyrx);
	// 	writedata(Re, nx, ny,hour,Rerx);
	cout<<"maxpf= "<<maxpf<<"   rightMid= "<<rpfdiff<<"   leftMid= "<<lpfdiff<<"  meanvx= "<<meanvx<<"   meanvy= "<<meanvy<<"  min_perm "<<min_perm <<"  max_perm "<<max_perm<<"   min_phi "<<min_phi <<"   max_phi "<<max_phi<<endl;
	cout<<"End Pf relaxation"<<endl;
	sec = 0;
	t = 0;

	//         readdata(Pw,nx, ny, "PwRx-0804-v0");
	// 	readdata(Pw0,nx, ny, "PwRx-0804-v0");
	//         readdata(Rhow,nx, ny, "RhowRx-0804-v0");
	//         readdata(Kappa,nx, ny, "KappaRx-0804-v0");
	//         readdata(Phi,nx, ny, "PhiRx-0804-v0");
	//         readdata(Cpw,nx, ny, "CpwRx-0804-v0");
	//         readdata(Viscw,nx, ny, "ViscRx-0804-v0");
	//         readdata(Vwx,nx, ny, "VwxRx-0804-v0");
	//         readdata(Vwy,nx, ny, "VwyRx-0804-v0");
	//         readdata(Re,nx, ny, "ReRx-0804-v0");




/*	for(int i=0; i< nx; i++)
	{										//hydrostatic pressure
		for(int j=0; j< ny; j++)
		{
			Pw0[i+j*nx] = Pw[i+j*nx] ;
		}
	}

	SaveDistAve(Pw0,Pw,Tw,Tr,Cs, Rhow,Phi,Cpw,Adv, distave , dx,dy, hour,  nx, ny);
	*/

	///main time loop
	while (t <= tmax)
	{
		min_perm = 1000;
		max_perm = 0;
		min_vx=1000;
		min_vy=1000;
		max_rho = 0;
		max_Re = 0;
		ave_vy = 0;
		ave_vx = 0;

		//TimeStepCalculation(dt, Phi, Kappa, Viscw, Betaphi, Betaw, Rhow, TCm, Dm,Cpw,  tcr, tcw, cpr, rhor0, dx, dy, nx,  ny);
		StateEqTemperetureDensityTransport(Viscw , Rhow, Tw, Pw,Cs, Betaw,Cpw, rhow0,nx,ny);
		MemoryCopy(iRhow, Rhow, nx, ny);
		//PorosityFunction(Phi,iPhi, Pw,Sxx, Szz, phiini,phir, phia, nx,  ny);
		//PermeabilitySinglePhase(Kappa, Phi, Kappa0,Pw, Sxx, Szz,  phiini, z_cap, k_cap,  c, h_hole, e_hole,  kappa_hole, dx,  dy,  nx, ny);
		EGSPressureFunction(Pw,Vwx, Vwy,iPw, Rhow, Phi, Kappa, Viscw ,Betaw, Betaphi, S, Re, Pw0,Y, rockr ,  h_hole, e_hole, pmax,pmin,vwx0 ,dx, dy, t, g_flag, dt, nx, ny);
		MemoryCopy(iPw, Pw, nx, ny);
		DensityTransportSML(Cs,iCs,Dm,Phi,Vwx,Vwy,X,Y, Rhow, h_hole, e_hole,c_fresh,dx, dy,  t, teq, dt, nx, ny);
		MemoryCopy(iCs, Cs, nx, ny);

		//     		ThermalConductivityFunction(TCm,Frac, Rcp, Phi, Rhow, Cpw,  tcr,  rhor0, cpr, tcw, nx, ny);
		//        		LTEDensityTransportSML(Tr, iTr, Tw, Vwx,  Vwy,Phi,  Rhow, Rcp,Cpw, TCm,X, Y, h_hole, e_hole,  tw0, dx, dy,  t, dt,nx, ny);
		//     		MemoryCopy(iTw, Tw, nx, ny);
		//     		MemoryCopy(iTr, Tr, nx, ny);
		//


		//LTNENonIsoThermal(Tw, Tr, iTw, iTr, Vwx, Vwy, Phi, Rhow, Cpw, tcw,   rhor0,   cpr,  tcr, rockr, X, Y, tw0,  dx, dy, dt,  nx, ny);
		//MemoryCopy(iTw, Tw, nx, ny);
		//MemoryCopy(iTr, Tr, nx, ny);

		t = t+dt;
		sec = sec+dt;
		if (sec >= seccounter)
		{
			min_perm = *min_element(Kappa,Kappa+(nx*ny));
			max_perm = *max_element(Kappa,Kappa+(nx*ny));
			min_vx = *min_element(Vwx,Vwx+(nx*ny));
			min_vy = *min_element(Vwy,Vwy+(nx*ny));
			max_Re = *max_element(Re,Re+(nx*ny));
			AverageFunction(ave_vx, Vwx, nx, ny);
			AverageFunction(ave_vy, Vwy, nx, ny);


			hourcounter += int(sec/seccounter);
			hour += int(sec/seccounter);
			coutcounter += int(sec/seccounter);
			savecounter += int(sec/seccounter);
			sec -=int(sec/seccounter)*seccounter;
			tag=int(hour);
			if (coutcounter >= tcout)
			{							//the txt file with radial ave. values
				cout<<"  min perm "<<min_perm<<"    max perm "<<max_perm<<"   Max Re "<<	max_Re<<endl;
				cout<<tag<<" day \t"<<setw(14)<<dt<<"\t" <<setw(14)<<pwleft<<"\t"<<setw(14)<<pwright<<"\t"<<setw(14)<<min_vx<<"\t"<<setw(14)<<min_vy<<endl;
				SaveDistAve(Pw0,Pw,Tw,Tr,Cs, Rhow,Phi,Cpw,Adv,distave ,dx,dy, tag,  nx, ny);
				coutcounter -=int(coutcounter/tcout)*tcout;
			}
			if (savecounter >= tsave)
			{							//saving in VTK format   and Pw(depth)
				cout<<endl<<"vtk time "<<"  day:  "<<tag<<", week :  "<<int(tag/7)+1<<endl;
				//velocity matrix
				SaveVtk(Pw,Cs,Tw,Tr,Kappa, Rhow,Vwx, Vwy,name, dx, dy,tag,nx,ny);
				//				PdWrite(Pd, Pw,Pw0, d_pump ,nx, ny, dy, hour,Pdepth);
				savecounter-=int(savecounter/tsave)*tsave;
			}
		}
	}

	writedata(Pw, nx, ny,hour,Prx);
	writedata(Rhow, nx, ny,hour,Rhorx);
	writedata(Kappa, nx, ny,hour,Kapparx);
	writedata(Phi, nx, ny,hour,Phirx);
	writedata(Cpw, nx, ny,hour,Cpwrx);
	writedata(Viscw, nx, ny,hour,Viscrx);
	writedata(Vwx, nx, ny,hour,Vwxrx);
	writedata(Vwy, nx, ny,hour,Vwyrx);
	writedata(Re, nx, ny,hour,Rerx);

	///free the pointer allocations
	DeleteHostVariable(X);
	DeleteHostVariable(Y);
	DeleteHostVariable(Phi);
	DeleteHostVariable(iPhi);
	DeleteHostVariable(Pw);
	DeleteHostVariable(Pw0);
	DeleteHostVariable(iPw);
	DeleteHostVariable(Kappa);
	DeleteHostVariable(Kappa0);
	DeleteHostVariable(KappaV);
	DeleteHostVariable(Sxx);
	DeleteHostVariable(Szz);
	DeleteHostVariable(iTr);
	DeleteHostVariable(Tr);
	DeleteHostVariable(TCm);
	DeleteHostVariable(iTw);
	DeleteHostVariable(Tw);
	DeleteHostVariable(Q);
	DeleteHostVariable(Vwx);
	DeleteHostVariable(Vwy);
	DeleteHostVariable(Viscw);
	DeleteHostVariable(Rhow);
	DeleteHostVariable( iRhow);
	DeleteHostVariable(Betaw);
	DeleteHostVariable(Betaphi);
	DeleteHostVariable(KappaTime);
	DeleteHostVariable(Cs);
	DeleteHostVariable(iCs);
	DeleteHostVariable(TCs);
	DeleteHostVariable(Dm);
	DeleteHostVariable(Hta);
	DeleteHostVariable(S);
	DeleteHostVariable(Cpw);
	DeleteHostVariable(Frac);
	DeleteHostVariable(Re);
	DeleteHostVariable(Adv);
	DeleteHostVariable(Pd);
	DeleteHostVariable(Rcp);
	cout<<"Alles klar!"<<endl;
}
