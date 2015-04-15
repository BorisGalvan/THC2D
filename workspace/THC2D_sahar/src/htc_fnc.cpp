



#include "egsheader.h"

template <typename T>
void MemoryCopy(T *in, T *out, int nx, int ny)
{
#pragma omp parallel for
	for(int j=0; j< ny; j++)
	{
		for(int i=0; i< nx; i++)
		{
			in[i+j*nx] = out[i+j*nx];

		}
	}
}

template <typename T>
void PorosityFunction( T *Phi,T *iPhi, T *Pw, T *Sxx, T *Szz, T phiini, T phir,  T phia,int nx, int ny)
{
#pragma omp parallel for
	for(int j=0; j< ny; j++)
	{
		for(int i=0; i< nx; i++)
		{
			Phi[i+j*nx] =phir+(iPhi[i+j*nx] -phir)* exp(phia * ( 0.5*(Sxx[i+j*nx]+Szz[i+j*nx]-2* Pw[i+j*nx]+(Szz[i+j*nx]-Sxx[i+j*nx])*cos(2*3.14159/9))));			//rutqvist tsang 2003 (9)
		}
	}
}

template <typename T>
void PermeabilitySinglePhase(T *Kappa,  T *Phi, T *Kappa0,T *Pw, T *Sxx, T *Szz,  T phiini, T cap_depth, T k_cap, T c,T h, T e,T PI, T dx, T dy,  int nx, int ny)
{
#pragma omp parallel for
	for(int j=0; j< ny; j++)
	{
		for(int i=0; i< nx; i++)
		{
			Kappa[i+j*nx] = Kappa0[i+j*nx] * exp(c * (Phi[i+j*nx]/phiini -1.0));				//Rutqvist-Tsang2002 , Kolditz
		}
	}
#pragma omp  for collapse(2)
	for (int j=0; j<2 ; j++)
	{
		for (int i=0; i<nx; i++)
		{
			Kappa0[i+j*nx] = k_cap;
			Kappa0[i+(ny-1-j)*nx] = k_cap;
		}
	}
}
template <typename T>
void ThermalConductivityFunction(T *TCm, T *Frac, T *Rcp, T *Phi, T *Rhow, T *Cpw, T tcr, T rhor0, T cpr, T  tcw, int nx, int ny)
{
#pragma omp parallel for
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			TCm[i+j*nx] =(1-Frac[i+j*nx] )*((1-Phi[i+j*nx])*tcr + Phi[i+j*nx]*tcw)+Frac[i+j*nx]*tcw ;
			Rcp[i+j*nx] = (1-Frac[i+j*nx] )*((1-Phi[i+j*nx])*rhor0*cpr + Phi[i+j*nx]*Cpw[i+j*nx]*Rhow[i+j*nx])+Frac[i+j*nx]*Rhow[i+j*nx]*Cpw[i+j*nx];
		}
	}
}

template <typename T>
void StateEqTemperetureDensityTransport(T *Viscw , T *Rhof, T *Tw, T *Pw, T *Cs, T *Betaw, T *Cpw, T rhof0, int nx, int ny)
{
	double a=rhof0;
#pragma omp parallel for
	for(int j=0; j< ny; j++)
	{
		for(int i=0; i< nx; i++)
		{
			//Viscw[i+j*nx] = 2.4e-5*pow(10,248.37/(Tw[i+j*nx]+133.15));						//Huyakorn & Pinder, 1977

			//Rhof[i+j*nx] = 9.992e2 + 9.539e-2 *Tw[i+j*nx]  -7.618e-3*Tw[i+j*nx]*Tw[i+j*nx]+ 3.130e-5*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx]  -6.174e-8*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] + 4.337e-1 * Pw[i+j*nx]*1e-6 +2.549e-5* Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx] - 2.899e-7*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] + 9.578e-10 *Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] +1.763e-3*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6 -1.231e-4 *Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx] +1.366e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]  +4.045e-9*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] -1.467e-5*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6 +8.839e-7 *Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx] -1.102e-9 *Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx] +4.247e-11 *Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] - 3.959e-14 *Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] *Tw[i+j*nx] + 7.999e-1 * Cs[i+j*nx]*1e3 - 2.409e-3* Cs[i+j*nx]*1e3*Tw[i+j*nx] +2.580e-5* Cs[i+j*nx]*1e3*Tw[i+j*nx]*Tw[i+j*nx]  -6.856e-8 *Cs[i+j*nx]*1e3*Tw[i+j*nx]*Tw[i+j*nx] *Tw[i+j*nx] -6.297e-4*Cs[i+j*nx]*1e3*Pw[i+j*nx]*1e-6 +9.363e-7*Pw[i+j*nx]*1e-6*Pw[i+j*nx]*1e-6*Cs[i+j*nx]*1e3;

			//Elders problem
			Rhof[i+j*nx] = rhof0+200*Cs[i+j*nx];// from  rhof0+((rhofs-rhof)/(Cs-Co))*Cs[i+j*nx]

			// Cpw[i+j*nx]  = (4.193 -2.273e-4 *Tw[i+j*nx] + 2.369e-6*Tw[i+j*nx]*Tw[i+j*nx] + 1.670e-10*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx]*Tw[i+j*nx] - 3.978e-3*Pw[i+j*nx] *1e-6 +3.229e-5*Pw[i+j*nx] *1e-6 *Tw[i+j*nx] -1.072e-9 *Pw[i+j*nx] *1e-6 *Tw[i+j*nx] *Tw[i+j*nx] *Tw[i+j*nx] + 1.913e-5 *Pw[i+j*nx]*Pw[i+j*nx]*1e-12 - 4.175e-7*Pw[i+j*nx]*Pw[i+j*nx]*1e-12*Tw[i+j*nx] +2.306e-9 *Pw[i+j*nx]*Pw[i+j*nx]*1e-12 *Tw[i+j*nx] *Tw[i+j*nx] )*1e3;  //Cs(Tw, Pw)


		}
	}
}

template <typename T>
void LTEDensityTransportSML(T *Tr, T *iTr, T *Tw, T *Vwx,  T *Vwy, T *Phi,  T *Rhow, T *Rcp, T *Cpw,T* TCr, T *x,T *y,T h,T e,T tw_fresh, T dx, T dy,T t,  T dt,int nx, int ny)
{
	T m, m1, m2, m3, m4;
	createHostVariable(T,tTr, nx*ny);

	T dtdy = dt/(dy*dy);
	T dtdx = dt/(dx*dx);

#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{
		for(int i=1; i<nx-1; i++)
		{
			tTr[i+j*nx] = iTr[i+j*nx] +dtdx/Rcp[i+j*nx]  * 0.5 * (TCr[i+j*nx] + TCr[(i+1)+j*nx])/*m1*/*(iTr[(i+1)+j*nx]-iTr[i+j*nx] ) +dtdx/ Rcp[i+j*nx] * 0.5 * (TCr[i+j*nx] + TCr[(i-1)+j*nx])/*m3*/*(iTr[(i-1)+j*nx]-iTr[i+j*nx] ) + dtdy/Rcp[i+j*nx] * 0.5 * (TCr[i+j*nx] + TCr[i+(j+1)*nx])/*m2*/*(iTr[i+(j+1)*nx]-iTr[i+j*nx] )+dtdy/Rcp[i+j*nx]* 0.5 * (TCr[i+j*nx] + TCr[i+(j-1)*nx])/*m4*/*(iTr[i+(j-1)*nx]-iTr[i+j*nx] );

		}
	}

	for(int j=1 ; j<ny-1; j++)
	{
		for(int i=1; i<nx-1; i++)
		{
			if (tTr[i+j*nx]+1 <tw_fresh)
			{
				tTr[i+j*nx] =iTr[i+j*nx];
				cout<<"Tr too cold!   "<<i<<"  "<<j<<"   tTr:  "<<tTr[i+j*nx]<<"   iTr:  "<<iTr[i+j*nx] <<endl;
			}
			if (tTr[i+j*nx] > 200 )
			{
				tTr[i+j*nx] =iTr[i+j*nx];
				cout<<"Tr too hot!   "<<i<<"  "<<j<<"   Tr:  "<<tTr[i+j*nx] <<"  "<<endl;
			}
		}
	}
	///// boundary conditions EGS first order Neumann
#pragma omp parallel for									//for bottom j=ny-1 ; top j=0 Neumann
	for(int i=1; i< nx-1; i++)
	{
		tTr[i+(ny-1)*nx] = tTr[i+(ny-2)*nx];
		tTr[i] = tTr[i+nx];
	}
#pragma omp parallel for								// boundary conditions for left i=0; right i=nx-1  Neumann
	for(int j=0; j< ny; j++)
	{
		tTr[j*nx]=  tTr[1+j*nx];
		tTr[nx-1+j*nx] =  tTr[nx-2+j*nx];
	}


	for(int j=int(h/dy); j<= int((e+h)/dy) ; j++)
	{			//bore hole   Dirichlet
		tTr[j*nx]=tw_fresh;
	}
	tTr[int(nx/2)+int(ny/2)*nx]=tw_fresh;
	T npox, npoy;
	int yind[4];
	int xind[4];
	T ivalx[4];
	T ivaly[4];
	T ival,valx,valy;
	T alphal;
#pragma omp parallel private(npox,npoy,yind,xind,ivalx,ivaly,ival,valx,valy,alphal)
	{
#pragma omp for
		for(int j=1; j<ny-1; j++)
		{
			for(int i=1; i<nx-1; i++)
			{

				valx=Vwx[i+j*nx]*Cpw[i+j*nx] * Rhow[i+j*nx]/Rcp[i+j*nx];
				valy=Vwy[i+j*nx]*Cpw[i+j*nx] * Rhow[i+j*nx]/Rcp[i+j*nx];

				//advection
				if(fabs(valx)>1.0e-10 || fabs(valy)>1.0e-10){
					//get new positions
					npox=x[i+j*nx]-dt*valx;
					npoy=y[i+j*nx]-dt*valy;

					//interpolation I need 16 points
					InterpolationPoints<T>(xind,yind,npox,npoy,x,y,nx,ny);
					for(int jj =0; jj<4;jj++)
					{
						for(int ii =0; ii<4;ii++)
						{
							ivalx[ii]=tTr[xind[ii]+yind[jj]*nx];
						}
						ivaly[jj]=spline_cubic(ivalx,(npox-x[xind[1]+yind[1]*nx])/dx);
					}
					ival=spline_cubic(ivaly,(npoy-y[xind[1]+yind[1]*nx])/dy);
					Tr[i+j*nx]=ival;
				}
			}
		}
	}
	///// boundary conditions EGS first order Neumann
#pragma omp parallel for
	for(int i=0; i< nx; i++)
	{
		Tr[i+(ny-1)*nx] = Tr[i+(ny-2)*nx] ;
		Tr[i] = Tr[i+nx];
	}
#pragma omp parallel for
	for(int j=1; j< ny-1; j++)
	{
		Tr[j*nx]=  Tr[1+j*nx];
		Tr[nx-1+j*nx] =  Tr[nx-2+j*nx];
	}

	Tr[int(nx/2)+int(ny/2)*nx]=tw_fresh;

#pragma omp parallel for
	for(int j=0; j<ny; j++)
	{
		for(int i=0; i<nx; i++)
		{
			Tw[i+j*nx] = Tr[i+j*nx];
		}
	}
	DeleteHostVariable(tTr);

}

template <typename T>
void LTNENonIsoThermal(T *Tw, T *Tr,  T *iTw, T *iTr, T *Vwx, T *Vwy, T *Phi, T *Rhow, T *Cpw, T tcw,  T rhor0, T  cpr, T tcr, T rockr, T *x, T *y,T tw_fresh, T dx, T dy, T dt, int nx, int ny)
{
	T  m1, m2, m3, m4, m5, m6, m7, m8, m; 	//fluid part
	T a1, a2;				//rock
	T dtdy = dt /(dy*dy)* tcw;
	T dtdx = dt/(dx*dx)* tcw;
	T hadT = 0;
	a1 = dt/(dx*dx)* tcr /rhor0 /cpr ;
	a2 = dt/(dy*dy) * tcr /rhor0 /cpr;
	char pause;
	createHostVariable(T,tTw, nx*ny);
	T Nu_fs= 0;
	T Pr = 7;					//Prandlt no. for water
	T beta = 10;				//for spherical particles
	T h_star = 1;
	T w =10;



#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{											//explicit temperature diffusion eq
		for(int i=1; i<nx-1; i++)
		{
			//				hadT=dt  *3*(1-Phi[i+j*nx])/rockr *w* (iTr[i+j*nx] -iTw[i+j*nx]); 			//~1.7e4 *w

			tTw[i+j*nx] = iTw[i+j*nx] +dtdx / (Rhow[i+j*nx] *Cpw[i+j*nx]) * (iTw[(i+1)+j*nx]+iTw[(i-1)+j*nx]-2*iTw[i+j*nx] )+dtdy/(Rhow[i+j*nx] *Cpw[i+j*nx] ) * (iTw[i+(j+1)*nx]+iTw[i+(j-1)*nx]-2*iTw[i+j*nx]) + dt  *3*(1-Phi[i+j*nx])/rockr*w * (iTr[i+j*nx] -iTw[i+j*nx])/(Rhow[i+j*nx]*Cpw[i+j*nx]*Phi[i+j*nx]);

			Tr[i+j*nx] = iTr[i+j*nx] + a1*(iTr[(i+1)+j*nx]-2*iTr[i+j*nx]+iTr[(i-1)+j*nx]) +a2*(iTr[i+(j+1)*nx]-2*iTr[i+j*nx] +iTr[i+(j-1)*nx]) - dt  *3*(1-Phi[i+j*nx])/rockr *w* (iTr[i+j*nx] -iTw[i+j*nx])/((1-Phi[i+j*nx])*rhor0 *cpr);
		}
	}

	// boundary conditions for bottom j=ny-1 ; top j=0
#pragma omp parallel for
	for(int i=1; i< nx-1; i++)
	{
		//bottom
		tTw[i+(ny-1)*nx] = iTw[i+(ny-1)*nx] +dtdx /Rhow[i+(ny-1)*nx] /Cpw[i+(ny-1)*nx] /*m1*/*(iTw[(i+1)+(ny-1)*nx]+iTw[(i-1)+(ny-1)*nx]-2*iTw[i+(ny-1)*nx] )+dtdy /Rhow[i+(ny-1)*nx] /Cpw[i+(ny-1)*nx] /*m2*/ *(2*iTw[i+(ny-2)*nx]-2*iTw[i+(ny-1)*nx]) +dt*3*(1-Phi[i+(ny-1)*nx])/rockr *w* (iTr[i+(ny-1)*nx] - iTw[i+(ny-1)*nx])/(Rhow[i+(ny-1)*nx]*Cpw[i+(ny-1)*nx]*Phi[i+(ny-1)*nx])/*Qdtw */;
		//top
		tTw[i] = iTw[i] +dtdx /Rhow[i] /Cpw[i] /*m1*/*(iTw[(i+1)]+iTw[(i-1)]-2*iTw[i] )+dtdy /Rhow[i] /Cpw[i] /*m2*/ *(2*iTw[i+nx]-2*iTw[i]) +dt*3*(1-Phi[i])/rockr *w*  (iTr[i] - iTw[i])/(Rhow[i]*Cpw[i]*Phi[i])/*Qdtw */;
	}
	// boundary conditions for left i=0; right i=nx-1
#pragma omp parallel for
	for(int j=0; j< ny-1; j++)
	{
		//left
		tTw[j*nx] = iTw[j*nx] +dtdx /Rhow[j*nx] /Cpw[j*nx] /*m1*/*(2*iTw[1+j*nx]-2*iTw[j*nx] )+dtdy/Rhow[j*nx] /Cpw[j*nx] /*m2*/ *(iTw[(j+1)*nx]+iTw[(j-1)*nx]-2*iTw[j*nx]) +dt *3*(1-Phi[j*nx])/rockr *w*  (iTr[j*nx] - iTw[j*nx])/(Rhow[j*nx]*Cpw[j*nx]*Phi[j*nx])/*Qdtw */;
		//right
		tTw[nx-1+j*nx] =  iTw[nx-1+j*nx] +dtdx /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m1*/*(2*iTw[(nx-2)+j*nx]-2*iTw[nx-1+j*nx] )+dtdy /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m2*/ *(iTw[nx-1+(j+1)*nx]+iTw[nx-1+(j-1)*nx]-2*iTw[nx-1+j*nx]) +dt  *3*(1-Phi[nx-1+j*nx])/rockr *w *  (iTr[nx-1+j*nx] - iTw[nx-1+j*nx])/(Rhow[nx-1+j*nx]*Cpw[nx-1+j*nx]*Phi[nx-1+j*nx])/*Qdtw */;

	}
	// boundary conditions for corners
	//top left
	tTw[0] =  iTw[0] +dtdx /Rhow[0] /Cpw[0] *2*(iTw[1]-iTw[0] )+dtdy /Rhow[0] /Cpw[0] *2 *(iTw[nx]-iTw[0]) +dt*3*(1-Phi[0])/rockr *w* (iTr[0] - iTw[0])/(Rhow[0]*Cpw[0]*Phi[0]);
	//top right
	tTw[nx-1] =  iTw[nx-1] +dtdx /Rhow[nx-1] /Cpw[nx-1] *2*(iTw[nx-2]-iTw[nx-1] )+dtdy /Rhow[nx-1] /Cpw[nx-1] *2*(iTw[nx-1+nx]-iTw[nx-1]) +dt*3*(1-Phi[nx-1])/rockr *w*(iTr[nx-1] - iTw[nx-1])/(Rhow[nx-1]*Cpw[nx-1]*Phi[nx-1]);
	//bottom left
	tTw[(ny-1)*nx] = iTw[(ny-1)*nx] + dtdx /Rhow[(ny-1)*nx] /Cpw[(ny-1)*nx]*2*(iTw[1+(ny-1)*nx]-iTw[(ny-1)*nx] )+dtdy /Rhow[(ny-1)*nx] /Cpw[(ny-1)*nx]  *2 *(iTw[(ny-2)*nx] - iTw[(ny-1)*nx]) +dt*3*(1-Phi[(ny-1)*nx])/rockr *w* (iTr[(ny-1)*nx] - iTw[(ny-1)*nx])/(Rhow[(ny-1)*nx]*Cpw[(ny-1)*nx]*Phi[(ny-1)*nx]);
	//bottom right
	tTw[nx-1+(ny-1)*nx] = iTw[nx-1+(ny-1)*nx] +dtdx /Rhow[nx-1+(ny-1)*nx] /Cpw[nx-1+(ny-1)*nx]*2 *(iTw[nx-2+(ny-1)*nx]-iTw[nx-1+(ny-1)*nx] )+dtdy /Rhow[nx-1+(ny-1)*nx] /Cpw[nx-1+(ny-1)*nx]*2 *(iTw[nx-1+(ny-2)*nx]-iTw[nx-1+(ny-1)*nx]) +dt*3*(1-Phi[nx-1+(ny-1)*nx])/rockr *w* (iTr[nx-1+(ny-1)*nx] - iTw[nx-1+(ny-1)*nx])/(Rhow[nx-1+(ny-1)*nx]*Cpw[nx-1+(ny-1)*nx]*Phi[nx-1+(ny-1)*nx]);


	tTw[int(nx/2)+int(ny/2)*nx]=tw_fresh;
	////////advection part
	T npox, npoy;
	int yind[4];
	int xind[4];
	T ivalx[4];
	T ivaly[4];
	T ival,valx,valy;
	T alphal;
#pragma omp parallel private(npox,npoy,yind,xind,ivalx,ivaly,ival,valx,valy,alphal)
	{
#pragma omp for
		for(int j=1; j<ny-1; j++)
		{
			for(int i=1; i<nx-1; i++){
				valx=Vwx[i+j*nx]/Phi[i+j*nx];
				valy=Vwy[i+j*nx]/Phi[i+j*nx];

				//get new positions
				npox=x[i+j*nx]-dt*valx;
				npoy=y[i+j*nx]-dt*valy;

				//interpolation I need 16 points
				InterpolationPoints<T>(xind,yind,npox,npoy,x,y,nx,ny);
				for(int jj =0; jj<4;jj++)
				{

					for(int ii =0; ii<4;ii++)
					{
						ivalx[ii]=tTw[xind[ii]+yind[jj]*nx];
					}
					ivaly[jj]=spline_cubic(ivalx,(npox-x[xind[1]+yind[1]*nx])/dx);
				}
				ival=spline_cubic(ivaly,(npoy-y[xind[1]+yind[1]*nx])/dy);
				Tw[i+j*nx]=ival;
			}
		}
	}

	// boundary conditions for left i=0; right i=nx-1
#pragma omp parallel for
	for(int j=1; j< ny-1; j++)
	{
		// left von Neumann i=0 fluid
		Tw[j*nx] =  iTw[j*nx] +dtdx /Rhow[j*nx] /Cpw[j*nx] *(2*iTw[1+j*nx]-2*iTw[j*nx] )+dtdy/Rhow[j*nx] /Cpw[j*nx] *(iTw[(j+1)*nx]+iTw[(j-1)*nx]-2*iTw[j*nx]) +dt *3*(1-Phi[j*nx])/rockr *w*  (iTr[j*nx] - iTw[j*nx])/(Rhow[j*nx]*Cpw[j*nx]*Phi[j*nx]);//Tw[1+j*nx];
		Tr[j*nx] = iTr[j*nx] + a1*2*(iTr[1+j*nx]-iTr[j*nx]) +a2*(iTr[(j+1)*nx]-2*iTr[j*nx] +iTr[(j-1)*nx]) +dt *3*(1-Phi[j*nx])/rockr *w * (iTw[j*nx] - iTr[j*nx])/((1-Phi[j*nx])*rhor0 *cpr);
		// right von Neumann i=nx-1 fluid
		Tw[nx-1+j*nx] =  iTw[nx-1+j*nx] +dtdx /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m1*/*(2*iTw[(nx-2)+j*nx]-2*iTw[nx-1+j*nx] )+dtdy /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m2*/ *(iTw[nx-1+(j+1)*nx]+iTw[nx-1+(j-1)*nx]-2*iTw[nx-1+j*nx]) +dt  *3*(1-Phi[nx-1+j*nx])/rockr *w *  (iTr[nx-1+j*nx] - iTw[nx-1+j*nx])/(Rhow[nx-1+j*nx]*Cpw[nx-1+j*nx]*Phi[nx-1+j*nx]);//Tw[nx-2+j*nx];
		Tr[nx-1+j*nx] =  iTr[nx-1+j*nx] + a1*2*(-iTr[nx-1+j*nx]+iTr[nx-2+j*nx]) +a2*(iTr[nx-1+(j+1)*nx]-2*iTr[nx-1+j*nx] +iTr[nx-1+(j-1)*nx]) +dt*3*(1-Phi[nx-1+j*nx])/rockr *w  * (iTw[nx-1+j*nx] - iTr[nx-1+j*nx])/((1-Phi[nx-1+j*nx])*rhor0 *cpr);
	}
	// boundary conditions for bottom j=ny-1 ; top j=0
#pragma omp parallel for
	for(int i=0; i< nx; i++)
	{
		// bottom von Neumann  j=ny-1 fluid
		Tw[i+(ny-1)*nx] = Tw[i+(ny-2)*nx] ;
		Tr[i+(ny-1)*nx] =Tr[i+(ny-2)*nx] ;
		// top von Neumann  j=0 fluid
		Tw[i] = Tw[i+nx];
		Tr[i] = iTr[i+1+nx];
	}

	// 		Tr[int(nx/2)+int(ny/2)*nx]=tw_fresh;
	Tw[int(nx/2)+int(ny/2)*nx]=tw_fresh;

	for(int j=0 ; j<ny; j++)
	{											//explicit temperature diffusion eq
		for(int i=0; i<nx; i++)
		{
			if (Tw[i+j*nx] <tw_fresh )
			{
				cout<<"Tw too cold!   "<<i<<"  "<<j<<"  "<<"   tTw:  "<<tTw[i+j*nx]<<"   iTw:  "<<iTw[i+j*nx] <<endl;
				pause=getchar();
			}
			if (Tw[i+j*nx]>200)
			{
				cout<<"Tw too hot!   "<<Tw[i+j*nx]-200<<"  "<<i<<"  "<<j<<"  "<<"    iTr-iTw:   "<< iTr[i+j*nx]-iTw[i+j*nx]<<"    Tr-tTw:   "<< Tr[i+j*nx]-tTw[i+j*nx]<<"  "<<"   Tr- Tw:   "<<Tr[i+j*nx]-Tw[i+j*nx] <<" Qdtw:   "<<"    cpw= "<<Cpw[i+j*nx]<<endl;
				pause=getchar();
			}
			if (Tr[i+j*nx] <tw_fresh)
			{
				cout<<"Tr too cold!   "<<i<<"  "<<j<<"   200-Tr:  "<<200-Tr[i+j*nx]<<"   iTr:  "<<iTr[i+j*nx]<<"   Tw:  "<<iTw[i+j*nx]<<"  " <<endl;
				pause=getchar();
			}
			if (Tr[i+j*nx] >200)
			{
				cout<<"Tr too hot!   "<<i<<"  "<<j<<"   Tr:  "<<Tr[i+j*nx]<<"   iTr:  "<<iTr[i+j*nx]<<"   Tw:  "<<iTw[i+j*nx]<<"  " <<endl;
				pause=getchar();
			}
		}
	}
	DeleteHostVariable(tTw);
}

template <typename T>
void LTNENonIsoThermalTF(T *Tw, T *iTw, T *iTr,T *Vwx,T *Vwy,T *Phi,T *Rhow,T *Hta,T*Cpw,T t,T tcw,T *x,T *y,T h,T e,T tw_fresh, T dx,T dy,T dt,int nx,int ny)
{
	T Qdtw, m1, m2, m3, m4, m5, m6, m7, m8, m; 	//fluid part
	T dtdy = dt /(dy*dy);
	T dtdx = dt/(dx*dx);
	Qdtw = 0;
	char pause;
	createHostVariable(T,tTw, nx*ny);


#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{											//explicit temperature diffusion eq
		for(int i=1; i<nx-1; i++)
		{
			tTw[i+j*nx] = iTw[i+j*nx] +dtdx* tcw /Rhow[i+j*nx] /Cpw[i+j*nx] /*m1*/*(iTw[(i+1)+j*nx]+iTw[(i-1)+j*nx]-2*iTw[i+j*nx] )+dtdy* tcw /Rhow[i+j*nx] /Cpw[i+j*nx] /*m2*/ *(iTw[i+(j+1)*nx]+iTw[i+(j-1)*nx]-2*iTw[i+j*nx]) +dt  *Hta[i+j*nx] *  (iTr[i+j*nx] - iTw[i+j*nx])/(Rhow[i+j*nx]*Cpw[i+j*nx]*Phi[i+j*nx])/*Qdtw */;
		}
	}
	// boundary conditions for bottom j=ny-1 ; top j=0
#pragma omp parallel for
	for(int i=1; i< nx-1; i++)
	{
		//bottom
		tTw[i+(ny-1)*nx] = iTw[i+(ny-1)*nx] +dtdx* tcw /Rhow[i+(ny-1)*nx] /Cpw[i+(ny-1)*nx] /*m1*/*(iTw[(i+1)+(ny-1)*nx]+iTw[(i-1)+(ny-1)*nx]-2*iTw[i+(ny-1)*nx] )+dtdy* tcw /Rhow[i+(ny-1)*nx] /Cpw[i+(ny-1)*nx] /*m2*/ *(2*iTw[i+(ny-2)*nx]-2*iTw[i+(ny-1)*nx]) +dt*Hta[i+(ny-1)*nx] *  (iTr[i+(ny-1)*nx] - iTw[i+(ny-1)*nx])/(Rhow[i+(ny-1)*nx]*Cpw[i+(ny-1)*nx]*Phi[i+(ny-1)*nx])/*Qdtw */;
		//top
		tTw[i] = iTw[i] +dtdx* tcw /Rhow[i] /Cpw[i] /*m1*/*(iTw[(i+1)]+iTw[(i-1)]-2*iTw[i] )+dtdy* tcw /Rhow[i] /Cpw[i] /*m2*/ *(2*iTw[i+nx]-2*iTw[i]) +dt*Hta[i]  *  (iTr[i] - iTw[i])/(Rhow[i]*Cpw[i]*Phi[i])/*Qdtw */;

	}

#pragma omp parallel for
	for(int j=1; j< ny-1; j++)
	{
		//left
		tTw[j*nx] = iTw[j*nx] +dtdx* tcw /Rhow[j*nx] /Cpw[j*nx] /*m1*/*(2*iTw[1+j*nx]-2*iTw[j*nx] )+dtdy* tcw /Rhow[j*nx] /Cpw[j*nx] /*m2*/ *(iTw[(j+1)*nx]+iTw[(j-1)*nx]-2*iTw[j*nx]) +dt *Hta[j*nx] *  (iTr[j*nx] - iTw[j*nx])/(Rhow[j*nx]*Cpw[j*nx]*Phi[j*nx])/*Qdtw */;
		//right
		tTw[nx-1+j*nx] =  iTw[nx-1+j*nx] +dtdx* tcw /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m1*/*(2*iTw[(nx-2)+j*nx]-2*iTw[nx-1+j*nx] )+dtdy* tcw /Rhow[nx-1+j*nx] /Cpw[nx-1+j*nx] /*m2*/ *(iTw[nx-1+(j+1)*nx]+iTw[nx-1+(j-1)*nx]-2*iTw[nx-1+j*nx]) +dt  *Hta[nx-1+j*nx] *  (iTr[nx-1+j*nx] - iTw[nx-1+j*nx])/(Rhow[nx-1+j*nx]*Cpw[nx-1+j*nx]*Phi[nx-1+j*nx])/*Qdtw */;

	}
	// boundary conditions for corners
	//top left
	tTw[0] =  iTw[0] +dtdx* tcw /Rhow[0] /Cpw[0] *2*(iTw[1]-iTw[0] )+dtdy* tcw /Rhow[0] /Cpw[0] *2 *(iTw[nx]-iTw[0]) +dt*Hta[0]  * (iTr[0] - iTw[0])/(Rhow[0]*Cpw[0]*Phi[0]);
	//top right
	tTw[nx-1] =  iTw[nx-1] +dtdx* tcw /Rhow[nx-1] /Cpw[nx-1] *2*(iTw[nx-2]-iTw[nx-1] )+dtdy* tcw /Rhow[nx-1] /Cpw[nx-1] *2*(iTw[nx-1+nx]-iTw[nx-1]) +dt*Hta[nx-1]* (iTr[nx-1] - iTw[nx-1])/(Rhow[nx-1]*Cpw[nx-1]*Phi[nx-1]);
	//bottom left
	tTw[(ny-1)*nx] = iTw[(ny-1)*nx] + dtdx* tcw /Rhow[(ny-1)*nx] /Cpw[(ny-1)*nx]*2*(iTw[1+(ny-1)*nx]-iTw[(ny-1)*nx] )+dtdy* tcw /Rhow[(ny-1)*nx] /Cpw[(ny-1)*nx]  *2 *(iTw[(ny-2)*nx] - iTw[(ny-1)*nx]) +dt*Hta[(ny-1)*nx] * (iTr[(ny-1)*nx] - iTw[(ny-1)*nx])/(Rhow[(ny-1)*nx]*Cpw[(ny-1)*nx]*Phi[(ny-1)*nx]);
	//bottom right
	tTw[nx-1+(ny-1)*nx] = iTw[nx-1+(ny-1)*nx] +dtdx* tcw /Rhow[nx-1+(ny-1)*nx] /Cpw[nx-1+(ny-1)*nx]*2 *(iTw[nx-2+(ny-1)*nx]-iTw[nx-1+(ny-1)*nx] )+dtdy* tcw /Rhow[nx-1+(ny-1)*nx] /Cpw[nx-1+(ny-1)*nx]*2 *(iTw[nx-1+(ny-2)*nx]-iTw[nx-1+(ny-1)*nx]) +dt*Hta[nx-1+(ny-1)*nx] * (iTr[nx-1+(ny-1)*nx] - iTw[nx-1+(ny-1)*nx])/(Rhow[nx-1+(ny-1)*nx]*Cpw[nx-1+(ny-1)*nx]*Phi[nx-1+(ny-1)*nx]);

	////////advection part
	T npox, npoy;
	int yind[4];
	int xind[4];
	T ivalx[4];
	T ivaly[4];
	T ival,valx,valy;
	T alphal;
#pragma omp parallel private(npox,npoy,yind,xind,ivalx,ivaly,ival,valx,valy,alphal)
	{
#pragma omp for
		for(int j=1; j<ny-1; j++)
		{
			for(int i=1; i<nx-1; i++){
				valx=Vwx[i+j*nx]/Phi[i+j*nx];
				valy=Vwy[i+j*nx]/Phi[i+j*nx];

				//         			//advection
				//         			if(fabs(valx)>1.0e-10 || fabs(valy)>1.0e-10)
				//         			{
				//get new positions
				npox=x[i+j*nx]-dt*valx;
				npoy=y[i+j*nx]-dt*valy;

				//interpolation I need 16 points
				InterpolationPoints<T>(xind,yind,npox,npoy,x,y,nx,ny);
				for(int jj =0; jj<4;jj++)
				{

					for(int ii =0; ii<4;ii++)
					{
						ivalx[ii]=tTw[xind[ii]+yind[jj]*nx];
					}
					ivaly[jj]=spline_cubic(ivalx,(npox-x[xind[1]+yind[1]*nx])/dx);
				}
				ival=spline_cubic(ivaly,(npoy-y[xind[1]+yind[1]*nx])/dy);
				Tw[i+j*nx]=ival;
				//         			}
		}
		}
	}

	// boundary conditions for bottom j=ny-1 ; top j=0
#pragma omp parallel for
	for(int i=0; i< nx; i++)
	{
		// bottom von Neumann  j=ny-1 fluid
		Tw[i+(ny-1)*nx] = Tw[i+(ny-2)*nx] ;
		// top von Neumann  j=0 fluid
		Tw[i] = Tw[i+nx];
	}
	// boundary conditions for left i=0; right i=nx-1
#pragma omp parallel for
	for(int j=1; j<= ny-1; j++)
	{
		// left von Neumann i=0 fluid
		Tw[j*nx] = Tw[1+j*nx];
		// right von Neumann i=nx-1 fluid
		Tw[nx-1+j*nx] = Tw[nx-2+j*nx];
	}

	//         	#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{											//explicit temperature diffusion eq
		for(int i=1; i<nx-1; i++)
		{
			//     		if(fabs(iTr[i+j*nx] - iTw[i+j*nx])<1e-18)
			//     		{}
			//     		else
			//     		{
			Qdtw =  dt  * Hta[i+j*nx] * (iTr[i+j*nx] -iTw[i+j*nx])/(Rhow[i+j*nx]*Cpw[i+j*nx]*Phi[i+j*nx]);
			// //     		}
			if(Qdtw<0)
			{
				//         			cout<<"Qdtw <0= "<<Qdtw<<"  "<<i<<"  "<<j<<"    iTr -  iTw= "<< iTr[i+j*nx]-iTw[i+j*nx]<<"     iTr- tTw:   "<<iTr[i+j*nx]-tTw[i+j*nx] <<"    t= "<<int(t/dt)<<"    cpw= "<<Cpw[i+j*nx]<<"   Hta= "<<Hta[i+j*nx]<<endl;
				//         			pause=getchar();
				Qdtw=0;
			}
			if (Tw[i+j*nx] <tw_fresh )
			{
				cout<<"Tw too cold!   "<<i<<"  "<<j<<"  "<<"   tTw:  "<<tTw[i+j*nx]<<"   iTw:  "<<iTw[i+j*nx] <<endl;
				pause=getchar();
			}
			if (Tw[i+j*nx]>200)
			{
				cout<<"Tw too hot!   "<<"  "<<i<<"  "<<j<<"  t:  "<<int(t/dt)<<"  "<<"    iTr-iTw:   "<< iTr[i+j*nx]-iTw[i+j*nx]<<"    iTr-tTw:   "<< iTr[i+j*nx]-tTw[i+j*nx]<<"  "<<"   iTr- Tw:   "<<iTr[i+j*nx]-Tw[i+j*nx] <<" Qdtw:   "<<"    cpw= "<<Cpw[i+j*nx]<< Qdtw<<endl;
				pause=getchar();
			}
		}
	}
	Tw[int(nx/2)+int(ny/2)*nx]=tw_fresh;

	DeleteHostVariable(tTw);

}

template <typename T>
void LTNENonIsoThermalTR(T *iTw, T *Tr,  T *iTr, T *Phi, T *Hta, T  t, T rhor0, T  cpr, T tcr, T h,T e,T tw_fresh,T dx, T dy, T dt, int nx, int ny)
{
	T Qdtr,a1, a2, a3, a4, a5, a6, a7, a8,a;				//rock
	char pause;

	a1 = dt/(dx*dx)* tcr /rhor0 /cpr ;
	a2 = dt/(dy*dy) * tcr /rhor0 /cpr;

#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{											//explicit temperature diffusion eq
		for(int i=1; i<nx-1; i++)
		{
			Tr[i+j*nx] = iTr[i+j*nx] + a1*(iTr[(i+1)+j*nx]-2*iTr[i+j*nx]+iTr[(i-1)+j*nx]) +a2*(iTr[i+(j+1)*nx]-2*iTr[i+j*nx] +iTr[i+(j-1)*nx]) +dt * Hta[i+j*nx] * (iTw[i+j*nx] - iTr[i+j*nx])/((1-Phi[i+j*nx])*rhor0 *cpr);
		}
	}
	// boundary conditions for bottom j=ny-1 ; top j=0
#pragma omp parallel for
	for(int i=0; i< nx; i++)
	{
		//bottom
		Tr[i+(ny-1)*nx] =Tr[i+(ny-2)*nx] ;
		//top
		Tr[i] = iTr[i+1+nx];// + a1*(iTr[(i+1)]-2*iTr[i]+iTr[(i-1)]) +a2*2*(iTr[i+nx]-iTr[i] ) +dt * Hta[i] * (iTw[i] - iTr[i])/((1-Phi[i])*rhor0 *cpr);

	}
	// boundary conditions for left i=0; right i=nx-1
#pragma omp parallel for
	for(int j=1; j< ny-1; j++)
	{
		//left
		Tr[j*nx] = iTr[j*nx] + a1*2*(iTr[1+j*nx]-iTr[j*nx]) +a2*(iTr[(j+1)*nx]-2*iTr[j*nx] +iTr[(j-1)*nx]) +dt * Hta[j*nx] * (iTw[j*nx] - iTr[j*nx])/((1-Phi[j*nx])*rhor0 *cpr);
		//right
		Tr[nx-1+j*nx] =  iTr[nx-1+j*nx] + a1*2*(-iTr[nx-1+j*nx]+iTr[nx-2+j*nx]) +a2*(iTr[nx-1+(j+1)*nx]-2*iTr[nx-1+j*nx] +iTr[nx-1+(j-1)*nx]) +dt * Hta[nx-1+j*nx] * (iTw[nx-1+j*nx] - iTr[nx-1+j*nx])/((1-Phi[nx-1+j*nx])*rhor0 *cpr);
	}

	for(int j=0 ; j<ny; j++)
	{											//explicit temperature diffusion eq
		for(int i=0; i<nx; i++)
		{
			//     			if(iTr[i+j*nx] - iTw[i+j*nx]<1e-20){Qdtr=0;}
			//     			else{
			Qdtr = dt * Hta[i+j*nx] * (iTw[i+j*nx] - iTr[i+j*nx])/((1-Phi[i+j*nx])*rhor0 *cpr);
			if(Qdtr>0)
			{
				//         				cout<<"Qdtr>0     "<<Qdtr<<"  "<<i<<"  "<<j<<"  "<<"   iTw- iTr:   "<<iTw[i+j*nx]- iTr[i+j*nx] <<"    Hta="<<Hta[i+j*nx]<<endl;
				//         				pause=getchar();
				Qdtr=0;
			}
			//         			}
			if (Tr[i+j*nx] <tw_fresh)
			{
				cout<<"Tr too cold!   "<<i<<"  "<<j<<"   200-Tr:  "<<200-Tr[i+j*nx]<<"   iTr:  "<<iTr[i+j*nx]<<"   Tw:  "<<iTw[i+j*nx]<<"  " <<Qdtr<<endl;
				pause=getchar();
			}
			if (Tr[i+j*nx] >200)
			{
				cout<<"Tr too hot!   "<<i<<"  "<<j<<"   Tr:  "<<Tr[i+j*nx]<<"   iTr:  "<<iTr[i+j*nx]<<"   Tw:  "<<iTw[i+j*nx]<<"  " <<Qdtr<<endl;
				pause=getchar();
			}
		}
	}
	Tr[int(nx/2)+int(ny/2)*nx]=tw_fresh;

}


template <typename T>
void Hfs(T *Hta ,T *TCr,T *Phi, T *Re, T tcw,T rockr ,int nx, int ny)
{
	char pause;
	T Nu_fs= 0;
	T Pr = 7;					//Prandlt no. for water
	T beta = 10;				//for spherical particles
	T h_star = 1;
	T w =1;
#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{											//explicit temperature diffusion eq
		for(int i=1; i<nx-1; i++)
		{
			Nu_fs=  (0.255/Phi[i+j*nx] ) *pow(Pr, 1/3)*pow(10,2/3);
			h_star=2*rockr*1/(0.255/Phi[i+j*nx] *pow(Pr, 1/3)*pow(10,2/3)*tcw) + 1/(beta*TCr[i+j*nx]);
			Hta[i+j*nx]=3*(1-Phi[i+j*nx])/rockr *w; 			//~1.7e4 *w
		}
	}

}



template <typename T>
void EGSPressureFunction(T *Pw,T *Vwx, T *Vwy,T *iPw, T *Rhof ,T *Phi, T *Kappa, T *Viscw ,T *Betaf, T *Betaphi, T *S, T *Re, T *Pw0, T * Y, T rockr , T h,  T e,T pmax, T pmin,T vwx0 ,T dx, T dy, T t, int g_flag, T dt,int nx, int ny)
{
	char pause;
	T m1, m2, m3, m4, m5, m6, m7, m8;
	T dtdy = dt/(dy*dy);
	T dtdx = dt/(dx*dx);
	T g = 9.806;
	T isomass =0;
	T isomass2 =0;

#pragma omp parallel for
for(int j=1; j< ny-1; j++)
{
	for(int i=1; i< nx-1; i++)
	{


		Pw[i+j*nx]= dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx] + Kappa[(i+1)+j*nx]/ Viscw[(i+1)+j*nx])*iPw[(i+1)+j*nx]+dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[(i-1)+j*nx]/ Viscw[(i-1)+j*nx])*iPw[(i-1)+j*nx]+ dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j+1)*nx]/ Viscw[i+(j+1)*nx])*iPw[i+(j+1)*nx]+dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j-1)*nx]/ Viscw[i+(j-1)*nx])*iPw[i+(j-1)*nx] +(1-dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx] + Kappa[(i+1)+j*nx]/ Viscw[(i+1)+j*nx])- dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j+1)*nx]/ Viscw[i+(j+1)*nx])-dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[(i-1)+j*nx]/ Viscw[(i-1)+j*nx])-dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j-1)*nx]/ Viscw[i+(j-1)*nx])) * iPw[i+j*nx] - dt/dy * 0.5 /(Viscw[i+j*nx] *S[i+j*nx])* g*g_flag * Kappa[i+j*nx] *(Rhof[i+(j+1)*nx] - Rhof[i+(j-1)*nx])  -dt/dy * 0.5 /(Viscw[i+j*nx] *S[i+j*nx])* g*g_flag * Rhof[i+j*nx] *( Kappa[i+(j+1)*nx] -  Kappa[i+(j-1)*nx]) -dt/dy * 0.5 /(S[i+j*nx])* g*g_flag * Rhof[i+j*nx]* Kappa[i+j*nx]  *(1/Viscw[i+(j+1)*nx]-1/Viscw[i+(j-1)*nx]);


		Vwx[i+j*nx] =-0.5 *Kappa[i+j*nx]/Viscw[i+j*nx]*(iPw[(i+1)+j*nx] - iPw[(i-1)+j*nx])/dx;
		Vwy[i+j*nx] =-0.5 *(Kappa[i+j*nx]/Viscw[i+j*nx])*(iPw[i+(j+1)*nx]-iPw[i+(j-1)*nx])/dy + g_flag *(Kappa[i+j*nx]/Viscw[i+j*nx])*Rhof[i+j*nx]*g ;
	}
}

for(int j=1; j< ny-1; j++)
{
	for(int i=1; i< nx-1; i++)
	{
		Re[i+j*nx] = Rhof[i+j*nx]*max(fabs(Vwy[i+j*nx]) ,fabs(Vwx[i+j*nx]) )* 2*rockr/(Viscw[i+j*nx] *Phi [i+j*nx]);
		if(Pw[i+j*nx]<0)
		{
			cout<<" day "<<t/(3600*24)<<" P =   "<<Pw[i+j*nx]<<"  ("<<i+1<<"  "<<j+1<<")"<<"  "<<"  "<<m1<<"  "<<m2<<"  "<<m3<<"  "<<m4<<"  "<<m5<<"  "<<m6<<"  "<<m7<<endl;
			pause=getchar();
			Pw[i+j*nx]=iPw[i+j*nx];
		}
	}
}
/// Pressure boundary conditions
#pragma omp parallel for
for(int i=1; i< nx-1; i++)
{
	//top
	Pw[i] = 0.0;// Elders problem Pw[i+nx];//-g_flag *Rhof[i+nx]*g*dy;
	// 		bottom
	Pw[i+(ny-1)*nx] = Pw[i+(ny-2)*nx];//+g_flag *Rhof[i+(ny-1)*nx]*g*dy;
}

#pragma omp parallel for

for(int j=1; j< ny-1 ; j++)
{

	Pw[j*nx] = 2*dtdx /S[j*nx] *0.5*(Kappa[j*nx]/Viscw[j*nx] + Kappa[1+j*nx]/Viscw[1+j*nx] )/*m1*/*iPw[1+j*nx]+ dtdy /S[j*nx] *0.5*(Kappa[j*nx]/ Viscw[j*nx]+ Kappa[(j+1)*nx]/Viscw[(j+1)*nx])/*m2*/*iPw[(j+1)*nx]+dtdy /S[j*nx] *0.5*(Kappa[j*nx]/Viscw[j*nx]+ Kappa[(j-1)*nx]/Viscw[(j-1)*nx])/*m4*/*iPw[(j-1)*nx] +(1-2*dtdx /S[j*nx] *0.5*(Kappa[j*nx]/Viscw[j*nx] + Kappa[1+j*nx]/Viscw[1+j*nx] )/*m1*/-dtdy /S[j*nx] *0.5*(Kappa[j*nx]/ Viscw[j*nx]+ Kappa[(j+1)*nx]/Viscw[(j+1)*nx])/*m2*/-dtdy /S[j*nx] *0.5*(Kappa[j*nx]/Viscw[j*nx]+ Kappa[(j-1)*nx]/Viscw[(j-1)*nx])/*m4*/) * iPw[j*nx] - dt/dy * 0.5 /(Viscw[j*nx] *S[j*nx])* g*g_flag * Kappa[j*nx] *(Rhof[(j+1)*nx] - Rhof[(j-1)*nx])/*m5*/ -dt/dy * 0.5 /(Viscw[j*nx] *S[j*nx])* g*g_flag * Rhof[j*nx] *( Kappa[(j+1)*nx] -  Kappa[(j-1)*nx])/*m6*/ -dt/dy * 0.5 /(S[j*nx])* g*g_flag * Rhof[j*nx]* Kappa[j*nx]  *(1/Viscw[(j+1)*nx]-1/Viscw[(j-1)*nx])/*m7*/;

	Pw[nx-1+j*nx]= 2*dtdx /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx]  + Kappa[nx-2+j*nx]/Viscw[nx-2+j*nx] )/*m3*/*iPw[nx-2+j*nx]+ dtdy /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx] + Kappa[nx-1+(j+1)*nx]/Viscw[nx-1+(j+1)*nx])/*m2*/*iPw[nx-1+(j+1)*nx]+dtdy /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx]  + Kappa[nx-1+(j-1)*nx]/Viscw[nx-1+(j-1)*nx] )/*m4*/*iPw[nx-1+(j-1)*nx] +(1-dtdy /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx] + Kappa[nx-1+(j+1)*nx]/Viscw[nx-1+(j+1)*nx])/*m2*/-2*dtdx /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx]  + Kappa[nx-2+j*nx]/Viscw[nx-2+j*nx] )/*m3*/-dtdy /S[nx-1+j*nx] *0.5*(Kappa[nx-1+j*nx]/Viscw[nx-1+j*nx]  + Kappa[nx-1+(j-1)*nx]/Viscw[nx-1+(j-1)*nx] )/*m4*/) * iPw[nx-1+j*nx] -dt/dy * 0.5 /(Viscw[nx-1+j*nx] *S[nx-1+j*nx])* g*g_flag * Kappa[nx-1+j*nx] *(Rhof[nx-1+(j+1)*nx] - Rhof[nx-1+(j-1)*nx])/* m5*/ -dt/dy * 0.5 /(Viscw[nx-1+j*nx] *S[nx-1+j*nx])* g*g_flag * Rhof[nx-1+j*nx] *( Kappa[nx-1+(j+1)*nx] -  Kappa[nx-1+(j-1)*nx])/*m6*/ -dt/dy * 0.5  /(S[nx-1+j*nx])* g*g_flag * Rhof[nx-1+j*nx]* Kappa[nx-1+j*nx] *(1/Viscw[nx-1+(j+1)*nx]-1/Viscw[nx-1+(j-1)*nx])/*m7*/;
}
///corners


Pw[0] =  Pw0[0] +pmin ;
Pw[nx-1] = Pw0[nx-1] +pmin;
Pw[(ny-1)*nx] = Pw0[(ny-1)*nx]+pmin;
Pw[nx-1+(ny-1)*nx] = Pw0[nx-1+(ny-1)*nx] +pmin;

///Velocity boundary conditions
#pragma omp parallel for
for(int i=1; i< nx-1; i++)
{
	//bottom
	Vwx[i+(ny-1)*nx] =  0.5 *  Kappa[i+(ny-1)*nx] /Viscw[i+(ny-1)*nx]* (iPw[(i+1)+(ny-1)*nx] - iPw[(i-1)+(ny-1)*nx])/dx;
	Vwy[i+(ny-1)*nx] =  -Kappa[i+(ny-1)*nx] /Viscw[i+(ny-1)*nx]* ((iPw[i+(ny-1)*nx] - iPw[i+(ny-2)*nx])/dy -Rhof[i+(ny-1)*nx] *g_flag * g);	//first order
	//top
	Vwx[i] = -0.5 *  Kappa[i] /Viscw[i] * (iPw[i+1] - iPw[i-1])/dx;
	Vwy[i] = -Kappa[i] /Viscw[i]* ((iPw[i+nx] - iPw[i])/dy - Rhof[i+nx]*g_flag* g); 	//first order
}
#pragma omp parallel for
for(int j=1; j< ny-1 ; j++)
{
	//left
	Vwx[j*nx] = -Kappa[j*nx] /Viscw[j*nx] * (iPw[1+j*nx] - iPw[j*nx])/dx;		//first order
	Vwy[j*nx] = - Kappa[j*nx] /Viscw[j*nx] * (0.5 * (iPw[(j+1)*nx]- iPw[(j-1)*nx])/dy - Rhof[j*nx]*g_flag* g);
	//right
	Vwx[nx-1+j*nx] =  - Kappa[nx-1+j*nx] /Viscw[nx-1+j*nx] * (iPw[nx-1+j*nx] - iPw[(nx-2)+j*nx] )/dx;
	Vwy[nx-1+j*nx] =  - Kappa[nx-1+j*nx] /Viscw[nx-1+j*nx]* (0.5 *(iPw[nx-1+(j+1)*nx]- iPw[nx-1+(j-1)*nx])/dy -Rhof[nx-1+j*nx]*g_flag * g);
}


///corners
///top left
Vwx[0] =  -Kappa[0] /Viscw[0]* (iPw[1] - iPw[0])/dx;		//first order
Vwy[0] =  -Kappa[0] /Viscw[0] * ((iPw[nx] - iPw[0])/dy - Rhof[0]*g_flag* g); 	//first order
///bottom left
Vwx[(ny-1)*nx] = -Kappa[(ny-1)*nx] /Viscw[(ny-1)*nx] * (iPw[1+(ny-1)*nx] - iPw[(ny-1)*nx])/dx;
Vwy[(ny-1)*nx] = -Kappa[(ny-1)*nx] /Viscw[(ny-1)*nx] * ((iPw[(ny-1)*nx] - iPw[(ny-2)*nx])/dy - Rhof[(ny-1)*nx]*g_flag* g);
///top right
Vwx[nx-1] =  -Kappa[nx-1] /Viscw[nx-1]* (iPw[nx-1] - iPw[(nx-2)])/dx;
Vwy[nx-1] =  -Kappa[nx-1] /Viscw[nx-1] * ((iPw[nx-1+nx] - iPw[nx-1])/dy - Rhof[nx-1]*g_flag* g);
///bottom right
Vwx[nx-1+(ny-1)*nx] =  -Kappa[nx-1+(ny-1)*nx] /Viscw[nx-1+(ny-1)*nx] * (iPw[nx-1+(ny-1)*nx] - iPw[(nx-2)+(ny-1)*nx] )/dx;
Vwy[nx-1+(ny-1)*nx] =  -Kappa[nx-1+(ny-1)*nx] /Viscw[nx-1+(ny-1)*nx] * ((iPw[nx-1+(ny-1)*nx] - iPw[nx-1+(ny-2)*nx])/dy - Rhof[nx-1+(ny-1)*nx]*g_flag* g);

///injection

Pw[int(nx/2)+int(ny/2)*nx]=Pw0[int(nx/2)+int(ny/2)*nx] + pmax ;
Vwx[int(nx/2)+int(ny/2)*nx] = -Kappa[int(nx/2)+int(ny/2)*nx] /Viscw[int(nx/2)+int(ny/2)*nx] *(Pw[int(nx/2)+1+int(ny/2)*nx]-Pw[int(nx/2)+int(ny/2)*nx])/dx;
Vwy[int(nx/2)+int(ny/2)*nx] = -Kappa[int(nx/2)+int(ny/2)*nx] /Viscw[int(nx/2)+int(ny/2)*nx] *(Pw[int(nx/2)+int(ny/2)*nx+nx]-Pw[int(nx/2)+int(ny/2)*nx])/dy;

}

template <typename T>
void TimeStepCalculation(T &dt,T *Phi, T *Kappa, T *Viscw, T *Betaphi, T *Betaw, T *Rhow,T *TCm,T *Dm,T *Cpw,T tcr,T tcw,T cpr,T rhor0,T dx,T dy,int nx, int ny)
{
	T kappapf,kappac,kappat,kappatr,kappatf;
	T maxkappa,maxkappa1,maxkappa2,maxkappa3,maxkappa4,maxkappa5;
	maxkappa=-10e10;
	maxkappa1=maxkappa2=maxkappa3=maxkappa4=maxkappa5=maxkappa;

	///determining dt
	for (int j=0; j<ny; j++)
	{
		for (int i=0; i<nx; i++)
		{
			//TCm[i+j*nx] = (1-Phi[i+j*nx])*tcr + Phi[i+j*nx]*tcw;
			kappapf = Kappa[i+j*nx]/(Viscw[i+j*nx]*Phi[i+j*nx]*(Betaw[i+j*nx]+Betaphi[i+j*nx]));					//pressure diffusivity
			//kappat = TCm[i+j*nx]/(Cpw[i+j*nx]  * Phi[i+j*nx] * Rhow[i+j*nx] + cpr * (1-Phi[i+j*nx]) * rhor0);			//thermal diffusivity
			kappac = Dm[i+j*nx]/Phi[i+j*nx];															//salinity diffusivity
			//kappatr = tcr/(cpr * rhor0);
			//kappatf = tcw /(Cpw[i+j*nx]  * Rhow[i+j*nx]);

			maxkappa1=max(maxkappa1,kappapf);
			//maxkappa2=max(maxkappa2,kappat);
			maxkappa3=max(maxkappa3,kappac);
			//maxkappa4=max(maxkappa4,kappatr);
			//maxkappa5=max(maxkappa5,kappatf);
		}
	}

	//maxkappa=max(maxkappa1,max(maxkappa2,max(maxkappa3,max(maxkappa4,maxkappa5))));

	maxkappa=max(maxkappa1,maxkappa3);

	dt = 0.1*min(dx,dy)*min(dx,dy)/maxkappa;
}

template <typename T>
void AdvectionGamma( T *Gamma, T *Alpha, T *Rhof ,T *Cpw, T *Phi,T rhor0, T cpr,T dx, T dy, int nx, int ny)
{
	for (int j=0; j<ny; j++)
	{
		for (int i=0; i<nx; i++)
		{
			Alpha [i+j*nx]= (rhor0*cpr)/(Rhof[i+j*nx]*Cpw[i+j*nx]);
			Gamma[i+j*nx] = Phi[i+j*nx]/(Alpha[i+j*nx]*(1-Phi[i+j*nx])+Phi[i+j*nx]);
		}
	}
}

template <typename T>
void AdvectionCoefficient( T *Adv, T *Rhof ,T *Cpw, T *Phi,T rhor0, T cpr, int lte_flag,int nx, int ny)
{
	if(lte_flag ==1)
	{
#pragma omp parallel for
		for (int j=0; j<ny; j++)
		{
			for (int i=0; i<nx; i++)
			{
				Adv[i+j*nx]= Rhof[i+j*nx]*Cpw[i+j*nx] /((1-Phi[i+j*nx])*rhor0*cpr + Phi[i+j*nx]*Rhof[i+j*nx]*Cpw[i+j*nx] );
			}
		}
	}
	else
	{
#pragma omp parallel for
		for (int j=0; j<ny; j++)
		{
			for (int i=0; i<nx; i++)
			{
				Adv[i+j*nx]= 1 / Phi[i+j*nx];
			}
		}

	}

}


template <typename T>
void Fracture(T *K, T *Phi,T *TCr, T *Pw, T *Szz, T *Sxx, T k_fracture, T perm_max,T phi_fracture,T depth_fracture, T width_fracture,T theta, T fracture_aperture, T l_fracture, T tcr_fracture,T phiini, T phir, T phia, T c, double PI, T dx, T dy, int nx, int ny)
{

	int fx1=0, fy1=0,fy2=0,fx2=0, fx_min, fy_min, fx_max, fy_max;
	double x1, x2, y1, y2;
	double m=0;			//the slope
	int lx=0, ly=0;					//x and y components of the fracture length
	double Mx=0, My=0;					//middle of the cell: x0 & y0
	double Dx, Dy,  d=0, r=0;							//distance
	theta*=PI/180;
	m=tan(theta);
	//initial point
	x1=width_fracture;
	y1=depth_fracture;
	fx1=int(width_fracture/dx);
	fy1=int(depth_fracture/dy);
	//length
	Dx=l_fracture*cos(theta);
	Dy=l_fracture*sin(theta);
	lx=int(Dx/dx);
	ly=int(Dy/dy);
	//final point
	x2=x1+Dx;
	y2=y1+Dy;
	fx2=lx+fx1;
	fy2=ly+fy1;
	fx_min=min(fx1,fx2);
	fx_max=max(fx1,fx2);
	fy_min=min(fy1,fy2);
	fy_max=max(fy1,fy2);

	//arbitrary distance: 1/2*diagonal
	r=0.5*sqrt(dx*dx+dy*dy);

	//cout<<"fx1:  "<<x1<<"  fy1:  "<<y1<<"  fx2:  "<<x2 <<"  fy2:  "<<y2 <<"  l:  "<<l_fracture <<"  slope:  "<<m<<"  aperture:  "<<fracture_aperture<<endl;
	for(int j=fy_min-2; j<=fy_max+2; j++){
		for(int i=fx_min-2; i<=fx_max+2; i++){
			Mx=(i+0.5)*dx;
			My=(j+0.5)*dy;
			d=abs(Dy*Mx-Dx*My-x1*y2+x2*y1)/l_fracture;
			if (d<=fracture_aperture){
				Phi[i+j*nx]*=3;//phir+(phi_fracture- phir) * exp(phia * ( 0.5*(Sxx[i+j*nx]+Szz[i+j*nx]) - Pw[i+j*nx]));
				if(Phi[i+j*nx]>=1){
					Phi[i+j*nx]=0.95;
				}
				K[i+j*nx] *= exp(c * (Phi[i+j*nx]/phiini - 1.0));
				if(K[i+j*nx]>perm_max){
					K[i+j*nx]=perm_max;
				}
			}

		}

	}


}
template <typename T>
void SaveDistAve(T *M0,T *M1, T *M2,T *M3,T *M4,T *M5, T *M6, T *M7, T *M8, string name, T dx,T dy,  int tim, int nx,int ny)
{
	char day[15];
	gcvt(tim,-1,day);
	string  txt;
	txt = ".txt";
	name.append(day);
	name.append(txt);
	double m0=0, m1 =0, m2 =0, m3 =0, m4 =0, m5 = 0, m6=0, m7=0, m8=0;
	const char *namenum;
	namenum=name.c_str();
	ofstream a_file (namenum);
	a_file<<"#Average value in distance"<<endl;
	a_file<<"#d\t \t P_over\t\t Pw(d)\t\t  Tw(d)\t\t Tr(d) \t \t Cs(d) \t\t Rhow(d) \t\t Phi(d) \t\t Cpf(d) \t\t AdvectionCoefficient  "<<endl;
	int x0 = int (nx/2);
	int y0 = int (ny/2);
	for (int  d = 0;d <x0; d++){
		m0=0.25*(M1[x0+d+(y0+d)*nx]+M1[x0-d+(y0+d)*nx]+M1[x0-d+(y0-d)*nx]+M1[x0+d+(y0-d)*nx]-(M0[x0+d+(y0+d)*nx]+M0[x0-d+(y0+d)*nx]+M0[x0-d+(y0-d)*nx]+M0[x0+d+(y0-d)*nx]));
		m1=0.25*(M1[x0+d+(y0+d)*nx]+M1[x0-d+(y0+d)*nx]+M1[x0-d+(y0-d)*nx]+M1[x0+d+(y0-d)*nx]);
		m2=0.25*(M2[x0+d+(y0+d)*nx]+M2[x0-d+(y0+d)*nx]+M2[x0-d+(y0-d)*nx]+M2[x0+d+(y0-d)*nx]);
		m3=0.25*(M3[x0+d+(y0+d)*nx]+M3[x0-d+(y0+d)*nx]+M3[x0-d+(y0-d)*nx]+M3[x0+d+(y0-d)*nx]);
		m4=0.25*(M4[x0+d+(y0+d)*nx]+M4[x0-d+(y0+d)*nx]+M4[x0-d+(y0-d)*nx]+M4[x0+d+(y0-d)*nx]);
		m5=0.25*(M5[x0+d+(y0+d)*nx]+M5[x0-d+(y0+d)*nx]+M5[x0-d+(y0-d)*nx]+M5[x0+d+(y0-d)*nx]);
		m6=0.25*(M6[x0+d+(y0+d)*nx]+M6[x0-d+(y0+d)*nx]+M6[x0-d+(y0-d)*nx]+M6[x0+d+(y0-d)*nx]);
		m7=0.25*(M7[x0+d+(y0+d)*nx]+M7[x0-d+(y0+d)*nx]+M7[x0-d+(y0-d)*nx]+M7[x0+d+(y0-d)*nx]);
		m8=0.25*(M8[x0+d+(y0+d)*nx]+M8[x0-d+(y0+d)*nx]+M8[x0-d+(y0-d)*nx]+M8[x0+d+(y0-d)*nx]);

		a_file<< d*sqrt(dx*dx+dy*dy) <<"\t\t"<<m0<<"\t\t"<<m1<<"\t\t"<<m2<<"\t\t"<<m3<<"\t\t"<<m4<<"\t\t"<<m5<<"\t\t"<<m6<<"\t\t"<<m7<<"\t\t"<<m8<<endl;
	}
	a_file.close();
}
template <typename T>
void SaveVector(T *Mx,  T *My, string name, T dx, T dy, int tim, int nx, int ny)
{
	char day[15];
	gcvt(tim,-1,day);
	string  txt;
	txt = ".txt";
	name.append(day);
	name.append(txt);

	const char *namenum;
	namenum=name.c_str();
	ofstream a_file (namenum);
	a_file<<"#Velocity field in 2D"<<endl;
	a_file<<"#X\t\t y\t\t Vx\t \tVy\t "<<endl;

	for (int  j = 0; j <ny; j+=3){
		for (int i = 0; i <nx ; i+=3){
			a_file<< i*dx <<"\t\t "<< j*dy<<"\t\t"<< Mx[i+j*nx]*2e3  <<"\t\t"<< My[i+j*nx]*2e3 <<endl;
		}
	}
	a_file.close();
}

template <typename T>
void SaveVtk(T *Pw,T *Cs,T *Tf, T *Tr, T *Kappa, T *Phi,T *Vx,  T *Vy, string name, T dx, T dy, int tim, int nx, int ny)
{
	char day[15];
	gcvt(tim,-1,day);
	string  txt;
	txt = ".vtk";
	name.append(day);
	name.append(txt);

	const char *namenum;
	namenum=name.c_str();
	ofstream a_file (namenum);
	a_file<< "# vtk DataFile Version 5.10"<<endl;
	a_file<<"Fluid Pressure, Temperature"<<endl;
	a_file<<"ASCII"<<endl;
	a_file<<"DATASET STRUCTURED_POINTS"<<endl;
	a_file<<"DIMENSIONS "<<nx<<" "<<ny<<" "<<"1 "<<endl;
	a_file<<"ORIGIN 0 0 0"<<endl;
	a_file<<"SPACING "<<dx<<" "<<dy<<" "<<"0 "<<endl;
	a_file<<"POINT_DATA "<<nx*ny<<endl;

	a_file<<"SCALARS Pressure double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Pw[i+j*nx]<<endl;
		}
	}

	a_file<<"SCALARS concentration double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Cs[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Tf double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Tf[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Tr double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Tr[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Kappa double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Kappa[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Rho_f double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Phi[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Vx double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Vx[i+j*nx] <<endl;
		}
	}

	a_file<<"SCALARS Vy double "<<endl;
	a_file<<"LOOKUP_TABLE default"<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Vy[i+j*nx] <<endl;
		}
	}


	a_file<<"VECTORS V double "<<endl;

	for (int  j = ny-1; j >=0; j--){
		for (int i = 0; i <nx ; i++){
			a_file<< Vx[i+j*nx] <<"\t"<<Vy[i+j*nx]<<"\t"<<"0.0"<<endl;
		}
	}
	a_file.close();
}


template <typename T>
void DensityTransportSML(T *Cs,T *iCs,T *Dm, T *Phi,T *Vwx,T *Vwy,T *x,T *y,T *Rhof,T h,T e, T c_fresh,T dx,T dy, T t, T teq,T dt,int nx,int ny)
{
	char pause;
	T rhocp, m1, m2, m3, m4;
	T dtdy = dt/(dy*dy);
	T dtdx = dt/(dx*dx);
	T val;
	createHostVariable(T,tCs, nx*ny);
	//         	for(int j=int(ny/2)-1; j<= int(ny/2)+1 ; j++)						//bore hole   Dirichlet
	//         	{
	//         		for(int i=int(nx/2)-1; i<=int(nx/2)+1; i++)
	//         		{
	//         			iCs[i+j*nx]=c_fresh;
	//
	//         		}
	//         	}
#pragma omp parallel for
	for (int j=0; j<ny; j++)
	{
		for (int i=0; i<nx; i++)
		{
			tCs[i+j*nx] = iCs[i+j*nx];
		}
	}

#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{
		for(int i=1; i<nx-1; i++)
		{
			// 			m1 = dtdx/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx] + Phi[(i+1)+j*nx]* Rhof[(i+1)+j*nx]);
			// 			m2 = dtdy/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+ Phi[i+(j+1)*nx]* Rhof[i+(j+1)*nx]);
			// 			m3 = dtdx/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+Phi[(i-1)+j*nx]* Rhof[(i-1)+j*nx]);
			// 			m4 = dtdy/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+ Phi[i+(j-1)*nx]* Rhof[i+(j-1)*nx]);

			tCs[i+j*nx] = iCs[i+j*nx] + dtdx/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx] + Phi[(i+1)+j*nx]* Rhof[(i+1)+j*nx])/*m1*/*(iCs[(i+1)+j*nx]-iCs[i+j*nx])+dtdx/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+Phi[(i-1)+j*nx]* Rhof[(i-1)+j*nx])/*m3*/*(iCs[(i-1)+j*nx]-iCs[i+j*nx])+dtdy/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+ Phi[i+(j+1)*nx]* Rhof[i+(j+1)*nx])/*m2*/*(iCs[i+(j+1)*nx]-iCs[i+j*nx])+dtdy/Phi[i+j*nx] /Rhof[i+j*nx] * 0.5 * Dm[i+j*nx]*(Phi[i+j*nx] * Rhof[i+j*nx]+ Phi[i+(j-1)*nx]* Rhof[i+(j-1)*nx])/*m4*/*(iCs[i+(j-1)*nx]-iCs[i+j*nx]);

		}
	}
#pragma omp parallel for
	for(int j=1 ; j<ny-1; j++)
	{
		for(int i=1; i<nx-1; i++)
		{
			if(tCs[i+j*nx]>0.021)
			{
				cout<<"Cs overdosed="<<tCs[i+j*nx]<<"  "<<i<<"  "<<j<<"  t:  "<<int(t/dt)<<"  "<<endl;//<<"    iTr-tTw:   "<< iTr[i+j*nx]-tTw[i+j*nx]<<"  "<<"   iTr- iTw:   "<<iTr[i+j*nx]-iTw[i+j*nx] <<" Qdtw:   "<< Qdtw<<"    m:   "<<m1<<endl;
				pause=getchar();
				tCs[i+j*nx]= 0.02;

			}
		}
	}

	// boundary conditions for bottom j=ny-1 ; top j=0
	// Elders problem
    #pragma omp parallel for
   	for(int i=1; i< nx-1; i++)
    {
        		tCs[i+(ny-1)*nx] = iCs[i+(ny-1)*nx] ;
        		tCs[i] = iCs[i] ;
     }
	/*        	for(int i=1; i< nx-1; i++)
        	{
        		m1 = dtdx * 0.5 /Phi[i+(ny-1)*nx]* (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx] + Dm[(i+1)+(ny-1)*nx]*Phi[(i+1)+(ny-1)*nx]);
        		m3 = dtdx * 0.5 /Phi[i+(ny-1)*nx]* (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx] + Dm[(i-1)+(ny-1)*nx]*Phi[(i-1)+(ny-1)*nx]);
        		m4 = dtdy * 0.5 /Phi[i+(ny-1)*nx]* (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx] + Dm[i+((ny-1)-1)*nx]*Phi[i+((ny-1)-1)*nx]);
        		tCs[i+(ny-1)*nx] = iCs[i+(ny-1)*nx] + m1*iCs[(i+1)+(ny-1)*nx] + m3*iCs[(i-1)+(ny-1)*nx] + 2*m4*iCs[i+(ny-2)*nx]-(m1+m3+2*m4)*iCs[i+(ny-1)*nx];

        		m1 = dtdx * 0.5 /Phi[i]* (Dm[i]*Phi[i] + Dm[(i+1)]*Phi[(i+1)]);
        		m2 = dtdy * 0.5 /Phi[i]* (Dm[i]*Phi[i] + Dm[i+nx]*Phi[i+nx]);
        		m3 = dtdx * 0.5 /Phi[i]* (Dm[i]*Phi[i] + Dm[(i-1)]*Phi[(i-1)]);
        		tCs[i] = iCs[i] + m1*iCs[(i+1)]+m3*iCs[(i-1)]+2*m2*iCs[i+nx]-(m1+2*m2+m3)*iCs[i];
        	}*/

	// boundary conditions for left i=0; right i=nx-1

	for(int j=1; j< ny-1; j++)
	{
		m1 = dtdx * 0.5 /Phi[j*nx]* (Dm[j*nx]*Phi[j*nx]  + Dm[1+j*nx]*Phi[1+j*nx]);
		m2 = dtdy * 0.5 /Phi[j*nx]* (Dm[j*nx]*Phi[j*nx] + Dm[(j+1)*nx]*Phi[(j+1)*nx]);
		m4 = dtdy * 0.5 /Phi[j*nx]* (Dm[j*nx]*Phi[j*nx] + Dm[(j-1)*nx]*Phi[(j-1)*nx]);
		tCs[j*nx] = iCs[j*nx] + m1*iCs[(1)+j*nx]+m2*iCs[(j+1)*nx]+2*m4*iCs[(j-1)*nx]-(m1+m2+2*m4)*iCs[j*nx];

		m2 = dtdy * 0.5 /Phi[(nx-1)+j*nx]* (Dm[(nx-1)+j*nx]*Phi[(nx-1)+j*nx] + Dm[(nx-1)+(j+1)*nx]*Phi[(nx-1)+(j+1)*nx]);
		m3 = dtdx * 0.5 /Phi[(nx-1)+j*nx]* (Dm[(nx-1)+j*nx]*Phi[(nx-1)+j*nx] + Dm[((nx-1)-1)+j*nx]*Phi[((nx-1)-1)+j*nx]);
		m4 = dtdy * 0.5 /Phi[(nx-1)+j*nx]* (Dm[(nx-1)+j*nx]*Phi[(nx-1)+j*nx] + Dm[(nx-1)+(j-1)*nx]*Phi[(nx-1)+(j-1)*nx]);
		tCs[(nx-1)+j*nx] = iCs[(nx-1)+j*nx]+2*m3*iCs[(nx-2)+j*nx]+m2*iCs[(nx-1)+(j+1)*nx]+m4*iCs[(nx-1)+(j-1)*nx]-(m2+2*m3+m4)*iCs[(nx-1)+j*nx];
	}
	//bore hole Dirichlet
	//         	for(int j=int(ny/2)-1; j<= int(ny/2)+1 ; j++)						//bore hole   Dirichlet
	//         	{
	//         		for(int i=int(nx/2)-1; i<=int(nx/2)+1; i++)
	//         		{
	//         			tCs[i+j*nx]=c_fresh;
	//         		}
	//         	}
	// 		 for(int j=int(h/dy); j<=int((e+h)/dy) ; j++)
	// 		{
	// 			tCs[j*nx]=0.0;
	// 		}
	//boundary conditions for corners
	//left bottom i=0, j =ny-1
	m1 = dtdx * 0.5 /Phi[(ny-1)*nx]* (Dm[(ny-1)*nx]*Phi[(ny-1)*nx] + Dm[1+(ny-1)*nx]*Phi[1+(ny-1)*nx]);
	m4 = dtdy * 0.5 /Phi[(ny-1)*nx]* (Dm[(ny-1)*nx]*Phi[(ny-1)*nx] + Dm[((ny-1)-1)*nx]*Phi[((ny-1)-1)*nx]);
	tCs[(ny-1)*nx] = iCs[(ny-1)*nx] + m1*iCs[(1)+(ny-1)*nx]+m4*iCs[(ny-2)*nx]-(2*m1+2*m4)*iCs[(ny-1)*nx];
	//right bottom i=nx-1,j =ny-1
	m3 = dtdx * 0.5 /Phi[(nx-1)+(ny-1)*nx]* (Dm[(nx-1)+(ny-1)*nx]*Phi[(nx-1)+(ny-1)*nx] + Dm[((nx-1)-1)+(ny-1)*nx]*Phi[((nx-1)-1)+(ny-1)*nx]);
	m4 = dtdy * 0.5 /Phi[(nx-1)+(ny-1)*nx]* (Dm[(nx-1)+(ny-1)*nx]*Phi[(nx-1)+(ny-1)*nx] + Dm[(nx-1)+((ny-1)-1)*nx]*Phi[(nx-1)+((ny-1)-1)*nx]);
	tCs[(nx-1)+(ny-1)*nx] = iCs[(nx-1)+(ny-1)*nx] + 2*m3*iCs[(nx-2)+(ny-1)*nx] + 2*m4*iCs[(nx-1)+(ny-2)*nx]-(2*m3+2*m4)*iCs[(nx-1)+(ny-1)*nx];
	// left top i=0,j=0 m1, m3
	m1 = dtdx * 0.5 /Phi[0]* (Dm[0]*Phi[0] + Dm[1]*Phi[1]);
	m2 = dtdy * 0.5 /Phi[0]* (Dm[0]*Phi[0] + Dm[nx]*Phi[nx]);
	tCs[0] = iCs[0] + 2*m1*iCs[1]+2*m2*iCs[nx]-(2*m1+2*m2)*iCs[0];
	// right top i=nx-1,j=0 m2, m3
	m2 = dtdy * 0.5 /Phi[(nx-1)]* (Dm[(nx-1)]*Phi[(nx-1)] + Dm[(nx-1)+nx]*Phi[(nx-1)+nx]);
	m3 = dtdx * 0.5 /Phi[(nx-1)]* (Dm[(nx-1)]*Phi[(nx-1)] + Dm[((nx-1)-1)]*Phi[((nx-1)-1)]);
	tCs[(nx-1)] = iCs[(nx-1)] + 2*m3*iCs[(nx-2)] + 2*m2*iCs[(nx-1)+nx]-(2*m2+2*m3)*iCs[(nx-1)];

	T npox, npoy;
	int yind[4];
	int xind[4];
	T ivalx[4];
	T ivaly[4];
	T ival,valx,valy;
	T alphal;
#pragma omp parallel private(npox,npoy,yind,xind,ivalx,ivaly,ival,valx,valy,alphal)
	{
#pragma omp for
		for(int j=1; j<ny-1; j++)
		{
			for(int i=1; i<nx-1; i++)
			{
				valx=Vwx[i+j*nx]/Phi[i+j*nx];
				valy=Vwy[i+j*nx]/Phi[i+j*nx];

				//advection
				if(fabs(valx)>1.0e-10 || fabs(valy)>1.0e-10)
				{
					//get new positions
					npox=x[i+j*nx]-dt*valx;
					npoy=y[i+j*nx]-dt*valy;

					//interpolation I need 16 points
					InterpolationPoints<T>(xind,yind,npox,npoy,x,y,nx,ny);
					for(int jj =0; jj<4;jj++){

						for(int ii =0; ii<4;ii++)
						{
							ivalx[ii]=tCs[xind[ii]+yind[jj]*nx];
						}
						ivaly[jj]=spline_cubic(ivalx,(npox-x[xind[1]+yind[1]*nx])/dx);
					}
					ival=spline_cubic(ivaly,(npoy-y[xind[1]+yind[1]*nx])/dy);
					Cs[i+j*nx]=ival;

				}
			}
		}
	}

	// Elders problem
	#pragma omp parallel for
	for(int i=0; i< nx; i++)
	{
			Cs[i+(ny-1)*nx] = iCs[i+(ny-1)*nx];
			Cs[i] = iCs[i];

	}
	/*
        #pragma omp parallel for
        	for(int i=0; i< nx; i++)
        	{
        			Cs[i+(ny-1)*nx] = Cs[i+(ny-2)*nx];
        			Cs[i] = Cs[i+nx];

        	}
	 */

	// boundary conditions for left i=0; right i=nx-1
#pragma omp parallel for
	for(int j=1; j< ny-1; j++)
	{
		Cs[j*nx] = Cs[1+j*nx];
		Cs[(nx-1)+j*nx] = Cs[(nx-2)+j*nx];
	}

	// 		for(int j=int(h/dy); j<=int((e+h)/dy) ; j++)
	// 		{			//bore hole Dirichlet
	// 			Cs[j*nx]=0.0;
	// 		}
	// #pragma omp parallel for
	// 	for(int j=int(ny/2)-1; j<= int(ny/2)+1 ; j++)						//bore hole   Dirichlet
	// 	{
	// 		for(int i=int(nx/2)-1; i<=int(nx/2)+1; i++)
	// 		{
	// 			Cs[i+j*nx]=c_fresh;
	// 		}
	// 	}
	Cs[int(nx/2)+int(ny/2)*nx]=c_fresh;

	DeleteHostVariable(tCs);

}


template <typename T>
void InterpolationPoints(int *xind,int *yind,T &npox,T &npoy,T *x, T *y,int nx,int ny)
{

	T xmin2=x[0];
	T xmax2=x[nx-1];
	T ymin2=y[0];
	T ymax2=y[(ny-1)*nx];

	T xmax=max(xmin2,xmax2);
	T xmin=min(xmin2,xmax2);
	T ymax=max(ymin2,ymax2);
	T ymin=min(ymin2,ymax2);

	T dx=fabs(x[1]-x[0]);
	T dy=fabs(y[nx]-y[0]);

	int indy1,indx1;

	if(npox>=xmax){npox=xmax-0.5*dx;xind[0]=nx-3;xind[1]=nx-2;xind[2]=nx-1;xind[3]=nx-1;}
	if(npox<=xmin){npox=xmin;xind[0]=0;xind[1]=0;xind[2]=1;xind[3]=2;}
	if(npox<xmax && npox>xmin){

		indx1=(int)((npox-xmin)/dx);


		xind[0]=indx1-1;xind[1]=indx1;xind[2]=indx1+1;xind[3]=indx1+2;

		if(indx1==0){
			xind[0]=indx1;xind[1]=indx1;xind[2]=indx1+1;xind[3]=indx1+2;
		}
		if(indx1==nx-2){
			xind[0]=indx1-1;xind[1]=indx1;xind[2]=indx1+1;xind[3]=indx1+1;
		}

	}

	if(npoy>=ymax){npoy=ymax-0.5*dy;yind[0]=0;yind[1]=0;yind[2]=1;yind[3]=2;}
	if(npoy<=ymin){npoy=ymin;yind[0]=ny-3;yind[1]=ny-2;yind[2]=ny-1;yind[3]=ny-1;}


	if(npoy<ymax && npoy>ymin){

		indy1=(int)((npoy-ymin)/dy);
		yind[0]=indy1-1;yind[1]=indy1;yind[2]=indy1+1;yind[3]=indy1+2;

		if(indy1==0){
			yind[0]=indy1;yind[1]=indy1;yind[2]=indy1+1;yind[3]=indy1+2;
		}
		if(indy1==ny-2){
			yind[0]=indy1-1;yind[1]=indy1;yind[2]=indy1+1;yind[3]=indy1+1;
		}
	}
}

template <typename T>
T spline_cubic(const T a[4],T x)
{
	if( a[0] == a[1] && a[0] == a[2] && a[0] == a[3] ) return a[0];

	T alpha1=3. * (a[2] - a[1]) - 3. * (a[1] - a[0]);
	T alpha2=3. * (a[3] - a[2]) - 3. * (a[2] - a[1]);

	T l2,z1,z2;

	z1=0.25*alpha1;
	l2= 4. - 0.25;
	z2=(alpha2-z1)/l2;

	T c1,b1,d1;

	c1= z1 - 0.25 * z2;
	b1 = a[2] - a[1] - (z2 + 2. * c1) / 3.;
	d1 = (z2 - c1) / 3.;

	return set2limits( a[1] + b1 * x + c1 * x * x + d1 * x * x * x, fmin(a[1],a[2]), fmax(a[1],a[2]) );
}

template <typename T>
void EGSMassFraction(T *Cs, T *iCs, T *Dm, T *Vwx, T *Vwy, T *Phi, T rhow0, T rho_max,T h,T e,T dx, T dy, T dt, int nx, int ny)
{

	T m1, m2, m3, m4, m5, m6, m7, m8;
	T vxp,vxm, vyp,vym;

	T dtdy = dt/(dy*dy);
	T dtdx = dt/(dx*dx);


#pragma omp parallel for									//Masstransport eq.
for(int j=1; j< ny-1; j++)
{
	for(int i=1; i< nx-1; i++)
	{
		m1 = dtdx * 0.5 /Phi[i+j*nx]  * (Dm[i+j*nx]*Phi[i+j*nx] + Dm[(i+1)+j*nx]*Phi[(i+1)+j*nx]);
		m2 = dtdy * 0.5 /Phi[i+j*nx]  * (Dm[i+j*nx]*Phi[i+j*nx] + Dm[i+(j+1)*nx]*Phi[i+(j+1)*nx]);
		m3 = dtdx * 0.5 /Phi[i+j*nx]  * (Dm[i+j*nx]*Phi[i+j*nx] + Dm[(i-1)+j*nx]*Phi[(i-1)+j*nx]);
		m4 = dtdy * 0.5 /Phi[i+j*nx]  * (Dm[i+j*nx]*Phi[i+j*nx] + Dm[i+(j-1)*nx]*Phi[i+(j-1)*nx]);
		vxp = dt/dx * max(Vwx[i+j*nx],T(0)) ;
		vxm = dt/dx * min(Vwx[i+j*nx],T(0)) ;
		vyp = dt/dy * max(Vwy[i+j*nx],T(0)) ;
		vym = dt/dy * min(Vwy[i+j*nx],T(0)) ;

		Cs[i+j*nx] = m1*iCs[(i+1)+j*nx]+ m3*iCs[(i-1)+j*nx]+ m2*iCs[i+(j+1)*nx]+m4*iCs[i+(j-1)*nx] +(1-m1-m2-m3-m4) * iCs[i+j*nx]
																														  -vxp*iCs[i+j*nx]+vxp*iCs[(i-1)+j*nx]-vxm*iCs[(i+1)+j*nx]+vxm*iCs[i+j*nx]-vyp*iCs[i+j*nx]+vyp*iCs[i+(j-1)*nx]-vym*iCs[i+(j+1)*nx]+vym*iCs[i+j*nx];

	}
}


#pragma omp parallel for											//boundary condition EGS problem
for(int i=1; i< nx-1; i++)
{
	//bottom
	m1 = dtdx * 0.5 /Phi[i+(ny-1)*nx] * (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx]+Dm[(i+1)+(ny-1)*nx]*Phi[(i+1)+(ny-1)*nx]);
	m3 = dtdx * 0.5 /Phi[i+(ny-1)*nx] * (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx]+Dm[(i-1)+(ny-1)*nx]*Phi[(i-1)+(ny-1)*nx]);
	m4 = dtdy * 0.5 /Phi[i+(ny-1)*nx] * (Dm[i+(ny-1)*nx]*Phi[i+(ny-1)*nx]+Dm[i+(ny-2)*nx]*Phi[i+(ny-2)*nx]);
	vxp = dt/dx * max(Vwx[i+(ny-1)*nx],T(0)) ;
	vxm = dt/dx * min(Vwx[i+(ny-1)*nx],T(0)) ;
	vyp = dt/dy * max(Vwy[i+(ny-1)*nx],T(0)) ;

	Cs[i+(ny-1)*nx]= m1*iCs[(i+1)+(ny-1)*nx]+ m3*iCs[(i-1)+(ny-1)*nx] +2*m4*Cs[i+(ny-2)*nx] +(1-m1-m3-2*m4) * iCs[i+(ny-1)*nx]
																												  -vxp*iCs[i+(ny-1)*nx]+vxp*iCs[(i-1)+(ny-1)*nx]-vxm*iCs[(i+1)+(ny-1)*nx]+vxm*iCs[i+(ny-1)*nx]-vyp*iCs[i+(ny-1)*nx]+vyp*iCs[i+(ny-2)*nx];			//vym=0
																												  //top
	m1 = dtdx * 0.5 /Phi[i]  * (Dm[i]*Phi[i] + Dm[i+1]*Phi[i+1]);
	m2 = dtdy * 0.5 /Phi[i]  * (Dm[i]*Phi[i] + Dm[i+nx]*Phi[i+nx]);
	m3 = dtdx * 0.5 /Phi[i]  * (Dm[i]*Phi[i] + Dm[i-1]*Phi[i-1]);
	vxp = dt/dx * max(Vwx[i],T(0)) ;
	vxm = dt/dx * min(Vwx[i],T(0)) ;
	vym = dt/dy * min(Vwy[i],T(0)) ;

	Cs[i]=m1*iCs[i+1]+ m3*iCs[i-1]+ 2*m2*iCs[i+nx] +(1-m1-m3-2*m2) * iCs[i]-vxp*iCs[i]+vxp*iCs[i-1]-vxm*iCs[i+1]+vxm*iCs[i]-vym*iCs[i+nx]+vym*iCs[i]; 	//vyp=0
}
#pragma omp parallel for
for(int j=1; j< ny-1 ; j++)
{
	//left i=0
			m1 = dtdx * 0.5 /Phi[j*nx]  * (Dm[j*nx]*Phi[j*nx] + Dm[1+j*nx]*Phi[1+j*nx]);
			m2 = dtdy * 0.5 /Phi[j*nx]  * (Dm[j*nx]*Phi[j*nx] + Dm[(j+1)*nx]*Phi[(j+1)*nx]);
			m4 = dtdy * 0.5 /Phi[j*nx]  * (Dm[j*nx]*Phi[j*nx] + Dm[(j-1)*nx]*Phi[(j-1)*nx]);
			vxm = dt/dx * min(Vwx[j*nx],T(0)) ;
			vyp = dt/dy * max(Vwy[j*nx],T(0)) ;
			vym = dt/dy * min(Vwy[j*nx],T(0)) ;
			Cs[j*nx]= 2*m1*iCs[1+j*nx]+ m2*iCs[(j+1)*nx]+m4*iCs[(j-1)*nx] +(1-2*m1-m2-m4) * iCs[j*nx]
																								-vxm*iCs[1+j*nx]+vxm*iCs[j*nx]-vyp*iCs[j*nx]+vyp*iCs[(j-1)*nx]-vym*iCs[(j+1)*nx]+vym*iCs[j*nx];				//vxp=0
																								//		if(j>ny-6){						//fresh water
			//		Cs[j*nx]=0;
			//		}

			//right
			m2 = dtdy * 0.5 /Phi[nx-1+j*nx]  * (Dm[nx-1+j*nx]*Phi[nx-1+j*nx] + Dm[nx-1+(j+1)*nx]*Phi[nx-1+(j+1)*nx]);
			m3 = dtdx * 0.5 /Phi[nx-1+j*nx]  * (Dm[nx-1+j*nx]*Phi[nx-1+j*nx] + Dm[(nx-2)+j*nx]*Phi[(nx-2)+j*nx]);
			m4 = dtdy * 0.5 /Phi[nx-1+j*nx]  * (Dm[nx-1+j*nx]*Phi[nx-1+j*nx] + Dm[nx-1+(j-1)*nx]*Phi[nx-1+(j-1)*nx]);
			vxp = dt/dx * max(Vwx[nx-1+j*nx],T(0)) ;
			vyp = dt/dy * max(Vwy[nx-1+j*nx],T(0)) ;
			vym = dt/dy * min(Vwy[nx-1+j*nx],T(0)) ;
			Cs[nx-1+j*nx]= 2*m3*iCs[(nx-2)+j*nx]+ m2*iCs[nx-1+(j+1)*nx]+m4*iCs[nx-1+(j-1)*nx] +(1-m2-2*m3-m4) * iCs[nx-1+j*nx]
																													-vxp*iCs[nx-1+j*nx]+vxp*iCs[(nx-2)+j*nx]-vyp*iCs[nx-1+j*nx]+vyp*iCs[nx-1+(j-1)*nx]-vym*iCs[nx-1+(j+1)*nx]+vym*iCs[nx-1+j*nx];		//vxm=0

}

////Cs corners
//top left  i=0 j=0
m1 = dtdx * 0.5 /Phi[0]  * (Dm[0]*Phi[0] + Dm[1]*Phi[1]);
m2 = dtdy * 0.5 /Phi[0]  * (Dm[0]*Phi[0] + Dm[nx]*Phi[nx]);
vxm = dt/dx * min(Vwx[0],T(0)) ;
vym = dt/dy * min(Vwy[0],T(0)) ;
Cs[0]= 2*m1*iCs[1]+ 2*m2*iCs[nx]+(1-2*m1-2*m2) * iCs[0]-vxm*iCs[1]+vxm*iCs[0]-vym*iCs[nx]+vym*iCs[0];				//vxp=0 & vyp=0

//bottom left i=0 j=ny-1
m1 = dtdx * 0.5 /Phi[(ny-1)*nx]  * (Dm[(ny-1)*nx]*Phi[(ny-1)*nx] + Dm[1+(ny-1)*nx]*Phi[1+(ny-1)*nx]);
m4 = dtdy * 0.5 /Phi[(ny-1)*nx]  * (Dm[(ny-1)*nx]*Phi[(ny-1)*nx] + Dm[(ny-2)*nx]*Phi[(ny-2)*nx]);
vxm = dt/dx * min(Vwx[(ny-1)*nx],T(0)) ;
vyp = dt/dy * max(Vwy[(ny-1)*nx],T(0)) ;
Cs[(ny-1)*nx]= 2*m1*iCs[1+(ny-1)*nx]+2*m4*iCs[(ny-2)*nx] +(1-2*m1-2*m4) * iCs[(ny-1)*nx]
																			  -vxm*iCs[1+(ny-1)*nx]+vxm*iCs[(ny-1)*nx]-vyp*iCs[(ny-1)*nx]+vyp*iCs[(ny-2)*nx];				//vxp=0, vym =0

//top right i=nx-1 j=0
m2 = dtdy * 0.5 /Phi[nx-1]  * (Dm[nx-1]*Phi[nx-1] + Dm[nx-1+nx]*Phi[nx-1+nx]);
m3 = dtdx * 0.5 /Phi[nx-1]  * (Dm[nx-1]*Phi[nx-1] + Dm[nx-2]*Phi[nx-2]);
vxp = dt/dx * max(Vwx[nx-1],T(0)) ;
vym = dt/dy * min(Vwy[nx-1],T(0)) ;

Cs[nx-1] = 2*m3*iCs[nx-2]+ 2*m2*iCs[nx-1+nx] +(1-2*m2-2*m3) * iCs[nx-1]-vxp*iCs[nx-1]+vxp*iCs[nx-2]-vym*iCs[nx-1+nx]+vym*iCs[nx-1];	//vxm=0, vyp=0

//bottom right i=nx-1 j=ny-1
m3 = dtdx * 0.5 /Phi[nx-1+(ny-1)*nx]  * (Dm[nx-1+(ny-1)*nx]*Phi[nx-1+(ny-1)*nx] + Dm[nx-2+(ny-1)*nx]*Phi[nx-2+(ny-1)*nx]);
m4 = dtdy * 0.5 /Phi[nx-1+(ny-1)*nx]  * (Dm[nx-1+(ny-1)*nx]*Phi[nx-1+(ny-1)*nx] + Dm[nx-1+(ny-2)*nx]*Phi[nx-1+(ny-2)*nx]);
vxp = dt/dx * max(Vwx[nx-1+(ny-1)*nx],T(0)) ;
vyp = dt/dy * max(Vwy[nx-1+(ny-1)*nx],T(0)) ;

Cs[nx-1+(ny-1)*nx] = 2*m3*iCs[nx-2+(ny-1)*nx]+ 2*+m4*iCs[nx-1+(ny-2)*nx] +(1-2*m3-2*m4) * iCs[nx-1+(ny-1)*nx]
																							  -vxp*iCs[nx-1+(ny-1)*nx]+vxp*iCs[nx-2+(ny-1)*nx]-vyp*iCs[nx-1+(ny-1)*nx]+vyp*iCs[nx-1+(ny-2)*nx]; //vxm=0, vym=0

//bore hole
for(int j=int(h/dy); j<=int((e+h)/dy) ; j++)
{
	Cs[j*nx]=0;
}
cout<<" Cs end     ";
}


template <typename T>
T weno5FluxM(T *iCs,T *V,int flag, T maxa)
{
	T C00,C01,C02;
	T C10,C11,C12;
	T C20,C21,C22;

	T p0x,p1x,p2x;
	T b0x,b1x,b2x;
	T epsilon=1e-7;
	T alpha0x,alpha1x,alpha2x,alphasumx;

	T w0x,w1x,w2x;

	T Fpx,Fpx1,Fpx2,Fpxm1,Fpxm2;

	switch(flag){
	case 0:
		Fpx2 =0.5*( V[4]*iCs[4] + fabs(V[4])*iCs[4] );
		Fpx1 =0.5*( V[3]*iCs[3] + fabs(V[3])*iCs[3] );
		Fpx  =0.5*( V[2]*iCs[2] + fabs(V[2])*iCs[2] );
		Fpxm1=0.5*( V[1]*iCs[1] + fabs(V[1])*iCs[1] );
		Fpxm2=0.5*( V[0]*iCs[0] + fabs(V[0])*iCs[0] );
		break;
	case 1:
		Fpx2 =0.5*( V[4]*iCs[4] + fabs(maxa)*iCs[4] );
		Fpx1 =0.5*( V[3]*iCs[3] + fabs(maxa)*iCs[3] );
		Fpx  =0.5*( V[2]*iCs[2] + fabs(maxa)*iCs[2] );
		Fpxm1=0.5*( V[1]*iCs[1] + fabs(maxa)*iCs[1] );
		Fpxm2=0.5*( V[0]*iCs[0] + fabs(maxa)*iCs[0] );
		break;
	default:
		Fpx2 =0.5*( V[4]*iCs[4] + fabs(V[4])*iCs[4] );
		Fpx1 =0.5*( V[3]*iCs[3] + fabs(V[3])*iCs[3] );
		Fpx  =0.5*( V[2]*iCs[2] + fabs(V[2])*iCs[2] );
		Fpxm1=0.5*( V[1]*iCs[1] + fabs(V[1])*iCs[1] );
		Fpxm2=0.5*( V[0]*iCs[0] + fabs(V[0])*iCs[0] );
		break;
	}


	/// right flux
	C00 =  1.0f/3.0f; C01 = -7.0f/6.0f; C02 = 11.0f/6.0f;
	C10 = -1.0f/6.0f; C11 =  5.0f/6.0f; C12 =  1.0f/3.0f;
	C20 =  1.0f/3.0f; C21 =  5.0f/6.0f; C22 = -1.0f/6.0f;

	p0x = C00*Fpxm2 + C01*Fpxm1 + C02*Fpx;
	p1x = C10*Fpxm1 + C11*Fpx + C12*Fpx1;
	p2x = C20*Fpx + C21*Fpx1 + C22*Fpx2;


	//Smooth Indicators, Beta factors

	b0x = (13.0f/12.0f)*(Fpxm2-2.0f*Fpxm1+Fpx)*(Fpxm2-2.0f*Fpxm1+Fpx) + 0.25*(Fpxm2-4.0f*Fpxm1+3.0f*Fpx)*(Fpxm2-4.0f*Fpxm1+3.0f*Fpx);
	b1x = (13.0f/12.0f)*(Fpxm1-2.0f*Fpx+Fpx1)*(Fpxm1-2.0f*Fpx+Fpx1) + 0.25*(Fpxm1-Fpx1)*(Fpxm1-Fpx1);
	b2x = (13.0f/12.0f)*(Fpx-2.0f*Fpx1+Fpx2)*(Fpx-2.0f*Fpx1+Fpx2) + 0.25*(3.0f*Fpx-4.0f*Fpx1+Fpx2)*(3.0f*Fpx-4*Fpx1+Fpx2);

	//
	alpha0x = 0.1f/((epsilon + b0x)*(epsilon + b0x));
	alpha1x = 0.6f/((epsilon + b1x)*(epsilon + b1x));
	alpha2x = 0.3f/((epsilon + b2x)*(epsilon + b2x));
	alphasumx = alpha0x + alpha1x + alpha2x;


	// ENO stencils
	w0x= alpha0x/alphasumx;
	w1x= alpha1x/alphasumx;
	w2x= alpha2x/alphasumx;

	return w0x*p0x + w1x*p1x + w2x*p2x;

}

template <typename T>
T weno5FluxP(T *iCs,T *V,int flag, T maxa)
{
	T C00,C01,C02;
	T C10,C11,C12;
	T C20,C21,C22;

	T p0x,p1x,p2x;
	T b0x,b1x,b2x;
	T epsilon=1e-7;
	T alpha0x,alpha1x,alpha2x,alphasumx;

	T w0x,w1x,w2x;


	T Fmx,Fmx1,Fmx2,Fmxp1,Fmxp2;

	// left fluxes
	C00 = -1.0f/6.0f; C01 =  5.0f/6.0f; C02 =  1.0f/3.0f;
	C10 =  1.0f/3.0f; C11 =  5.0f/6.0f; C12 = -1.0f/6.0f;
	C20 = 11.0f/6.0f; C21 = -7.0f/6.0f; C22 =  1.0f/3.0f;
	switch(flag){
	case 0:
		Fmx2 =0.5*( V[0]*iCs[0] - fabs(V[0])*iCs[0] );
		Fmx1 =0.5*( V[1]*iCs[1] - fabs(V[1])*iCs[1] );
		Fmx	 =0.5*( V[2]*iCs[2] - fabs(V[2])*iCs[2] );
		Fmxp1=0.5*( V[3]*iCs[3] - fabs(V[3])*iCs[3] );
		Fmxp2=0.5*( V[4]*iCs[4] - fabs(V[4])*iCs[4] );
		break;
	case 1:
		Fmx2 =0.5*( V[0]*iCs[0] - fabs(maxa)*iCs[0] );
		Fmx1 =0.5*( V[1]*iCs[1] - fabs(maxa)*iCs[1] );
		Fmx	 =0.5*( V[2]*iCs[2] - fabs(maxa)*iCs[2] );
		Fmxp1=0.5*( V[3]*iCs[3] - fabs(maxa)*iCs[3] );
		Fmxp2=0.5*( V[4]*iCs[4] - fabs(maxa)*iCs[4] );
		break;
	default:
		Fmx2 =0.5*( V[0]*iCs[0] - fabs(V[0])*iCs[0] );
		Fmx1 =0.5*( V[1]*iCs[1] - fabs(V[1])*iCs[1] );
		Fmx	 =0.5*( V[2]*iCs[2] - fabs(V[2])*iCs[2] );
		Fmxp1=0.5*( V[3]*iCs[3] - fabs(V[3])*iCs[3] );
		Fmxp2=0.5*( V[4]*iCs[4] - fabs(V[4])*iCs[4] );
		break;
	}

	p0x = C00*Fmx2 + C01*Fmx1 + C02*Fmx;
	p1x = C10*Fmx1 + C11*Fmx + C12*Fmxp1;
	p2x = C20*Fmx + C21*Fmxp1 + C22*Fmxp2;

	//Smooth Indicators, Beta factors

	b0x = (13.0f/12.0f)*(Fmx2-2.0f*Fmx1+Fmx)*(Fmx2-2.0f*Fmx1+Fmx) + 0.25*(Fmx2-4.0f*Fmx1+3.0f*Fmx)*(Fmx2-4.0f*Fmx1+3.0f*Fmx);
	b1x = (13.0f/12.0f)*(Fmx1-2.0f*Fmx+Fmxp1)*(Fmx1-2.0f*Fmx+Fmxp1) + 0.25*(Fmx1-Fmxp1)*(Fmx1-Fmxp1);
	b2x = (13.0f/12.0f)*(Fmx-2.0f*Fmxp1+Fmxp2)*(Fmx-2.0f*Fmxp1+Fmxp2) + 0.25*(3.0f*Fmx-4.0f*Fmxp1+Fmxp2)*(3*Fmx-4.0f*Fmxp1+Fmxp2);

	//
	alpha0x = 0.3f/((epsilon + b0x)*(epsilon + b0x));
	alpha1x = 0.6f/((epsilon + b1x)*(epsilon + b1x));
	alpha2x = 0.1f/((epsilon + b2x)*(epsilon + b2x));
	alphasumx = alpha0x + alpha1x + alpha2x;

	// ENO stencils
	w0x= alpha0x/alphasumx;
	w1x= alpha1x/alphasumx;
	w2x= alpha2x/alphasumx;

	//total fluxes
	return w0x*p0x + w1x*p1x + w2x*p2x;

}

template <typename T>
T set2limits( T value, T min, T max )
{
	if( value < min ) return min;
	if( value > max ) return max;
	return value;
}


template <typename T>
void BottomAverageFunction(T &ave, T *M,int nx, int ny)
{
	ave = 0;
	for(int i=0; i<nx; i++)
	{
		ave += M[i+(ny-2)*nx];
	}

	ave = ave/nx;

}

template <typename T>
void TopAverageFunction(T &ave, T *M, int nx, int ny)
{
	ave = 0;
	for(int i=0; i<nx; i++)
	{
		ave += M[i];
	}
	ave = ave/nx;
}

template <typename T>
void LeftAverageFunction(T &ave, T *M,int nx, int ny)
{
	ave = 0;
	for(int j=0; j<ny; j++)
	{
		ave += M[1+j*nx];
	}
	ave = ave/ny;
}

template <typename T>
void RightAverageFunction(T &ave, T *M, int nx, int ny)
{
	ave = 0;
	for(int j=0; j<ny; j++)
	{
		ave += M[nx-2+j*nx];
	}
	ave = ave/ny;
}

template <typename T>
void MidAverageFunction(T &ave, T *M, int nx, int ny)
{
	ave = 0;
	for(int j=0; j<ny; j++)
	{
		ave += M[int(nx/2)-1+j*nx];
	}
	ave = ave/ny;
}

template <typename T>
void AverageFunction(T &ave, T *M, int nx, int ny)
{
	ave = 0;
	for(int i=0; i<nx; i++)
	{
		for(int j=0; j<ny; j++)
		{
			ave += M[i+j*nx];
		}
	}
	ave = ave/(ny*nx);
}

template <typename T>
void InjectionAverageFunction(T &ave, T *M, T h, T e, T dy, int nx, int ny)
{
	ave = 0;
	int n=0;
	for(int j=int(h/dy); j<=int((e+h)/dy) ; j++)
	{			//bore hole Dirichlet
		ave += M[2+j*nx];
	ave += M[1+j*nx];
	n+=2;
	}
	ave = ave/n;
}


template <typename T>
void readdata(T *odata,int nx,int ny,const char* s)
{

	ostringstream osfilename;
	osfilename <<  "./fracture_thomas/" << s <<".txt";
	string filename = osfilename.str();
	ifstream myfile( filename.c_str() );
	for( int j = 0; j <ny; j++ )
	{
		for( int i = 0; i < nx; i++ )
		{
			myfile >> odata[i+j*nx];
		}
	}
}



template <typename T>
void writedata(T *idata,int nx,int ny, int tim, string name)
{
	char t[15];
	gcvt(tim,-1,t);
	string  txt;
	txt = ".txt";

	name.append(t);
	name.append(txt);

	const char *namenum;
	namenum=name.c_str();
	ofstream a_file (namenum);

	for (int  j = 0; j <ny; j+=1)
	{
		for (int i = 0; i <nx ; i+=1)
		{
			a_file<< idata[i+j*nx]<<"\t";
		}
		a_file<<endl;
	}
	a_file.close();
}

template <typename T>
void PdWrite(T *Pd, T *Pw, T *Pw0,T d ,int nx,int ny, T dy, int tim, string name)
{
	char t[15];
	gcvt(tim,-1,t);
	string  txt;
	txt = ".txt";
	name.append(t);
	name.append(txt);
	const char *namenum;
	namenum=name.c_str();
	ofstream a_file (namenum);
	int di = int(d);
	int xp = int (nx/2/di);
	a_file<<"#Pressure in depth"<<endl;
	a_file<<"'# depth (m) \t pressure at d=0m (Pa) \t pressure at d=200m (Pa) \t pressure at d=400m (Pa) \t pressure at d=600m (Pa)"<<endl;
	for (int  j = 0; j <ny; j++)
	{
		a_file<<endl<< j*dy;
		for(int i=0; i<=xp ; i++)
		{
			Pd[j] = 0.5 * (Pw[int (nx/2-i*di)+j*nx]+Pw[int (nx/2-i*di)+j*nx]);
			a_file<<"\t\t"<<Pd[j] ;
		}
	}

	a_file.close();
}

template <typename T>
void PressureRelaxFunction(T *Pw,T *Vwx, T *Vwy,T *iPw, T *Rhof ,T *Phi, T *Kappa, T *Viscw ,T *Betaf, T *Betaphi, T *S,T vwx0 ,T dx, T dy, T t, int g_flag, T dt,int nx, int ny)
{
	char pause;
	T m1, m2, m3, m4, m5, m6, m7, m8;
	T dtdy = dt/(dy*dy);
	T dtdx = dt/(dx*dx);
	T g = 9.806;
	T isomass=0;
	T isomass2=0;
#pragma omp parallel for
for(int j=1; j< ny-1; j++)
{
	for(int i=1; i< nx-1; i++)
	{
		Pw[i+j*nx]= dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx] + Kappa[(i+1)+j*nx]/ Viscw[(i+1)+j*nx])*iPw[(i+1)+j*nx]+dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[(i-1)+j*nx]/ Viscw[(i-1)+j*nx])*iPw[(i-1)+j*nx]+ dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j+1)*nx]/ Viscw[i+(j+1)*nx])*iPw[i+(j+1)*nx]+dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j-1)*nx]/ Viscw[i+(j-1)*nx])*iPw[i+(j-1)*nx] +(1-dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx] + Kappa[(i+1)+j*nx]/ Viscw[(i+1)+j*nx])- dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j+1)*nx]/ Viscw[i+(j+1)*nx])-dtdx /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[(i-1)+j*nx]/ Viscw[(i-1)+j*nx])-dtdy /S[i+j*nx]  *0.5*(Kappa[i+j*nx]/ Viscw[i+j*nx]  + Kappa[i+(j-1)*nx]/ Viscw[i+(j-1)*nx])) * iPw[i+j*nx] - dt/dy * 0.5 /(Viscw[i+j*nx] *S[i+j*nx])* g*g_flag * Kappa[i+j*nx] *(Rhof[i+(j+1)*nx] - Rhof[i+(j-1)*nx])  -dt/dy * 0.5 /(Viscw[i+j*nx] *S[i+j*nx])* g*g_flag * Rhof[i+j*nx]
																																																																																																																																																																																																																																																													 *( Kappa[i+(j+1)*nx] -  Kappa[i+(j-1)*nx]) -dt/dy * 0.5 /(S[i+j*nx])* g*g_flag * Rhof[i+j*nx]* Kappa[i+j*nx]  *(1/Viscw[i+(j+1)*nx]-1/Viscw[i+(j-1)*nx]);

		Vwx[i+j*nx] =-0.5 *Kappa[i+j*nx]/Viscw[i+j*nx]*(iPw[(i+1)+j*nx] - iPw[(i-1)+j*nx])/dx;
		Vwy[i+j*nx] =-0.5 *(Kappa[i+j*nx]/Viscw[i+j*nx])*(iPw[i+(j+1)*nx]-iPw[i+(j-1)*nx])/dy + g_flag *(Kappa[i+j*nx]/Viscw[i+j*nx])*Rhof[i+j*nx]*g ;
	}
}

for(int j=1; j< ny-1; j++)
{
	for(int i=1; i< nx-1; i++)
	{
		if(Pw[i+j*nx]<0)
		{
			cout<<"P negative!   "<<i<<"  "<<j<<"  "<<Pw[i+j*nx]<<"  "<<t<<endl;
			pause=getchar();
			Pw[i+j*nx]=iPw[i+j*nx];
		}
	}
}

/// Pressure boundary conditions
#pragma omp parallel for
for(int i=0; i< nx; i++)
{
	//top
	Pw[i] = Pw[i+nx]- g_flag *Phi[i]*Rhof[i]*g*dy;// - g_flag *(1-Phi[i])*2820.0*g*dy;
	//bottom
	Pw[i+(ny-1)*nx] = Pw[i+(ny-2)*nx] + g_flag *Rhof[i+(ny-1)*nx]*Phi[i+(ny-1)*nx]*g*dy;// + g_flag *2820.0*(1-Phi[i+(ny-1)*nx])*g*dy;
}
// 	//bottom
// 		for(int i=0; i< nx; i++)
// 		{
// 			for(int j=1; j< ny-1 ; j++)
// 			{
// 			  isomass += g_flag *Rhof[i+j*nx]*g*dy;
// 			}
// 			Pw[i+(ny-1)*nx] = isomass;
// 			 isomass = 0;
// 		}


for(int j=1; j< ny-1 ; j++)
{

	//left
	Pw[j*nx] = Pw[1+j*nx];
	//right
	Pw[nx-1+j*nx] = Pw[nx-2+j*nx];

}

///Velocity boundary conditions
#pragma omp parallel for
for(int i=1; i< nx-1; i++)
{
	//bottom
	Vwx[i+(ny-1)*nx] =  0.5 *  Kappa[i+(ny-1)*nx] /Viscw[i+(ny-1)*nx]* (iPw[(i+1)+(ny-1)*nx] - iPw[(i-1)+(ny-1)*nx])/dx;
	Vwy[i+(ny-1)*nx] =  -Kappa[i+(ny-1)*nx] /Viscw[i+(ny-1)*nx]* ((iPw[i+(ny-1)*nx] - iPw[i+(ny-2)*nx])/dy -Rhof[i+(ny-1)*nx] *g_flag * g);	//first order
	//top
	Vwx[i] = -0.5 *  Kappa[i] /Viscw[i] * (iPw[i+1] - iPw[i-1])/dx;
	Vwy[i] = -Kappa[i] /Viscw[i]* ((iPw[i+nx] - iPw[i])/dy - Rhof[i+nx]*g_flag* g); 	//first order
}
#pragma omp parallel for
for(int j=1; j< ny-1 ; j++)
{
	//left
	Vwx[j*nx] = -Kappa[j*nx] /Viscw[j*nx] * (iPw[1+j*nx] - iPw[j*nx])/dx;		//first order
	Vwy[j*nx] = - Kappa[j*nx] /Viscw[j*nx] * (0.5 * (iPw[(j+1)*nx]- iPw[(j-1)*nx])/dy - Rhof[j*nx]*g_flag* g);
	//right
	Vwx[nx-1+j*nx] =  - Kappa[nx-1+j*nx] /Viscw[nx-1+j*nx] * (iPw[nx-1+j*nx] - iPw[(nx-2)+j*nx] )/dx;
	Vwy[nx-1+j*nx] =  - Kappa[nx-1+j*nx] /Viscw[nx-1+j*nx]* (0.5 *(iPw[nx-1+(j+1)*nx]- iPw[nx-1+(j-1)*nx])/dy -Rhof[nx-1+j*nx]*g_flag * g);
}


///corners
///top left
Vwx[0] =  -Kappa[0] /Viscw[0]* (iPw[1] - iPw[0])/dx;		//first order
Vwy[0] =  -Kappa[0] /Viscw[0] * ((iPw[nx] - iPw[0])/dy - Rhof[0]*g_flag* g); 	//first order
///bottom left
Vwx[(ny-1)*nx] = -Kappa[(ny-1)*nx] /Viscw[(ny-1)*nx] * (iPw[1+(ny-1)*nx] - iPw[(ny-1)*nx])/dx;
Vwy[(ny-1)*nx] = -Kappa[(ny-1)*nx] /Viscw[(ny-1)*nx] * ((iPw[(ny-1)*nx] - iPw[(ny-2)*nx])/dy - Rhof[(ny-1)*nx]*g_flag* g);
///top right
Vwx[nx-1] =  -Kappa[nx-1] /Viscw[nx-1]* (iPw[nx-1] - iPw[(nx-2)])/dx;
Vwy[nx-1] =  -Kappa[nx-1] /Viscw[nx-1] * ((iPw[nx-1+nx] - iPw[nx-1])/dy - Rhof[nx-1]*g_flag* g);
///bottom right
Vwx[nx-1+(ny-1)*nx] =  -Kappa[nx-1+(ny-1)*nx] /Viscw[nx-1+(ny-1)*nx] * (iPw[nx-1+(ny-1)*nx] - iPw[(nx-2)+(ny-1)*nx] )/dx;
Vwy[nx-1+(ny-1)*nx] =  -Kappa[nx-1+(ny-1)*nx] /Viscw[nx-1+(ny-1)*nx] * ((iPw[nx-1+(ny-1)*nx] - iPw[nx-1+(ny-2)*nx])/dy - Rhof[nx-1+(ny-1)*nx]*g_flag* g);


#pragma omp parallel for
for(int i=0; i< nx; i++)
{
	//bottom
	Vwx[i+(ny-1)*nx] =  Vwx[i+(ny-2)*nx];
	Vwy[i+(ny-1)*nx] =  Vwy[i+(ny-2)*nx] ;
	//top
	Vwx[i] = Vwx[i+nx] ;
	Vwy[i] = Vwy[i+nx] ;
}
#pragma omp parallel for
for(int j=1; j< ny-1 ; j++)
{
	//left
	Vwx[j*nx] =Vwx[1+j*nx] ;
	Vwy[j*nx] = Vwy[1+j*nx];
	//right
	Vwx[nx-1+j*nx] =  Vwx[nx-2+j*nx];
	Vwy[nx-1+j*nx] =  Vwy[nx-2+j*nx];
}

}



template void readdata(double *odata, int nx, int ny, const char* s);

template void writedata(double *idata,int nx,int ny, int tim, string name);

template void PdWrite(double *Pd, double *Pw, double *Pw0,double d ,int nx,int ny,double dy, int tim, string name);

template void MemoryCopy(double *in, double *out, int nx, int ny);

template void PorosityFunction( double  *Phi,double  *iPhi, double  *Pw, double  *Sxx, double  *Szz, double  phiini, double  phir,  double  phia, int nx, int ny);

template void PermeabilitySinglePhase(double *Kappa, double *Phi, double *Kappa0, double *Pw,double *Sxx, double *Szz, double phiini, double cap_depth, double k_cap,double c,double h, double e, double kappa_hole,double dx, double dy,  int nx, int ny);

template void ThermalConductivityFunction(double *TCm, double *Frac,double *Rcp,double *Phi, double  *Rhow, double *Cpw, double tcr, double  rhor0, double cpr, double  tcw, int nx, int ny);

template void StateEqTemperetureDensityTransport(double *Viscw , double *Rhof, double *Tw, double *Pw, double *Cs, double *Betaw,double *Cpw, double rhof0,  int nx, int ny);

template void LTEDensityTransportSML(double  *Tr, double  *iTr,double  *Tw, double  *Vwx,  double  *Vwy,  double  *Phi,  double  *Rhow, double *Rcp, double  *Cpw, double  *TCr, double *x,double *y,  double  h, double e, double  tw_fresh, double  dx, double  dy, double t,double  dt,int nx, int ny);

template void LTNENonIsoThermal(double *Tw, double *Tr, double *iTw, double *iTr, double *Vwx, double *Vwy, double *Phi, double *Rhow, double *Cpw, double tcw, double rhor0, double cpr, double tcr, double rockr ,double *x, double *y,double tw_fresh, double dx, double dy, double dt, int nx, int ny);

template void LTNENonIsoThermalTR(double *iTw, double *Tr, double *iTr,double *Phi,double *Hta,double t,double rhor0,double cpr,double tcr,double h,double e,double tw_fresh, double dx,double dy,double dt, int nx, int ny);

template void LTNENonIsoThermalTF(double *Tw,double *iTw,double *iTr,double *Vwx,double *Vwy,double *Phi,double *Rhow,double *Hta,double *Cpw,double t,double tcw,double *x,double *y,double h,double e,double tw_fresh, double dx,double dy,double dt,int nx,int ny);

template void Hfs(double *Hta ,double *TCr,double *Phi, double *Re, double tcw,double rockr ,int nx, int ny);

template void EGSPressureFunction(double *Pw,double *Vwx, double *Vwy,double *iPw, double *Rhof, double *Phi, double *Kappa, double *Viscw ,double *Betaf, double *Betaphi,double *S,double *Re,double *Pw0,double *Y, double rockr ,  double h,  double e, double pmax, double pmin, double vwx0 ,double dx, double dy, double t ,int g_flag, double dt,int nx, int ny);

template void EGSMassFraction(double *Cs, double *iCs, double *Dm, double *Vwx, double *Vwy, double *Phi, double rhow0, double rho_max,double h,double e, double dx, double dy, double dt, int nx, int ny);

template void TimeStepCalculation(double &dt,double *Phi, double *Kappa, double *Viscw, double *Betaphi, double *Betaw, double *Rhow,double *TCm,double *Dm,double *Cpw,double tcr,double tcw,double cpr,double rhor0,double dx,double dy,int nx, int ny);

template  void AdvectionGamma( double *Gamma, double *Alpha, double *Rhof ,double *Cpw, double *Phi,double rhor0, double cpr,double dx, double dy, int nx, int ny);

template void AdvectionCoefficient(double *Adv, double *Rhof ,double *Cpw, double *Phi,double rhor0, double cpr, int lte_flag,int nx, int ny);

template void Fracture(double *K, double *Phi,double *TCr, double *Pw, double *Szz, double *Sxx, double k_fracture, double perm_max,double phi_fracture,double depth_fracture, double width_fracture,double theta, double fracture_aperture, double l_fracture, double tcr_fracture,double phiini, double phir, double phia, double c, double PI, double dx, double dy, int nx, int ny);

template  void SaveDistAve(double *M0,double *M1,double *M2,double *M3,double *M4,double *M5,double *M6,double *M7,double *M8, string name, double dx,double dy,  int tim, int nx,int ny);

template void SaveVector(double *Mx, double *My, string name, double dx, double dy, int tim, int nx, int ny);

template void SaveVtk(double *Pw,double *Cs,double *Tf,double *Tr,double *Kappa,double *Phi,double *Vx,double *Vy, string name,double dx, double dy, int tim, int nx, int ny);

template void DensityTransportSML(double *Cs,double *iCs,double *Dm,double *Phi,double *Vwx,double *Vwy,double *x,double *y,double *Rhof, double h,double e,double c_fresh, double dx,double dy,double t, double teq,double dt,int nx,int ny);

template double spline_cubic(const double a[4],double x);

template void InterpolationPoints(int *xind,int *yind,double &npox,double &npoy,double *x,double *y,int nx,int ny);

template double weno5FluxP(double *iC,double *V, int flag, double maxa);

template double weno5FluxM(double *iC,double *V, int flag, double maxa);

template double set2limits( double value, double min, double max );

template void BottomAverageFunction(double  &ave, double  *M,int nx, int ny);

template void TopAverageFunction(double  &ave, double  *M, int nx, int ny);

template void LeftAverageFunction(double  &ave, double  *M,int nx, int ny);

template void RightAverageFunction(double  &ave, double  *M, int nx, int ny);

template void MidAverageFunction(double  &ave, double  *M, int nx, int ny);

template void AverageFunction(double &ave,double *M, int nx, int ny);

template void InjectionAverageFunction(double &ave,double *M, double h, double e, double dy, int nx, int ny);

template void PressureRelaxFunction(double *Pw,double *Vwx,double *Vwy,double *iPw, double *Rhof ,double *Phi, double *Kappa,double *Viscw ,double *Betaf, double *Betaphi, double *S, double vwx0 ,double dx, double dy, double t, int g_flag, double dt,int nx, int ny);

