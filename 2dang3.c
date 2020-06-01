/* New version of 2dang.c to include cell proliferation and death and diffusion of fibronectin (i.e. Crank-Nicholson, Gauss-Siedel iteration) */
/* NOTE: cell meshpoints m=0 represent the capillary, m>0 are the ECM */
/* NOTE: substrate meshpoints m=0 are in the capillary, m=1 are the capillary wall and m>1 are in the ECM */

/* NB Increments of m and n in subecm() are increased to reduce the amount of data sent to file */

/* NOTE: because x-range is [0,1] and y-range is [0,l/L], we only use the first (l/L)*nn or (l/L)*vv meshpoints in the y-direction */
/* y-range in procedures subecm and trailsend is [0,l/L] - alter manually */

/* Transition probabilities can depend on VEGF if required (gamma5, gamma6 != 0) */

/* EC diff coeff x100  (changes many values) */
/* Transition probability exponents are x25 */
/* Transition probabilities are in terms of protease and fibronectin, but protease decay is also included */
/* Protease production multiplication factor is 100 (was 10 but had to increase because of decay term) */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#define pi 3.141592654		/* Constant Pi */
#define TRUE 1
#define FALSE 0


/* Set parameters */

#define uu 100		/* Upper bound for number of cells */
#define jj 21600		/* Number of time steps (standard is 900 for 12 hours when nn=51, 14400 for 12 hours when nn=201) */
#define nn 201		/* Number of space steps in each direction (must be odd, h=2/(nn-1) )     (std. 51) */
#define vv 401		/* Number of horizontal (y=const) VEGF meshlines  =2nn-1 */

#define lnn 101		/* Lowered nn for use in y-direction =(l/L)*nn      (std. 26) */
#define lvv 201		/* Lowered vv for use in y-direction =2*lnn-1 (must be odd) */

#define subinc 4

/* Dimensionless parameters */
#define hours 694.4444444	/* Time scale in hours, i.e. this is the number of hours equivalent to a dimensionless time of 1.0 */
#define length 0.06912		/* Time length of simulation: time in hours is length*hours = length*L^2/D_p */
#define fi 1.0			/* fi/f0 (Scaled) Initial fibronectin */
#define f1 0.6			/* f1/f0 (Scaled) Threshold fibronectin for penetration into ECM */
#define gamma1 100.0		/* Exponent for protease in capillary (4.0) */
#define gamma2 100.0		/* Exponent for fibronectin in capillary (4.0) */
#define gamma3 50.0		/* Exponent for protease in ECM (2.0) */
#define gamma4 37.5		/* Exponent for fibronectin in ECM (1.5) */
#define gamma5 0.0		/* Expnent for VEGF in capillary (40.0 in runs) */
#define gamma6 0.0		/* Exponent for VEGF in ECM (20.0 in runs) */
#define v1 0.1			/* 0.1 in runs */
#define v2 10.0			/* 10.0 in runs */
#define v3 0.1			/* 0.1 in runs */
#define v4 10.0			/* 10.0 in runs */
#define m0 12.0			/* Exponent in VEGF source term */
#define m1 2.0
#define a1 0.0			/* A1 coefficient of rate of increase of VEGF at capillary in capillary source term */
#define k1 51923.077
#define k2 694.4444444		/* 694.44444444444 */
#define k3 3166.6666667		/* 3166.66666667 */
#define k4 154.32099999
#define k5 1740294.5		/* 1740294.5 */
#define k6 0.0120481
#define k7 1.0			/* NB No account is currently taken of this value in the program: would have to modify lambda in the ECM if k7 is not equal to 1 */
#define k8 10.000
#define k9 0.0001
#define k10 154.32099
#define k11 0.00076923077
#define k12 0.0076923077		/* Note k12>k11 so EC are attracted to areas of HIGH protease */
#define k13 100.0
#define k14 10.0		/* Note k14<k13 so EC are attracted to areas of LOW fibronectin */
#define k15 0.00076923077
#define k16 0.0076923077
#define k17 100.0
#define k18 50.0
#define k19 0.00036
#define k20 1.0
#define k21 1.294620097
#define k22 5736.9
#define k23 0.00000076923077
#define k24 0.00001859
#define k25 3.4722222
#define k26 38.88888888
#define pmax 6.0		/* Maximum diemsnionless cell density before logistic proliferation stops (NB initial density = 1) */
#define q 0.71428571		/* 0.71428571 */
#define profact 1.0		/* Factor of increased protease production */
#define child 0.25		/* (12 hours) Length of time (as a proportion of the simulation length) for which a cell is prevented from dividing after its last division */

#define peaks 1.0		/* Number of peaks (integer) in the VEGF source function */
#define relax 1.3		/* Relaxation parameter for Gauss-Siedel SOR for VEGF: 1.3 for nn=51, 1.45 for nn=201 */
#define relax2 1.0		/* Relaxation parameter for Gauss-Siedel SOR for fibronectin */
#define tolerance 0.000001	/* Tolerance for Gauss-Siedel iteration */


/* Declare global variables */

/* Total number of cells */
int ii;

/* Horizontal and vertical position (cell meshpoint number) of cell i at timestep j */
int hpos[uu][jj],vpos[uu][jj];

/* Fibronectin concentration at VEGF meshpoint (n,m)
Here y=m/(nn-1)-1 and when m is even, x=(2n+1)/(nn-1)-1 (nn-1 meshpoints); when m is odd, x=2n/(nn-1)-1 (nn meshpoints). */
double vegf[lvv][nn];
double pro[lvv][nn];
double fib[lvv][nn];		/* Note arguments are in reverse order to increase running speed */
double vegfold[lvv][nn];
double proold[lvv][nn];
double fibold[lvv][nn];
double maxsub[6];	/* max[0]=vegfmax,  max[1]=fibmax,  max[2]=promax,  max[3]=vegfcapmax,  max[4]=fibcapmax,  max[5]=procapmax  */

double denscale;

/* Status of (no of cells occupying) cell meshpoint (n,m) [ x=2n/(nn-1) - 1 , y=2m/(nn-1) - 1 ] at time j */
int occ[nn][nn];
int occold[nn][nn];

/* Time step at which cell i leaves the simulation (by death or by reaching the tumour) (=jj-1 by default), time step at which cell i last divided and time step at which cell i was born */
int deathtime[uu],divtime[uu],birthtime[uu];

int capin,capout,up,down,in,out,births,deaths,breach,maxcell;
int vegfwarn,fibwarn,prowarn;




/* Declare functions */
void writeparams(int,int,int),writecellscap(int),writetrails(),writesub(int),writebatch();
void movecell(int,int,double,double,double,double,double),updatevegf(),updatefib(),updatepro(int);
void trailsintro(FILE *,int,char *),trailsend(FILE *,float);
void cellscap(int *,FILE *,int,char *),subecm(int,double *,FILE *,int,char *),subcap(int,double *,FILE *,int,char *);
double pstill(int),pmove(int,int,int),lambda(),tau(double,double,double,int),h(),k(),vegf0(double,double),cellxy(int),vy(int),vx(int,int),fib0(double,double),pro0(double,double);
int go(int,int,int *),heaviside(double),max(int,int),min(int,int),optatoi(char *),init();
float ran0(int *);


int main(int argc,char *argv[])
{
	int starttime=(int)(unsigned)time(NULL),idum=-starttime;
	int option;
	int run,runs;
	FILE *outp;
	
	ran0(&idum);		/* Seed random number generator */
	
	if(argc==3)
	{
		runs=atoi(argv[1]);
		option=optatoi(argv[2]);
		
		if(option!=2)
		{
			/* Initiate params.dat file */
			if((outp=fopen("params.dat","w"))==NULL)
				printf("Cannot open output file params.dat\n");		
			else
			{	
				fprintf(outp,"params.dat\n\n");
				fclose(outp);
			}
			
			for(run=0;run<runs;run++)
			{
				writeparams(go(init(),option,&idum),starttime,run);		/* This executes init and then go before write and passes the result of go as the first argument to write */
				printf("\n");
			}
			if(option==1)			/* Skip if no files option was specified */
			{
				writetrails();
				writebatch();
			}		
			printf("\n%s executed %i runs successfully (output option %s)\n\n",argv[0],runs,argv[2]);
		}
		else
		{
			printf("ERROR! Program requires two arguments:\n");
			printf("1. Number of runs to be performed\n");
			printf("2. Trails or no files (t/n)\n");
		}
	}
	else
	{
		printf("ERROR! Program requires two arguments:\n");
		printf("1. Number of runs to be performed\n");
		printf("2. Trails or no files (t/n)\n");
	}
	

	
	return(0);
}




/* INITIAL DATA */
int init()
{	
	int n,m;
	double x,y;
	FILE *outp;
		
	ii=0;		/* Start with no cells */
	
	/* Reset variables */
	capin=0;
	capout=0;
	up=0;
	down=0;
	in=0;
	out=0;
	births=0;
	deaths=0;
	breach=FALSE;
	maxcell=0;
	vegfwarn=0;
	fibwarn=0;
	prowarn=0;
	for(n=0;n<6;n++)
		maxsub[n]=0.0;

	for(n=0;n<nn;n++)
	{
		for(m=0;m<nn;m++)
			occ[m][n]=0;
	}
	
	/* Cycle through capillary cell meshpoints */
	for(n=19;n<nn;n=n+40)
	{
		/* Seed cell at meshpoint (n,0) */
		hpos[ii][0]=n;
		vpos[ii][0]=0;
		occ[0][n]++;
		deathtime[ii]=jj-1;	/* jj-1 because on the jjth time step, j=jj-1.   If it was jj it would be out of range! */
		divtime[ii]=-jj;		/* Allows cells to start proliferating straight away (as soon as out of capillary) */
		birthtime[ii]=0;
		ii++;	/* Count total number of cells */
	}
	
	printf("Cells ii = %i \n",ii);

	
	/* Cycle through substrate meshpoints */
	for(m=0;m<lvv;m++)
	{
		y=vy(m);
	
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
			{
				x=vx(n,m);
				vegf[m][n]=vegf0(x,y);
				fib[m][n]=fib0(x,y);
				pro[m][n]=pro0(x,y);
				vegfold[m][n]=vegf0(x,y);
				fibold[m][n]=fib0(x,y);
				proold[m][n]=pro0(x,y);
				
			}
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
			{
				x=vx(n,m);
				vegf[m][n]=vegf0(x,y);
				fib[m][n]=fib0(x,y);
				pro[m][n]=pro0(x,y);
				vegfold[m][n]=vegf0(x,y);
				fibold[m][n]=fib0(x,y);
				proold[m][n]=pro0(x,y);
			}
		}
	}

	/* Initialise totprocap.dat file */
	if((outp=fopen("totprocap.dat","w"))==NULL)
		printf("Cannot open output file totprocap.dat\n");		
	else
	{	
		fprintf(outp,"totprocap.dat\n\n");
		fclose(outp);
	}
	
	
	denscale=(double)nn/ii;
	return(ii);
}
	



/* RUN SIMULATION */
int go(int cellsincap,int option,int *idum)
{
	
	/* Declare variables */
	int i,j,n,m;
	int count=0;
	double a,b,c,d,r,fibcap,p1,p2,cm0,cm1,logistic,g;
	char z[2];

	/* Cycle through time steps */
	for(j=0;j<jj-1;j++)	/* Only go to j<jj-1 as routine writes information hpos[i][j+1] */
	{
		/* printf("Time = %i  %i  %i\n",j,ii,deaths); */
		
		if(j%(int)(jj/10)==0) printf("Time = %i\n",j);
	
		if(strcmp(z,"k")!=0)
			strcpy(z,"");
			
			
		/* Copy occ variable to occold */
		for(n=0;n<nn;n++)
		{
			for(m=0;m<nn;m++)
				occold[m][n]=occ[m][n];
		}

		
		/* Cycle through cells */
		for(i=0;i<ii;i++)
		{

			/* Assume no cell movement for now */
			n=hpos[i][j];
			m=vpos[i][j];
			hpos[i][j+1]=n;
			vpos[i][j+1]=m;
			
			/* This section only executed if cell is still in simulation and has left the capillary */
			if(deathtime[i]==jj-1 && m>0)
			{
			
		
				/* Cell leaves simulation when it reaches maximum value of y (but it stays counted in the occ variable?) */
				if(m==lnn-1)
				{
					deathtime[i]=j;
					occ[m][n]--;
				}
					

				/* Calculate proliferation and death probabilities */
				/* First calculate average protease at the surrounding substrate meshpoints at times j and j-1 */
				/* count is how many surrounding meshpoints there are - usually 4 but less for a boundary point */
				cm0=0;
				cm1=0;
				count=0;

				if(n>0)
				{
					cm0=cm0+pro[2*m][n-1];	/* pro to the left */
					cm1=cm1+proold[2*m][n-1];
					count++;
				}	
				if(n<nn-1)
				{
					cm0=cm0+pro[2*m][n];		/* pro to the right */
					cm1=cm1+proold[2*m][n];
					count++;
				}
				if(m>0)
				{
					cm0=cm0+pro[2*m-1][n];	/* pro below */
					cm1=cm1+proold[2*m-1][n];
					count++;
				}
				if(m<lnn-1)
				{
					cm0=cm0+pro[2*m+1][n];	/* pro above */
					cm1=cm1+proold[2*m+1][n];
					count++;
				}
				
				cm0=cm0/count;
				cm1=cm1/count;

				if(cm1>cm0)
					cm1=cm0;		/* If protease level has dropped, no proliferation contribution */
					
				if(denscale*occ[m][n]>pmax)		/* Logistic term in expression for p: but do not allow to be <0 */
					logistic=0;
				else
					logistic=(double)1-denscale*occ[m][n]/pmax;
					
				g=k22*exp(-k24*pow(cm0,m1))*(1-k24*m1*pow(cm0,m1))/(1+k22*cm0*exp(-k24*pow(cm0,m1)));
				if(g<0) g=0;			
				
				/* Second heaviside function doesn't allow proliferation if last proliferation was less than (child*jj) timesteps ago */
				p1=q * heaviside(cm0-profact*k23) * heaviside((double)j-divtime[i]-child*jj) * (k()*k26*logistic + (cm0-cm1)*g);		/* Probability of dividing */
				p2=k()*q*k25;		/* Probability of dying */
				
				if(p1<(double)0)
				{
					printf("n %i m %i  cm0 %9.6f  cm1 %9.6f   p1 %9.6f\n",n,m,cm0,cm1,p1);
					gets(z);
					p1=(double)0;
				}
				if(p1>(double)1-p2)
				{
					printf("n %i m %i  cm0 %9.6f  cm1 %9.6f   p1 %9.6f\n",n,m,cm0,cm1,p1);
					gets(z);
					p1=(double)1;
				}
				
				/* if(strcmp(z,"s")!=0 && strcmp(z,"k")!=0)		
					printf("Birth %9.6f   Death %9.6f\n",p1,p2);
				*/

				/* Generate random number to simulate proliferation/death */
				r=(double)ran0(idum);

				if(r<p2)		/* Cell dies */
				{
					deathtime[i]=j;
					occ[m][n]--;
					deaths++;
					printf("Death: time %i, cell %i, n=%i, m=%i\n",j,i,hpos[i][j],vpos[i][j]);
				}
				else
				if(r<p1+p2)		/* Cell divides */
				{
					if(ii==uu)
					{
						printf("Max cell number reached!\n");
						break;
					}
					else
					{
					hpos[ii][j]=n;
					vpos[ii][j]=m;		/* NOTE: Time step = j because cell has opportunity to move after this code (but in same time step) and input values are hpos[i][j] etc. */
					occ[m][n]++;		
					deathtime[ii]=jj-1;
					divtime[ii]=j;
					birthtime[ii]=j;
					divtime[i]=j;
					ii++;
					births++;
					printf("Birth: time %i, cell %i, n=%i, m=%i\n",j,i,hpos[i][j],vpos[i][j]);
					}
				}
			}

			if(deathtime[i]==jj-1)		/* Only proceed if cell is still in the simulation */
			{
				
				
	
		
				/* Generate transition probabilities and random number between 0 and 1 */
				a=pstill(m);		/* Probability of staying still */
				b=pmove(n,m,0);	/* Probability of moving left */
				c=pmove(n,m,1);	/* Probability of moving right */
				d=pmove(n,m,2);	/* Probability of moving down */
				r=(double)ran0(idum);
				
				if(a<0)
					printf("Warning: negative probabilities\n");

				/* Print probabilities: s skips printout for current timestep, k skips for all */
				/* if(strcmp(z,"s")!=0 && strcmp(z,"k")!=0)		
				{
					printf("%2d %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",i,a,b,c,d,(double)1-a-b-c-d,r);
					gets(z);
				}
				*/


				if(m==0)	/* Check if cell can be released from the capillary into the ECM */
				{
					if(n==0)
						fibcap=fib[0][0];
					else
					if(n==nn-1)
						fibcap=fib[0][nn-2];
					else
						fibcap=(fib[0][n-1]+fib[0][n])/2;								

					if(fibcap<f1)		/* If average fibronectin value in capillary is below threshold value f1, cell moves up into ECM  */
					{
						r=(double)2;
						breach=TRUE;
						cellsincap--;
						printf("BREACH: j=%i n=%i cellsincap=%i\n",j,n,cellsincap);
					}
				}

				movecell(i,j,a,b,c,d,r);
			
			}
		}
		
		updatevegf();
		updatefib();
		updatepro(j);
		
		if( ((5*j)%jj==0 || j==jj-2) && option==1 )
		{
			writesub(j);
			if(cellsincap>0)
				writecellscap(j);
		}
				

	}

	/* Count and return number of cells left in simulation at end */
	count=0;
	for(i=0;i<ii;i++)
		if(deathtime[i]==jj-1)
			count++;
	return(count);		
}




/* MOVE CELL */
void movecell(int i,int j,double a,double b,double c,double d,double r)
{
	/* Cell moves up */
	if(r>a+b+c+d)
	{
		up++;
		vpos[i][j+1]=vpos[i][j]+1;
	}

	/* Cell moves down */
	else
	if(r>a+b+c)
	{
		down++;
		vpos[i][j+1]=vpos[i][j]-1;
	}
	
	/* Cell moves right */
	else
	if(r>a+b)
	{
		if(vpos[i][j]==0)
		{
			if(hpos[i][j]<(float)nn/2)
				capin++;
			else
				capout++;
		}
		else
		{
			if(hpos[i][j]<(float)nn/2)
				in++;
			else
				out++;
		}
		hpos[i][j+1]=hpos[i][j]+1;
	}
	
	/* Cell moves left */
	else
	if(r>a)
	{
		if(vpos[i][j]==0)
		{
			if(hpos[i][j]<(float)nn/2)
				capout++;
			else
				capin++;
		}
		else
		{
			if(hpos[i][j]<(float)nn/2)
				out++;
			else
				in++;
		}
		hpos[i][j+1]=hpos[i][j]-1;
	}
	
	if(r>a)		/* If cell has moved (in any direction) */
	{
		occ[vpos[i][j]][hpos[i][j]]--;
		occ[vpos[i][j+1]][hpos[i][j+1]]++;
	}

}



/* GAUSS-SIEDEL ITERATION: CALCULATES VEGF VALUES AT TIME STEP j+1 */
void updatevegf()
{
	double v[lvv][nn];
	double vold;
	double p;
	int n,m;
	int intol=0;
	int count=0;
	double cvj,cvjm1,vdiff;
	
	/* Cycle through substrate meshpoints */
	/* First do m=0, capillary. m is even so number of substrate meshpoint in x-direction is nn-1 */
	for(n=0;n<nn-1;n++)
	{
		
		/* Average density at cell meshpoints to the right and left of substrate meshpoint at time j, use occold and not occ because occ has already been updated this timestep */
		/* Scaled by denscale because p should really be cell density, not just number of cells at a meshpoint */
		p=(double)denscale*(occold[0][n]+occold[0][n+1])/2;	
		cvj=(vegf[1][n]+vegf[1][n+1])/2;			/* Average VEGF level on capillary wall (m=1) at time j */
		cvjm1=(vegfold[1][n]+vegfold[1][n+1])/2;		/* Average VEGF level on capillary wall (m=1) at time j-1 */
		
		vdiff=cvj-vegf[0][n];					/* Difference between VEGF level on capillary wall (which is classed as ECM) and in capillary */
		if(vdiff<0)
			vdiff=0;

		vegfold[0][n]=vegf[0][n];			
		vegf[0][n] = vegf[0][n] + k()*(-k1*vegf[0][n]*p/(1+vegf[0][n]) + k2*vdiff) + a1*(cvj-cvjm1);
		if(vegf[0][n]<(double)0)
		{
			vegf[0][n]=(double)0;
			vegfwarn++;
		}
	}
	
	
	/* Cycle through VEGF meshpoints and set initial guess for v to be value for previous time step (j) */
	for(m=1;m<lvv;m++)
	{
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
				v[m][n]=vegf[m][n];
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
				v[m][n]=vegf[m][n];
		}
	}


	while(intol==0)	/* Keep iterating until the value of v at each meshpoint changes by less than a certain tolerance */
	{

		intol=1;
		count++;
		
		/* Cycle through VEGF meshpoints and make improved estimate for v - treating boundary meshpoints separately - and checking tolerance at each update */
		/* NOTE: Rows are treated in reverse order (and comments read in reverse order) */
		
		/* Same routine for m=lvv-1 */
		vold=v[lvv-1][0];
		v[lvv-1][0]=v[lvv-1][1];
		if(v[lvv-1][0]-vold>tolerance || v[lvv-1][0]-vold<-tolerance)
			intol=0;
		for(n=1;n<nn-2;n++)
		{
			vold=v[lvv-1][n];
			v[lvv-1][n]=k21*h()*pow(1-cos(2*peaks*pi*vx(n,lvv-1)),m0) + v[lvv-3][n];
			if(v[lvv-1][n]-vold>tolerance || v[lvv-1][n]-vold<-tolerance)
				intol=0;
		}
		vold=v[lvv-1][nn-2];
		v[lvv-1][nn-2]=v[lvv-1][nn-3];
		if(v[lvv-1][nn-2]-vold>tolerance || v[lvv-1][nn-2]-vold<-tolerance)
			intol=0;

		/* And last of all treat final two rows, m=lvv-2 and m=lvv-1, in similar way to first two (lvv is odd so lvv-2 is odd, lvv-1 is even) */
		vold=v[lvv-2][0];
		v[lvv-2][0]=v[lvv-2][1];
		if(v[lvv-2][0]-vold>tolerance || v[lvv-2][0]-vold<-tolerance)
			intol=0;
		for(n=1;n<nn-1;n++)
		{
			vold=v[lvv-2][n];
			v[lvv-2][n]=k21*h()*pow(1-cos(2*peaks*pi*vx(n,lvv-2)),m0) + v[lvv-4][n];	/* Note source VEGF, function of x and t */
			if(v[lvv-2][n]-vold>tolerance || v[lvv-2][n]-vold<-tolerance)
				intol=0;
		}
		vold=v[lvv-2][nn-1];
		v[lvv-2][nn-1]=v[lvv-2][nn-2];
		if(v[lvv-2][nn-1]-vold>tolerance || v[lvv-2][nn-1]-vold<-tolerance)
			intol=0;

		/* Now cycle through interior rows */
		for(m=lvv-3;m>2;m--)
		{
			/* First meshpoint on row requires BC at x=0 */
			vold=v[m][0];
			v[m][0]=v[m][1];
			if(v[m][0]-vold>tolerance || v[m][0]-vold<-tolerance)
				intol=0;
		
			if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 (but don't include last meshpoint in loop - boundary element) */
			{
				for(n=1;n<nn-2;n++)
				{
					p=(double)pow(denscale,2)*(occold[m/2][n]+occold[m/2][n+1])/2;	/* Average density at cell meshpoints to the right and left of VEGF meshpoint at time j. Here scale by denscale^2 because density in ECM is cells per unit AREA not LENGTH as in capillary */
					vold=v[m][n];
					v[m][n]=relax/(h()*h()+2*k8*k()) * (0.5*k8*k()*(v[m][n+1]+v[m][n-1]+v[m+2][n]+v[m-2][n]+vegf[m][n+1]+vegf[m][n-1]+vegf[m+2][n]+vegf[m-2][n]) + (h()*h()-2*k8*k()-h()*h()*k()*k1*p/(1+vegf[m][n]))*vegf[m][n]) + (1-relax)*v[m][n];
					if(v[m][n]-vold>tolerance || v[m][n]-vold<-tolerance)
						intol=0;
				}
				vold=v[m][nn-2];
				v[m][nn-2]=v[m][nn-3];
				if(v[m][nn-2]-vold>tolerance || v[m][nn-2]-vold<-tolerance)
					intol=0;
			}
			
			else		/* If m is odd, number of subsrate meshpoints in x-direction is nn (but don't include last meshpoint in loop) */
			{
				for(n=1;n<nn-1;n++)
				{
					p=(double)pow(denscale,2)*(occold[(m-1)/2][n]+occold[(m+1)/2][n])/2;	/* Average density at cell meshpoints above and below VEGF meshpoint at time j */
					vold=v[m][n];
					v[m][n]=relax/(h()*h()+2*k8*k()) * (0.5*k8*k()*(v[m][n+1]+v[m][n-1]+v[m+2][n]+v[m-2][n]+vegf[m][n+1]+vegf[m][n-1]+vegf[m+2][n]+vegf[m-2][n]) + (h()*h()-2*k8*k()-h()*h()*k()*k1*p/(1+vegf[m][n]))*vegf[m][n]) + (1-relax)*v[m][n];
					if(v[m][n]-vold>tolerance || v[m][n]-vold<-tolerance)
						intol=0;
				}
				vold=v[m][nn-1];
				v[m][nn-1]=v[m][nn-2];
				if(v[m][nn-1]-vold>tolerance || v[m][nn-1]-vold<-tolerance)
					intol=0;
			}
		}


		/* Same routine for m=2, which is also classed as a boundary row */
		vold=v[2][0];
		v[2][0]=v[2][1];
		if(v[2][0]-vold>tolerance || v[2][0]-vold<-tolerance)
			intol=0;
		for(n=1;n<nn-2;n++)
		{
			vold=v[2][n];
			v[2][n]=1/(k19+h()*k20) * (k19*v[4][n]+h()*vegfold[0][n]);	/* Use vegfold because capillary values have already been updated */
			if(v[2][n]-vold>tolerance || v[2][n]-vold<-tolerance)
				intol=0;
		}
		vold=v[2][nn-2];
		v[2][nn-2]=v[2][nn-3];
		if(v[2][nn-2]-vold>tolerance || v[2][nn-2]-vold<-tolerance)
			intol=0;
		
		
		/* First treat m=1 (capillary wall), n=0 using BC at x=0 (simpler condition than capillary BC) */
		vold=v[1][0];
		v[1][0]=v[1][1];
		if(v[1][0]-vold>tolerance || v[1][0]-vold<-tolerance)	/* If the value of v at any meshpoint changes by more than tolerance, then iterate again */
			intol=0;

		
		/* Now treat rest of m=1, using the capillary BC */
		for(n=1;n<nn-1;n++)
		{
			vold=v[1][n];
			v[1][n]=1/(k19+h()*k20) * (k19*v[3][n]+h()*(vegfold[0][n-1]+vegfold[0][n])/2);	/* Use vegfold because capillary values have already been updated. Take average because there is no meshpoint in the capillary directly below */
			if(v[1][n]-vold>tolerance || v[1][n]-vold<-tolerance)
				intol=0;
		}
		
		/* And finally treat m=1, n=nn-1 (final meshpoint on m=1 row) using BC at x=1 */
		vold=v[1][nn-1];
		v[1][nn-1]=v[1][nn-2];
		if(v[1][nn-1]-vold>tolerance || v[1][nn-1]-vold<-tolerance)
			intol=0;
		
		

	}
	
	/* Cycle through VEGF meshpoints and set VEGF at time step j+1 */
	for(m=1;m<lvv;m++)
	{
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
			{
				vegfold[m][n]=vegf[m][n];
				vegf[m][n]=v[m][n];
				if(vegf[m][n]<(double)0)
				{
					vegf[m][n]=(double)0;
					vegfwarn++;
				}
			}
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
			{
				vegfold[m][n]=vegf[m][n];
				vegf[m][n]=v[m][n];
				if(vegf[m][n]<(double)0)
				{
					vegf[m][n]=(double)0;
					vegfwarn++;
				}
			}
		}
	}
	/* printf("GS iterations %4d",count); */
}






/* GAUSS-SIEDEL ITERATION: CALCULATES FIBRONECTIN VALUES AT TIME STEP j+1 */
void updatefib()
{
	double f[lvv][nn];
	double fold;
	double p;
	int n,m;
	int intol=0;
	int count=0;


	/* Cycle through substrate meshpoints */
	/* First do m=0, capillary. m is even so number of substrate meshpoint in x-direction is nn-1 */
	for(n=0;n<nn-1;n++)
	{
		p=(double)denscale*(occold[0][n]+occold[0][n+1])/2;		/* Average density at cell meshpoints to the right and left of substrate meshpoint at time j */
		fibold[0][n]=fib[0][n];
		fib[0][n]=fib[0][n] + k()*(k4*fib[0][n]*(1-fib[0][n])*p - k5*pro[0][n]*fib[0][n]/(1+k6*fib[0][n]));
		if(fib[0][n]<(double)0)
		{
			fib[0][n]=(double)0;
			fibwarn++;
		}
		if(fib[0][n]>(double)1)
		{
			fib[0][n]=(double)1;
			fibwarn++;
		}
	}

	
	
	/* Cycle through substrate meshpoints and set initial guess for f to be value for previous time step (j) */
	for(m=1;m<lvv;m++)
	{
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
				f[m][n]=fib[m][n];
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
				f[m][n]=fib[m][n];
		}
	}


	while(intol==0)	/* Keep iterating until the value of v at each meshpoint changes by less than a certain tolerance */
	{

		intol=1;
		count++;
		
		/* Cycle through substrate meshpoints and make improved estimate for f - treating boundary meshpoints separately - and checking tolerance at each update */
		/* NOTE: Rows are treated in reverse order (and comments read in reverse order) */
		
		/* Same routine for m=lvv-1 */
		fold=f[lvv-1][0];
		f[lvv-1][0]=f[lvv-1][1];
		if(f[lvv-1][0]-fold>tolerance || f[lvv-1][0]-fold<-tolerance)
			intol=0;
		for(n=1;n<nn-2;n++)
		{
			fold=f[lvv-1][n];
			f[lvv-1][n]=f[lvv-3][n];
			if(f[lvv-1][n]-fold>tolerance || f[lvv-1][n]-fold<-tolerance)
				intol=0;
		}
		fold=f[lvv-1][nn-2];
		f[lvv-1][nn-2]=f[lvv-1][nn-3];
		if(f[lvv-1][nn-2]-fold>tolerance || f[lvv-1][nn-2]-fold<-tolerance)
			intol=0;

		/* And last of all treat final two rows, m=lvv-2 and m=lvv-1, in similar way to first two (lvv is odd so lvv-2 is odd, lvv-1 is even) */
		fold=f[lvv-2][0];
		f[lvv-2][0]=f[lvv-2][1];
		if(f[lvv-2][0]-fold>tolerance || f[lvv-2][0]-fold<-tolerance)
			intol=0;
		for(n=1;n<nn-1;n++)
		{
			fold=f[lvv-2][n];
			f[lvv-2][n]=f[lvv-4][n];
			if(f[lvv-2][n]-fold>tolerance || f[lvv-2][n]-fold<-tolerance)
				intol=0;
		}
		fold=f[lvv-2][nn-1];
		f[lvv-2][nn-1]=f[lvv-2][nn-2];
		if(f[lvv-2][nn-1]-fold>tolerance || f[lvv-2][nn-1]-fold<-tolerance)
			intol=0;

		/* Now cycle through interior rows */
		for(m=lvv-3;m>2;m--)
		{
			/* First meshpoint on row requires BC at x=0 */
			fold=f[m][0];
			f[m][0]=f[m][1];
			if(f[m][0]-fold>tolerance || f[m][0]-fold<-tolerance)
				intol=0;
		
			if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 (but don't include last meshpoint in loop - boundary element) */
			{
				for(n=1;n<nn-2;n++)
				{
					fold=f[m][n];
					f[m][n]=relax2*k()/(h()*h()+2*k9*k()) * (0.5*k9*(f[m][n+1]+f[m][n-1]+f[m+2][n]+f[m-2][n]+fib[m][n+1]+fib[m][n-1]+fib[m+2][n]+fib[m-2][n]) + (h()*h()/k()-2*k9+k10*h()*h()*(1-fib[m][n])-k5*h()*h()*pro[m][n]/(1+k6*fib[m][n]))*fib[m][n]) + (1-relax2)*f[m][n];
					if(f[m][n]-fold>tolerance || f[m][n]-fold<-tolerance)
						intol=0;
				}
				fold=f[m][nn-2];
				f[m][nn-2]=f[m][nn-3];
				if(f[m][nn-2]-fold>tolerance || f[m][nn-2]-fold<-tolerance)
					intol=0;
			}
			
			else		/* If m is odd, number of subsrate meshpoints in x-direction is nn (but don't include last meshpoint in loop) */
			{
				for(n=1;n<nn-1;n++)
				{
					fold=f[m][n];
					f[m][n]=relax2*k()/(h()*h()+2*k9*k()) * (0.5*k9*(f[m][n+1]+f[m][n-1]+f[m+2][n]+f[m-2][n]+fib[m][n+1]+fib[m][n-1]+fib[m+2][n]+fib[m-2][n]) + (h()*h()/k()-2*k9+k10*h()*h()*(1-fib[m][n])-k5*h()*h()*pro[m][n]/(1+k6*fib[m][n]))*fib[m][n]) + (1-relax2)*f[m][n];
					if(f[m][n]-fold>tolerance || f[m][n]-fold<-tolerance)
						intol=0;
				}
				fold=f[m][nn-1];
				f[m][nn-1]=f[m][nn-2];
				if(f[m][nn-1]-fold>tolerance || f[m][nn-1]-fold<-tolerance)
					intol=0;
			}
		}


		/* Same routine for m=2, which is also classed as a boundary row */
		fold=f[2][0];
		f[2][0]=f[2][1];
		if(f[2][0]-fold>tolerance || f[2][0]-fold<-tolerance)
			intol=0;
		for(n=1;n<nn-2;n++)
		{
			fold=f[2][n];
			f[2][n]=f[4][n];
			if(f[2][n]-fold>tolerance || f[2][n]-fold<-tolerance)
				intol=0;
		}
		fold=f[2][nn-2];
		f[2][nn-2]=f[2][nn-3];
		if(f[2][nn-2]-fold>tolerance || f[2][nn-2]-fold<-tolerance)
			intol=0;
		
		
		/* First treat m=1 (capillary wall), n=0 using BC at x=0 (simpler condition than capillary BC) */
		fold=f[1][0];
		f[1][0]=f[1][1];
		if(f[1][0]-fold>tolerance || f[1][0]-fold<-tolerance)	/* If the value of f at any meshpoint changes by more than tolerance, then iterate again */
			intol=0;

		
		/* Now treat rest of m=1, using the capillary BC */
		for(n=1;n<nn-1;n++)
		{
			fold=f[1][n];
			f[1][n]=f[3][n];
			if(f[1][n]-fold>tolerance || f[1][n]-fold<-tolerance)
				intol=0;
		}
		
		/* And finally treat m=1, n=nn-1 (final meshpoint on m=1 row) using BC at x=1 */
		fold=f[1][nn-1];
		f[1][nn-1]=f[1][nn-2];
		if(f[1][nn-1]-fold>tolerance || f[1][nn-1]-fold<-tolerance)
			intol=0;
		
		

	}
	
	/* Cycle through substrate meshpoints and set fibronectin at time step j+1 */
	for(m=1;m<lvv;m++)
	{
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
			{
				fibold[m][n]=fib[m][n];
				fib[m][n]=f[m][n];
				if(fib[m][n]<(double)0)
				{
					fib[m][n]=(double)0;
					fibwarn++;
				}
				if(fib[m][n]>(double)1)
				{
					fib[m][n]=(double)1;
					fibwarn++;
				}
			}
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
			{
				fibold[m][n]=fib[m][n];
				fib[m][n]=f[m][n];
				if(fib[m][n]<(double)0)
				{
					fib[m][n]=(double)0;
					fibwarn++;
				}
				if(fib[m][n]>(double)1)
				{
					fib[m][n]=(double)1;
					fibwarn++;
				}
			}
		}
	}
	/* printf(" and %4d\n",count); */
}



/* UPDATE PROTEASE LEVELS AT TIMESTEP j */
void updatepro(int j)
{
	int n,m;
	double p;
	double totprocap=0.0;
	FILE *outp;

	/* Cycle through substrate meshpoints */
	/* First do m=0, capillary. m is even so number of substrate meshpoint in x-direction is nn-1 */
	for(n=0;n<nn-1;n++)
	{
		p=(double)denscale*(occold[0][n]+occold[0][n+1])/2;		/* Average density at cell meshpoints to the right and left of substrate meshpoint at time j */
		proold[0][n]=pro[0][n];
		pro[0][n]=pro[0][n] + k()*(profact*k1*vegfold[0][n]*p/(vegfold[0][n]+1) - k3*pro[0][n]);	/* Use vegfold because vegf values have already been updated (as have fib values) */
		if(pro[0][n]<(double)0)
		{
			pro[0][n]=(double)0;
			prowarn++;
		}
		totprocap=totprocap+pro[0][n];
	}
	totprocap=totprocap/nn;
	
	/* Also treat m=1 seperately - shouldn't take account of cell denisty in the capillary, so just use occold[1][n], not average of occold[0][n] and occold[1][n] */
	for(n=0;n<nn;n++)
	{	p=(double)pow(denscale,2)*occold[1][n];
		proold[1][n]=pro[1][n];
		pro[1][n]=pro[1][n] + k()*(profact*k1*vegfold[1][n]*p/(vegfold[1][n]+1) - k3*pro[1][n]);
		if(pro[1][n]<(double)0)
		{
			pro[1][n]=(double)0;
			prowarn++;
		}
	}

		
	for(m=2;m<lvv;m++)
	{
		if(m%2==0)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			for(n=0;n<nn-1;n++)
			{
				p=(double)pow(denscale,2)*(occold[m/2][n]+occold[m/2][n+1])/2;		/* Average density at cell meshpoints to the right and left of substrate meshpoint at time j */
				proold[m][n]=pro[m][n];
				pro[m][n]=pro[m][n] + k()*(profact*k1*vegfold[m][n]*p/(vegfold[m][n]+1) - k3*pro[m][n]);
				if(pro[m][n]<(double)0)
				{
					pro[m][n]=(double)0;
					prowarn++;
				}
			}
		}
		
		else		/* If m is odd, number of substrate meshpoints in x-direction is nn */
		{
			for(n=0;n<nn;n++)
			{
				p=(double)pow(denscale,2)*(occold[(m-1)/2][n]+occold[(m+1)/2][n])/2;	/* Average density at cell meshpoints above and below substrate meshpoint at time j */
				proold[m][n]=pro[m][n];
				pro[m][n]=pro[m][n] + k()*(profact*k1*vegfold[m][n]*p/(vegfold[m][n]+1) - k3*pro[m][n]);
				if(pro[m][n]<(double)0)
				{
					pro[m][n]=(double)0;
					prowarn++;
				}
			}
		}
	}
	
	if(j%10==0)
	{
		if((outp=fopen("totprocap.dat","a"))==NULL)
			printf("Cannot open output file totprocap.dat\n");		
		else
		{	
			fprintf(outp,"j=%i  totprocap=%9.6f\n",j,totprocap);
			fclose(outp);
		}
	}
}






/* WRITE PARAMETERS TO FILE */
void writeparams(int cellsinsim,int starttime,int run)
{
	int i;
	double avginv=(double)0;
	int cellsincap=0;
	FILE *outp;
	int mins,secs;
	int endtime=(int)(unsigned)time(NULL);
	
	/* Count number of cells in capillary and calculate average invasion distance of live cells  */
	for(i=0;i<ii;i++)
	{
		if(deathtime[i]==jj-1 && vpos[i][jj-1]==0)
			cellsincap++;
		avginv=avginv+cellxy(vpos[i][deathtime[i]]);
	}
	avginv=avginv/ii;
	
	mins=(endtime-starttime)/60;
	secs=(endtime-starttime)%60;
	
	
	/* Write parameters */
	if((outp=fopen("params.dat","a"))==NULL)
	{
		printf("Cannot open output file params.dat\n");
	}
	else
	{
		if(run==0)
		{
			fprintf(outp,"fi=%9.6f\n",fi);
			fprintf(outp,"f1=%9.6f\n",f1);
			fprintf(outp,"gamma1=%9.6f\n",gamma1);
			fprintf(outp,"gamma2=%9.6f\n",gamma2);
			fprintf(outp,"gamma3=%9.6f\n",gamma3);
			fprintf(outp,"gamma4=%9.6f\n",gamma4);
			fprintf(outp,"gamma5=%9.6f\n",gamma5);
			fprintf(outp,"gamma6=%9.6f\n",gamma6);
			fprintf(outp,"v1=%9.6f\n",v1);
			fprintf(outp,"v2=%9.6f\n",v2);
			fprintf(outp,"v3=%9.6f\n",v3);
			fprintf(outp,"v4=%9.6f\n",v4);
			fprintf(outp,"m0=%9.6f\n",m0);	
			fprintf(outp,"m1=%9.6f\n",m1);
			fprintf(outp,"a1=%9.6f\n",a1);
			fprintf(outp,"k1=%9.6f\n",k1);
			fprintf(outp,"k2=%9.6f\n",k2);
			fprintf(outp,"k3=%9.6f\n",k3);
			fprintf(outp,"k4=%9.6f\n",k4);
			fprintf(outp,"k5=%9.6f\n",k5);
			fprintf(outp,"k6=%9.6f\n",k6);
			fprintf(outp,"k7=%9.6f\n",k7);
			fprintf(outp,"k8=%9.6f\n",k8);
			fprintf(outp,"k9=%9.6f\n",k9);
			fprintf(outp,"k10=%9.6f\n",k10);
			fprintf(outp,"k11=%9.6f\n",k11);
			fprintf(outp,"k12=%9.6f\n",k12);
			fprintf(outp,"k13=%9.6f\n",k13);
			fprintf(outp,"k14=%9.6f\n",k14);
			fprintf(outp,"k15=%9.6f\n",k15);
			fprintf(outp,"k16=%9.6f\n",k16);
			fprintf(outp,"k17=%9.6f\n",k17);
			fprintf(outp,"k18=%9.6f\n",k18);
			fprintf(outp,"k19=%9.6f\n",k19);
			fprintf(outp,"k20=%9.6f\n",k20);
			fprintf(outp,"k21=%9.6f\n",k21);
			fprintf(outp,"k22=%9.6f\n",k22);
			fprintf(outp,"k23=%9.6f\n",k23);
			fprintf(outp,"k24=%9.6f\n",k24);
			fprintf(outp,"k25=%9.6f\n",k25);
			fprintf(outp,"k26=%9.6f\n",k26);
			fprintf(outp,"pmax=%9.6f\n",pmax);
			fprintf(outp,"Q=%9.6f\n",q);
			fprintf(outp,"profact=%9.6f\n",profact);
			fprintf(outp,"child=%9.6f\n",child);
			fprintf(outp,"\nUpper bound for number of cells         uu = %i\n",uu);
			fprintf(outp,"Dmnls length of time                length = %12.9f\n",length);
			fprintf(outp,"Actual length of time               length = %12.9fh\n",length*hours);
			fprintf(outp,"Number of timesteps                     jj = %4d\n",jj);
			fprintf(outp,"Timestep length                          k = %15.12f\n",k());
			fprintf(outp,"Number of cell meshpoints               nn = %4d\n",nn);
			fprintf(outp,"Mesh size                                h = %9.6f\n",h());
			fprintf(outp,"Number of substrate meshpoints          vv = %4d\n",vv);
			fprintf(outp,"No of cell meshpoints used in y        lnn = %4d\n",lnn);
			fprintf(outp,"No of substrate meshpoints used in y   lvv = %4d\n",lvv);
			fprintf(outp,"GS tolerance                     tolerance = %15.12f\n",tolerance);
			fprintf(outp,"Relaxation parameter VEGF            omega = %9.6f\n",relax);
			fprintf(outp,"Relaxation parameter fibronectin    omega2 = %9.6f\n",relax2);
			fprintf(outp,"Waiting time parameter              lambda = %9.6f\n",lambda());
			fprintf(outp,"Prob of staying still (cap)   1-2*lambda*k = %9.6f\n",1-2*lambda()*k());
			fprintf(outp,"Prob of staying still (ecm)   1-4*lambda*k = %9.6f\n\n\n",1-4*lambda()*k());
		}
		
		fprintf(outp,"Run %i\n",run);
		fprintf(outp,"======\n\n");
		fprintf(outp,"Number of cells                         ii = %4d\n",ii);
		fprintf(outp,"Density scaling                   denscale = %9.6f\n",denscale);
		fprintf(outp,"Cells remaining in simulation        insim = %4d\n",cellsinsim);
		fprintf(outp,"Cells remaining in capillary         incap = %4d\n",cellsincap);
		fprintf(outp,"VEGF warnings                     vegfwarn = %i\n",vegfwarn);
		fprintf(outp,"Fibronectin warnings               fibwarn = %i\n",fibwarn);
		fprintf(outp,"Protease warnings                  prowarn = %i\n\n",prowarn);
		fprintf(outp,"Average vertical invasion distance  avginv = %9.6f\n",avginv);
		fprintf(outp,"Capillary moves in                   capin = %i\n",capin);
		fprintf(outp,"Capillary moves out                 capout = %i\n",capout);
		fprintf(outp,"Moves in                                in = %i\n",in);
		fprintf(outp,"Moves out                              out = %i\n",out);
		fprintf(outp,"Moves up                                up = %i\n",up);
		fprintf(outp,"Moves down                            down = %i\n",down);
		fprintf(outp,"Total moves                          moves = %i\n",capin+capout+in+out+up+down);
		fprintf(outp,"Births                              births = %i\n",births);
		fprintf(outp,"Deaths                              deaths = %i\n",deaths);
		fprintf(outp,"Running time: %i minutes and %i seconds\n",mins,secs);
		fprintf(outp,"Cumulative running/real time ratio: %4.1f\%\n\n",100*mins/(60*hours*length*(run+1)));
		fclose(outp);
	}
	

	/* Write parameters */
	if((outp=fopen("max.dat","w"))==NULL)
	{
		printf("Cannot open output file max.dat\n");
	}
	else
	{
		fprintf(outp,"%i\n",maxcell);
		for(i=0;i<6;i++)
			fprintf(outp,"%9.6f\n",maxsub[i]);
		fclose(outp);
	}
	
}

/* WRITE CELLS.PRO, TRAILS DATA */
void writetrails()
{
	char psname[10];
	char *f;
	int i,j,a,b,jp1;
	int move,lasthpos,lastvpos,npoints;
	int stopmove[6],file;
	int count=0;
	FILE *outp;
	int n[jj],m[jj];
	
	f=psname;
	
	/* Write trails */
	
	strcpy(psname,"cells");

	/* Write cell position data as IDL file */
	if((outp=fopen("cells.pro","w"))==NULL)
		printf("Cannot open output file cells.pro\n");
	else
	{
		npoints=jj*(1-pstill(1))*1.3;
		trailsintro(outp,npoints,f);

		for(i=0;i<ii;i++)
		{
			for(count=1;count<5;count++)
				stopmove[count]=npoints-1;
			count=1;
			move=1;
			n[0]=hpos[i][birthtime[i]];
			m[0]=vpos[i][birthtime[i]];
			
			for(j=birthtime[i]+1;j<=deathtime[i];j++)
			{
				jp1=min(j+1,jj-1);
				if((hpos[i][j]!=n[move-1] || vpos[i][j]!=m[move-1]) && vpos[i][jp1]>0)		/* Only write position data once out of capillary and each time the cell moves */
				{
					n[move]=hpos[i][j];
					m[move]=vpos[i][j];
					/* Record move number for which the cell exceeds count*jj/5 (to enable trails only upto certains times) (while instead of if allows for multiple 'snapshot points' to be passed on a single move) */
					while(j>count*jj/5)
					{
						stopmove[count]=move;
						count++;
					}
					
					a=move-1;
					while(a!=move && a>0)	/* Cycle through previous movements to check if an oxbow lake has been formed */
					{
						a--;
						if(n[move]==n[a] && m[move]==m[a])	/* If so, set the move counter back to the point where the oxbow started and thus subsequent data will overwrite the oxbow data */
						{
							move=a;
							for(b=1;b<5;b++)		/* Also pull stopmoves back to a+1 so they are stopped at the oxbow point */
							{
								if(stopmove[b]>a+1) stopmove[b]=a+1;
							}
					
						}
					}
					
					
					move++;
				}
			}
			
			for(a=0;a<move;a++)
			{
				fprintf(outp,"x(%d,%d)=%12.9f & y(%d,%d)=%12.9f\n",a,i,cellxy(n[a]),a,i,cellxy(m[a]));
			}				
					
       			if(move>npoints)
      				printf("WARNING: excess cell trails data  %i moves    npoints=%i\n",move,npoints);
      			else if(move<npoints)	/* Fill remaining array space with the final position of the cell */
      				fprintf(outp,"x(%i:%i,%i)=%12.9f & y(%i:%i,%i)=%12.9f\n",move,npoints-1,i,cellxy(n[move-1]),move,npoints-1,i,cellxy(m[move-1]));
  			for(file=2;file<6;file++)	/* For partial trails, fill unwanted array space with last wanted position */
  			{
  				if(stopmove[file-1]==0)
  				{
	 				fprintf(outp,"if (file eq %i) then begin\n",file);
 					fprintf(outp,"x(0:%i,%i)=0.0 & y(0:%i,%i)=0.0\n",npoints-1,i,npoints-1,i);
   					fprintf(outp,"endif\n");
   				}
 				else
 				if(stopmove[file-1]!=npoints-1)
 				{
	  				fprintf(outp,"if (file eq %i) then begin\n",file);
 					fprintf(outp,"x(%i:%i,%i)=x(%i,%i) & y(%i:%i,%i)=y(%i,%i)\n",stopmove[file-1],npoints-1,i,stopmove[file-1]-1,i,stopmove[file-1],npoints-1,i,stopmove[file-1]-1,i);
  					fprintf(outp,"endif\n");
  				}
  			}
  				
  							
    			fprintf(outp,"\n");
    		}
    	  		
    		trailsend(outp,(float)jj*k());
		fclose(outp);
	}

	/* For old format output including oxbow lakes, just remove the oxbow lake ckecking lines */


}

/* WRITE CELL CAPILLARY SNAPSHOT */
void writecellscap(int j)
{
	char path[25],psname[13],proend[6],psend[5];
	char *f;
	int *celldata;
	int count=0;
	FILE *outp;
	
	f=psname;

	/* Set celldata pointer to appropriate substrate variable */
	celldata=occ;


	/* Create IDL path+filename and postscipt filename */
	strcpy(path,"cellscap");
	strcpy(psname,"cellscap");
	count=(5*j)/jj+1;
	if(j==jj-2) count=6;
	switch(count)
	{
		case 1:
			strcpy(proend,"1.pro");
			strcpy(psend,"1.ps");
			break;
		case 2:
			strcpy(proend,"2.pro");
			strcpy(psend,"2.ps");
			break;
		case 3:
			strcpy(proend,"3.pro");
			strcpy(psend,"3.ps");
			break;
		case 4:
			strcpy(proend,"4.pro");
			strcpy(psend,"4.ps");
			break;
		case 5:
			strcpy(proend,"5.pro");
			strcpy(psend,"5.ps");
			break;
		default:
			strcpy(proend,"6.pro");
			strcpy(psend,"6.ps");
			break;
	}
		
	strcat(path,proend);
	strcat(psname,psend);
	
	/* Write data */
	if((outp=fopen(path,"w"))==NULL)
		printf("Cannot open output file %s\n",path);		
	else
	{
		cellscap(celldata,outp,j,f);
		fclose(outp);
	}
}



/* WRITE SUBSTRATE RESULTS TO FILE */
void writesub(int j)
{
	char path[22],pathcap[25],psname[10],psnamecap[13],proend[6],psend[5];
	char *f,*fcap;
	double *data;
	int type;
	int count=0;
	FILE *outp;
	
	f=psname;
	fcap=psnamecap;



	for(type=1;type<4;type++)	/* Type1=VEGF, Type2=Fibronectin, Type3=Protease */
	{
	
		/* Set data pointer to appropriate substrate variable */
		switch(type)
		{
			case 1:
				data=vegf;
				break;
			case 2:
				data=fib;
				break;
			default:
				data=pro;
				break;
		}


		/* Create substrate IDL path+filename and postscipt filename */
		
		switch(type)
		{
			case 1:
				strcpy(path,"vegf");
				break;
			case 2:
				strcpy(path,"fib");
				break;
			default:
				strcpy(path,"pro");
				break;
		}
		strcpy(psname,path);
		
		/* Also create capillary filenames */
		strcpy(pathcap,path);
		strcat(pathcap,"cap");
		strcpy(psnamecap,psname);
		strcat(psnamecap,"cap");					
		
		count=(5*j)/jj+1;
		if(j==jj-2) count=6;
		switch(count)
		{
			case 1:
				strcpy(proend,"1.pro");
				strcpy(psend,"1.ps");
				break;
			case 2:
				strcpy(proend,"2.pro");
				strcpy(psend,"2.ps");
				break;
			case 3:
				strcpy(proend,"3.pro");
				strcpy(psend,"3.ps");
				break;
			case 4:
				strcpy(proend,"4.pro");
				strcpy(psend,"4.ps");
				break;
			case 5:
				strcpy(proend,"5.pro");
				strcpy(psend,"5.ps");
				break;
			default:
				strcpy(proend,"6.pro");
				strcpy(psend,"6.ps");
				break;
		}
		
		strcat(path,proend);
		strcat(psname,psend);
		strcat(pathcap,proend);
		strcat(psnamecap,psend);
		
		if(breach || type==1)	/* Only write ECM after cells break through (except for VEGF) */
		{
			/* Write data in ECM  */
			if((outp=fopen(path,"w"))==NULL)
				printf("Cannot open output file %s\n",path);		
			else
			{			
				subecm(type,data,outp,j,f);
				fclose(outp);
			}
		}
		
		/* Write data in capillary */
		if((outp=fopen(pathcap,"w"))==NULL)
			printf("Cannot open output file %s\n",pathcap);		
		else
		{
			subcap(type,data,outp,j,fcap);
			fclose(outp);
		}

	}

}


/* CREATE IDL BATCH FILES AND WRITE RUNNING TIME */
void writebatch()
{
	char path[22],proend[6];
	int i,type;
	FILE *outp;
	
	if((outp=fopen("runcaps.pro","w"))==NULL)
		printf("Cannot open output file runcaps.pro\n");		
	else
	{	
		for(i=1;i<7;i++)
		{
			fprintf(outp,".run cellscap%i.pro\n",i);
			fprintf(outp,"$ rm cellscap%i.pro\n",i);
			fprintf(outp,".run vegfcap%i.pro\n",i);
			fprintf(outp,"$ rm vegfcap%i.pro\n",i);
			fprintf(outp,".run fibcap%i.pro\n",i);
			fprintf(outp,"$ rm fibcap%i.pro\n",i);
			fprintf(outp,".run procap%i.pro\n",i);
			fprintf(outp,"$ rm procap%i.pro\n",i);
		}
		fclose(outp);
	}


	for(type=1;type<4;type++)
	{
		switch(type)
		{
			case 1:
				strcpy(proend,"vegf");
				break;
			case 2:
				strcpy(proend,"fib");
				break;
			default:
				strcpy(proend,"pro");
				break;
		}
			
		strcpy(path,"run");
		strcat(path,proend);
		strcat(path,".pro");
		
		if((outp=fopen(path,"w"))==NULL)
			printf("Cannot open output file %s\n",path);		
		else
		{	
			for(i=1;i<7;i++)
			{		
				fprintf(outp,".run %s%i.pro\n",proend,i);
				fprintf(outp,"$ rm %s%i.pro\n",proend,i);
			}
			fclose(outp);
		}
	}
}












/* PROBABILITY OF STAYING STILL */
double pstill(int m)
{
	double prob;
	
	/* Assume waiting time is constant */
	
	if(m==0)
		prob=(double)1-2*lambda()*k();
	else
		prob=(double)1-4*lambda()*k();
	
	return(prob);
}



/* PROBABILITY OF MOVING FROM CELL MESHPOINT (n,m) ACROSS BARRIER c (0=left, 1=right, 2=down, 3=up) */
double pmove(int n,int m,int c)
{
	double prob;
	double bar[4];
	
	if(n>0)		/* If possible to move left */
		bar[0]=tau(pro[2*m][n-1],fib[2*m][n-1],vegf[2*m][n-1],m);	/* Tau to the left */
	else
		bar[0]=0.0;
		
	if(n<nn-1)	/* If possible to move right */
		bar[1]=tau(pro[2*m][n],fib[2*m][n],vegf[2*m][n],m);		/* Tau to the right */
	else
		bar[1]=0.0;
		
	if(m>1)		/* If possible to move down (assume cannot move back into capillary from ECM) */
		bar[2]=tau(pro[2*m-1][n],fib[2*m-1][n],vegf[2*m-1][n],m);	/* Tau below */
	else
		bar[2]=0.0;
	
	if(m<lnn-1 && m>0)	/* If possible to move up (assume cannot move out of capillary - this happens by a separate mechanism) */
		bar[3]=tau(pro[2*m+1][n],fib[2*m+1][n],vegf[2*m+1][n],m);	/* Tau above */
	else
		bar[3]=0.0;
		
	if(m==0)
		prob=(double)2*lambda()*k()*bar[c]/(bar[0]+bar[1]);
	else
		prob=(double)4*lambda()*k()*bar[c]/(bar[0]+bar[1]+bar[2]+bar[3]);
		
	return(prob);
}


/* TAU FUNCTION OF PROTEASE AND FIBROECTIN */
double tau(double c,double f,double v,int m)
{
	if(m==0)	/* Cell is in capillary */
		return(pow((c+k11)/(c+k12),gamma1) * pow((f+k13)/(f+k14),gamma2) * pow((v+v1)/(v+v2),gamma5));
	else		/* Cell is in ECM */
		return(pow((c+k15)/(c+k16),gamma3) * pow((f+k17)/(f+k18),gamma4) * pow((v+v3)/(v+v4),gamma6));
}



/* CALCULATE SPACE STEP LENGTH h */
double h()
{
	return((double)1/(nn-1));
}

/* CALCULATE TIME STEP LENGTH k */
double k()
{
	return((double)length/jj);
}


/* CALCULATE PARAMETER lambda */
double lambda()
{
	/* Lambda=D/h^2 */
	return(1/(h()*h()));
}


void trailsintro(FILE *outp,int a,char *f)
{	
	fprintf(outp,"set_plot, 'ps'\n");
	fprintf(outp,"!P.MULTI = [0,1,1]\n");
	fprintf(outp,"npoints = %d\n",a);
	fprintf(outp,"ntracks = %d\n",ii);
	fprintf(outp,"x = fltarr(npoints,ntracks) & y = x\n");
	fprintf(outp,"file=0\n");
	fprintf(outp,"read,file,prompt='Enter file number (2-6) '\n");
	fprintf(outp,"fname=strcompress('%s'+strlowcase(file),/REMOVE_ALL)+'.ps'\n",f);
	fprintf(outp,"device, /portrait, filename = fname, $\n");
	fprintf(outp,"ysize=20, yoff=5\n");

	fprintf(outp,"\n");						
}


void trailsend(FILE *outp,float t)
{
	fprintf(outp,"\n");
	fprintf(outp,"!p.psym = -3\n");
	fprintf(outp,"!p.thick = 1.5\n");
	fprintf(outp,"!p.charsize = 2.0\n");
	fprintf(outp,"\n");
	fprintf(outp,"for k = 0, ntracks-1 do begin\n");
	fprintf(outp,"hours=%6.3f*(file-1)/5\n",t*hours);
	fprintf(outp,"ttl=strcompress('t='+strlowcase(hours),/REMOVE_ALL)+'h'\n");
	fprintf(outp,"plot, x(*,k), y(*,k), $\n");
	fprintf(outp,"Title=ttl,$\n");
	fprintf(outp,"XTITLE = 'x', YTITLE = 'y', xrange = [0,1], $\n");
	fprintf(outp,"yrange = [%12.9f,0.5], /noerase\n",h()/2);
	fprintf(outp,";device, /close & spawn, 'ghostview fname & \n"); 
	fprintf(outp,"endfor\n");
	fprintf(outp,"set_plot, 'X'\n");
	fprintf(outp,"END\n");	
}


/* INITIAL VEGF FUNCTION */
double vegf0(double x,double y)
{
	return((double)0);
}

/* X- OR Y-COORDINATE FOR CELL MESHPOINT (N,?) OR (?,M) */
double cellxy(int n)
{
	return((double)n/(nn-1));
}

/* Y-COORDINATE FOR VEGF MESHPOINT (?,M) */
double vy(int m)
{
	return((double)m/(2*(nn-1)));
}

/* X-COORDINATE FOR VEGF MESHPOINT (N,M) */
double vx(int n,int m)
{
	if(m%2==0)
		return((double)(n+0.5)/(nn-1));
	else
		return((double)n/(nn-1));
}

/* INITIAL FIBRONECTIN FUNCTION */
double fib0(double x,double y)
{
	if(y==0)
		return((double)fi);
	else
		return((double)fi);
}

/* INITIAL PROTEASE FUNCTION */
double pro0(double x,double y)
{
	return((double)0);
}

/* HEAVISIDE FUNCTION */
int heaviside(double x)
{
	if(x>=0)
		return(1);
	else
		return(0);
}

/* WRITE SUBSTRATE ECM DATA FOR A GIVEN TIMESTEP */
void subecm(int type,double *data,FILE *outp,int j,char *f)
{
	char label[20];
	int n,m,nsize,msize;
	
	nsize=(int)((nn-2)/subinc)+1;
	msize=(int)((lvv-1)/(2*subinc))+1;

	fprintf(outp,"openr,1,'max.dat'\n");
	fprintf(outp,"for k=0,%i do begin\n",type);	/* type here and type-1 below because maxcell comes first in max.dat */
	fprintf(outp,"readf,1,max\n");
	fprintf(outp,"endfor\n");
	fprintf(outp,"close,1\n");
	fprintf(outp,"Set_Plot, 'PS'\n");
	fprintf(outp,"Device, Filename=\"%s\"\n",f);
	fprintf(outp,"z=fltarr(%d,%d)\n",nsize,msize);
	fprintf(outp,"\n");
	
	for(m=2*subinc;m<lvv;m=m+2*subinc)	/* Write data for m even only, and not including capillary data */
	{
		for(n=0;n<nn-1;n=n+subinc)	/* If m is even, number of substrate meshpoints in x-direction is nn-1 */
		{
			fprintf(outp,"z(%d,%d)=%12.9f\n",n/subinc,m/(2*subinc),*(data + m*nn + n));	/* Pointer points to data[m][n] */
			if(*(data+m*nn+n)>maxsub[type-1])
				maxsub[type-1]=*(data+m*nn+n);
		}
	}
	
	fprintf(outp,"\n");
	fprintf(outp,"x=FIndGen(%d)/%d\n",nsize,(nn-1)/subinc-1);
	fprintf(outp,"y=FIndGen(%d)/%d\n",msize,(nn-1)/subinc);
	fprintf(outp,"surface,z,x,y,$\n");
	fprintf(outp,"Title='t=%6.3fh',$\n",(double)j*k()*hours);
	switch(type)
	{
		case 1:
			strcpy(label,"VEGF conc.");
			break;
		case 2:
			strcpy(label,"Fibronectin conc.");
			break;
		default:
			strcpy(label,"Protease conc.");
			break;
	}
	fprintf(outp,"XTitle='x', YTitle='y', ZTitle='%s',$\n",label);
	
	fprintf(outp,"ZRange=[0,max],$\n");		/* b is top end of z=range */
	fprintf(outp,"CharSize=2.0\n");
	fprintf(outp,"Device, /Close\n");
	fprintf(outp,"Set_Plot,'X'\n");
	fprintf(outp,"end\n");
}



/* WRITE SUBSTRATE CAPILLARY DATA FOR A GIVEN TIMESTEP */
void subcap(int type,double *data,FILE *outp,int j,char *f)
{
	char label[20];
	int n;

	fprintf(outp,"openr,1,'max.dat'\n");
	fprintf(outp,"for k=0,%i do begin\n",type+3);
	fprintf(outp,"readf,1,max\n");
	fprintf(outp,"endfor\n");
	fprintf(outp,"close,1\n");
	fprintf(outp,"Set_Plot, 'PS'\n");
	fprintf(outp,"Device, Filename=\"%s\"\n",f);
	fprintf(outp,"z=fltarr(%d)\n",nn-1);
	fprintf(outp,"\n");
	
	for(n=0;n<nn-1;n++)	/* m=0, even, so number of substrate meshpoints in x-direction is nn-1 */
	{
		fprintf(outp,"z(%d)=%12.9f\n",n,*(data + n));	/* Pointer points to data[0][n] */
		if(*(data+n)>maxsub[type+2])
			maxsub[type+2]=*(data+n);
	}
	fprintf(outp,"\n");
	fprintf(outp,"x=FIndGen(%d)/%d\n",nn-1,nn-2);
	switch(type)
	{
		case 1:
			strcpy(label,"VEGF conc.");
			break;
		case 2:
			strcpy(label,"Fibronectin conc.");
			break;
		default:
			strcpy(label,"Protease conc.");
			break;
	}

	fprintf(outp,"plot,x,z,Title='t=%6.3fh',XTitle='x', YTitle='%s',YRange=[0,max]\n",(double)j*k()*hours,label);	/* b is top end of y-range */
	fprintf(outp,"CharSize=2.0\n");
	fprintf(outp,"Device, /Close\n");
	fprintf(outp,"Set_Plot,'X'\n");
	fprintf(outp,"end\n");
}



/* WRITE CELL CAPILLARY DATA FOR A GIVEN TIMESTEP */
void cellscap(int *data,FILE *outp,int j,char *f)
{
	int n;
	int count=0;
	
	fprintf(outp,"openr,1,'max.dat'\n");
	fprintf(outp,"readf,1,max\n");
	fprintf(outp,"close,1\n");
	fprintf(outp,"Set_Plot, 'PS'\n");
	fprintf(outp,"Device, Filename=\"%s\"\n",f);
	fprintf(outp,"z=intarr(%d)\n",nn);
	fprintf(outp,"\n");
	
	for(n=0;n<nn;n++)	/* m=0, even, so number of substrate meshpoints in x-direction is nn-1 */
	{
		fprintf(outp,"z(%d)=%i\n",n,*(data + n));	/* Pointer points to data[0][n] */
		count = count + *(data + n);			/* Find total number of cells in capillary */
		if(*(data+n)>maxcell)
			maxcell=*(data+n);
		
	}
	fprintf(outp,"\n");
	fprintf(outp,"x=FIndGen(%d)/%d\n",nn,nn-1);
	fprintf(outp,"plot,x,z,Title='t=%6.3fh',XTitle='x', YTitle='Cells (total %i)',YRange=[0,max]\n",(double)j*k()*hours,count);
	fprintf(outp,"CharSize=2.0\n");
	fprintf(outp,"Device, /Close\n");
	fprintf(outp,"Set_Plot,'X'\n");
	fprintf(outp,"end\n");
}


/* RETURNS MAXIMUM OF TWO INTEGERS */
int max(int a,int b)
{
	if(a>b)
		return(a);
	else
		return(b);
}
	
/* RETURNS MINIMUM OF TWO INTEGERS */
int min(int a,int b)
{
	if(a<b)
		return(a);
	else
		return(b);
}


/* RETURNS 0 IF OPTION IS "n", 1 IF OPTION IS "t", AND 2 IF NEITHER (->ERROR) */
int optatoi(char *opt)
{
	if(strcmp(opt,"n")==0) return(0);
	else if(strcmp(opt,"t")==0) return(1);
	else return(2);
}

/* RAN0 improved random number generator: firstcall should have a pointer to a negative seed (-starttime) as argument; subsequent calls have the same pointer as argument (which will now point to the integer 1) */
float ran0(int *idum)
{
	static float y,maxran,v[98];
	float dum;
	static int iff=0;
	int j;
	unsigned i,k;
	
	if(*idum<0 || iff==0)
	{
		iff=1;
		i=2;
		maxran=RAND_MAX+1.0;
		srand(*idum);
		*idum=1;
		for(j=1;j<97;j++) dum=rand();
		for(j=1;j<97;j++) v[j]=rand();
		y=rand();
	}
	j=1+97.0*y/maxran;
	y=v[j];
	v[j]=rand();
	return y/maxran;
}
