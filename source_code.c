
//*********STRING OF CELLS EXCHANGING mtDNA with LINEAR FEEDBACK CONTROL FOR THE BIRTH RATE *****************
//Ferdinando Insalata
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define RND drand48() //random number generator, numbers in [0,1[

//functions, see each function for an explanation of what it does
double rep_rate(double w,double m, double delta);
double sum_array(double array[], int num_elements);
double partsum(double tosum[],int i);

//This code simulates the birth-deat process of mtDNA in a string of n cells that exchange mtDNA molecules as well.
//To change the number of cells it is necessary to change the value of the parameter n and to initialize the 
//population vector nitial_pop[4]={450.0,50.0,500.0,0} accordingly.

//The degradation rate \mu is constant, the birth rate is controlled by a LINEAR FEEDBACK CONTROL SYSTEM,
//implemented by the function "rep_rate(double w,double m, double delta)"
//The birth-death process is a Poisson point process simulated by means of the Gillespie algorithm 


//This version of the source code reproduces the result of FIg. S1.D of the manuscript entitled
// "Survival of the densest accounts for the expansion of mitochondrial mutations in ageing", specifically of the 
//second version of the manuscript submitted and uplodaded to bioRxiv in May 2021.
//The initial parameter configuration reproduces the points in the red line of panel S1.D, that allows the observation of Survival of the Densest 
//in the simplest possible systems. Green and blue can be reproduced
//by adding preferential elimination of mutants, as explained in section SI.8. 

//It might take some time (30 mins) to reproduce the results, since the final time of the simulation (final_t) is large.
// Setting the variable final_t to a smaller value (like 1000) should be enough to see the increase in mean heteroplasmy in only a few minutes.

int main()
{

/*definition of parameters */
int realizations; //number of samples
int nprint; // prints every nprint realizations
double final_t;
int n; /*number of cells*/
int M; /*number of mesaurements*/
double epsilon; /*selection against mutants, IN PERCENTAGE*/
double k; /*hopping rate between cells 1 and 2*/
double delta; //control on mutants
double mu; /*degradation rate of wildtypes*/


//initialization of parameters
n=2;  
M=30;
final_t=30000;
realizations=1000;
nprint= 100;
mu=0.07; //common degradation rate
k=mu; //hopping rate !!CHANGE ALSO IN THE FUNCTION THAT GIVES THE REPLICATION RATE!!
epsilon=0.0;  //selection against mutants, set to zero for Fig. S1.D red line. Add selection for green and blue line as explained in SI.8 of the paper.
delta =0.1; //pstrength of control on mutants. delta= 1 for same control on mutants and wildtypes
srand48(12127); //set the seed for the random number generator

FILE * fp;

FILE * prg;
//file to check progress

FILE * pop;


//***************************************************definition of variables **************************
int i,j,h,rl,z; //indexes
double xi; //random number
double wait_time; //time increment
double t; //time variable

double mum; /*degradation rate of mutants*/
int V[2*n][8*n-4]; //stoichiometry matrix S
double x[2*n]; //vector of copy numbers, wildtype and mutant
double ratev[8*n-4], rate; //vector of rates and its sum
double initial_pop[4]={455.0,0.0,450.0,50.0};
double temp_mean;

//create array of times
double times[M+1];
times[0]=0;
for( i = 1; i <= M; i = i + 1 ){
	times[i]=(final_t/M)*i;    
	}

//matrices to store data on population
double populations[M+1][2*n]; //matrix to record copy numbers of all 2n species,for the M+1 measurement times
double popsquared[M+1][2*n]; //matrix to record squares of copy numbers of all 2n species,  for the M+1 measurement times
double popvars[M+1][2*n];//variances of populations, for the M+1 measurement times
double popsigmas[M+1][2*n];//std of populations, for the M+1 measurement times

//matrices to store data on heteroplasmies
double het[M+1][n];//heteroplasmies, for the M+1 measurement times
double hetsquared[M+1][n]; //heteroplasmies squared, for the M+1 measurement times
double hvars[M+1][n]; //heteroplasmies variances, for the M+1 measurement times
double hsigmas[M+1][n]; //heteroplasmies std, for the M+1 measurement times

////////////////////////////////////////////////////////////////////////////////////////
////initialization of variables///////////////
temp_mean=0;
mum=mu*(1+epsilon);

for( i=0; i<M+1; i = i + 1 ){   
	for (j=0; j<n; j=j+1) {		
		het[i][j]=0;
		hetsquared[i][j]=0;
		hvars[i][j]=0;
		hsigmas[i][j]=0;
	}	
	for (j=0; j<2*n; j=j+1) {		
		populations[i][j]=0;
		popsquared[i][j]=0;
		popvars[i][j]=0;
		popsigmas[i][j]=0;
	}	
}

//here the stoichiometry matrix is built #################################################
for (i=0; i<2*n; i++) {
	for(j=0; j<8*n-4; j++) {
		V[i][j]=0;
	}
}
for (i=0; i<2*n; i++) {
	V[i][i]=-1;
	V[i][2*n+i]=1;
}
for (h=0; h<n-1; h++) {
	for (i=2*h; i<2*(h+2); i++) {
		for (j=4*(n+h);j<4*(n+h+1); j++) {
			if ((i+j)%2 == 0) V[i][j]=-1+abs(i-2*h-j+4*(h+n));			
		}
	}
}
//the stoichiometry matrix has been built ###############################

printf("\n");
printf("\n");
printf("The stoichiometry matrix is :");
printf("\n");

//print the stoichiometry matrix, just to check
for (i=0; i<2*n; i++) {
	for(j=0; j<8*n-4; j++) {
		printf("%i  ", V[i][j]); 
		
	}
	printf("\n");
}

if (n==1)
{
	k=0;
}
//////*********************************************************LOOP OVER REALIZATIONS **************************************************
for( rl = 1; rl <= realizations; rl++ )  {
	
	for (i=0;i<2*n;i++) 
	x[i]=initial_pop[i];


	for(i=0; i<2*n; i++) {
		populations[0][i]=populations[0][i]+x[i];
		popsquared[0][i]=popsquared[0][i]+(x[i])*x[i];
	}
	
	for(i=0; i<n; i++) {
		het[0][i]=het[0][i]+x[2*i+1]/(x[2*i]+x[2*i+1]);
		hetsquared[0][i] = hetsquared[0][i] +(x[2*i+1]/(x[2*i]+x[2*i+1]))*(x[2*i+1]/(x[2*i]+x[2*i+1]));
	}

	
	//////**************************************LOOP OVER TIME ******************************
	t=0;
	j=1; //put it to 1 because j=0 corresponds to time zero
	while(t<final_t) {
		//rate vector construction (or update after the first realization):
		//1-Cycle for the degradation rates of the 2n species; I distinguish between wildtype and mutants 
		//degradation rates
		for (i=0; i<2*n; i++) {
			if (i%2==0) ratev[i]=mu*x[i];
			else ratev[i] = mum*x[i];
		}

		//2-Cycle for the replicaton rates of the 2n species
		for (i=2*n; i<4*n; i++) {
			ratev[i]=x[i-2*n]*rep_rate(x[i-2*n-i%2],x[i-2*n-i%2+1],delta);
			//printf("%f \n",rep_rate(x[i-2*n-i%2],x[i-2*n-i%2+1],delta));
		
		}
		if (n>1) {
			//3-Cycle for the exchange of first cell of the chain toward the second
			for (i=4*n; i<4*n+2; i++) {
				ratev[i]=x[i-4*n]*k;
			}

			//4-Cycle for the exchange of the central cells (including of the second cell toward the first and of the  second last toward the last)
			for (i=4*n+2; i<8*n-6; i++) {
				ratev[i]=k*x[2+2*((i-4*n-2)/4)+(i-4*n-2)%2];
			}

			//5-Cycle for the exchange last one toward the second last
			for (i=8*n-6; i<8*n-4; i++) {
				ratev[i]=k*x[2*n-(8*n-5-i)-1];
		}
		}

		//rate vector has been built/updated

		//print the rate vector
		// for (i=0;i<8*n-4;i++){
		// 	printf("%f    ",ratev[i]);
		// }
		// printf("\n");
		
		
		rate =sum_array(ratev, 8*n-4);
		
      //-------------------------------------------------------------Gillespie's algorithm-----------------------------------------------------------
	    //extract a random number exponentially distributed with rate given by "rate" :
	    xi=1-RND;
		wait_time=-log(xi)/(float)rate; //inverse CDF method to get exponentially distributed number
		//this is the time after which the next event happens


		//now figure out which event has happened:
		
		xi=RND*rate;//extract uniform random number between 0 and rate;
		//printf("%f\n",xi);
		//printf("%f  %f\n",rate, xi);
		for ( i = 0; i < 8*n-4; i = i + 1 ){
			if (partsum(ratev,i) >= xi) {
				for (h=0; h<2*n; h++) {  
					x[h]=x[h]+V[h][i]; //update the vector X according to the event that occurred
					//printf("%f ", x[h]);
				}
				
				break;
			}
		}
		
		t=t+wait_time; //update time
		
       //-----------------------end of Gillespie's algorithm -------------------------------------
		

		//calculate heteroplasmies and populations for the various times
		if (t>times[j]) { //j is the time index
            for(i=0; i<2*n; i++) {
				populations[j][i] = populations[j][i]+x[i];
				popsquared[j][i]  = popsquared[j][i]+x[i]*x[i];

			}
				
			for (i=0; i<n; i++) {
				het[j][i]=het[j][i]+x[2*i+1]/(x[2*i]+x[2*i+1]);
				hetsquared[j][i] = hetsquared[j][i] +(x[2*i+1]/(x[2*i]+x[2*i+1]))*(x[2*i+1]/(x[2*i]+x[2*i+1]));

			}
           
			j++;
		}
		
					
		
		
	} //****************************************END OF LOOP OVER TIME*******************************************
	


	if (rl%nprint==0) {
		prg = fopen ("progress.txt", "w"); //just a file to check the progress of the simulation
		fp = fopen ("heteroplasmies.dat", "w+");
		pop = fopen ("populations.dat", "w+");

		// write on file to check progress
		fprintf(prg,"%i realizations out of %i \n",rl, realizations);
		fflush(prg);
		fclose(prg);
		printf("%i\n", rl);

		//Heading of file with heteroplasmies
		fprintf(fp,"#Heteroplasmies over time and their standard deviations: t, hmean, h1, \xCF\x83h1,  h2, \xCF\x83h2  and so on. \n \n");
		fprintf(fp,"#Initial condition for heteroplasmies (h1, h2 and so on) : ");
		for(z=0;z<n;z++) {
			fprintf(fp, "%5.3f  ",het[0][z]/(float)(rl));
		}
		fprintf(fp," \n");
		fprintf(fp,"#Initial condition for population (w1, m1, w2, m2  and so on) : ");
		for(z=0;z<2*n;z++) {
			fprintf(fp, "%5.1f  ",populations[0][z]/(float)(rl));
		}
		fprintf(fp," \n");
		fprintf(fp,"#Hopping rate \xC6\x94=%10.9f; selection against mutants of %8.7f\xC2\xB5; \xCE\xB4 =%5.3f \n", k,epsilon, delta );
		fprintf(fp,"#Final time=%6.0f,  Current number of realizations= %i; Total number of realizations= %i;\n", final_t, rl,realizations );


		//Heading of file with populations
		fprintf(pop,"#Populations of cells with stand dev, i.e. w1, \xCF\x83w1,  m1, \xCF\x83m1 etc. \n \n");
		fprintf(pop,"#Initial condition for heteroplasmies (h1, h2 and so on) : ");
		for(z=0;z<n;z++) {
			fprintf(pop, "%5.3f  ",het[0][z]/(float)(rl));
		}
		fprintf(pop," \n");
		fprintf(pop,"#Initial condition for population (w1, m1, w2, m2  and so on) : ");
		for(z=0;z<2*n;z++) {
			fprintf(pop, "%5.1f  ",populations[0][z]/(float)(rl));
		}
		fprintf(pop," \n");
		fprintf(pop,"#Hopping rate \xC6\x94=%10.9f; selection against mutants of %8.7f\xC2\xB5; \xCE\xB4=%5.3f \n", k,epsilon, delta );
		fprintf(pop,"#Final time=%6.0f,  Current number of realizations= %i; Total number of realizations= %i; \n", final_t, rl,realizations );

		

		//Now Write on files

		//here calculating the mean heteroplasmy of the n cells ...
		for (j=0;j<M+1;j++) { 
			temp_mean=0;
			for (z=0; z<n; z++) {
				temp_mean = temp_mean + het[j][z];
			}
			temp_mean = temp_mean/(float)n;	
			fprintf(fp, "%10.4f %12.6f ", times[j],temp_mean/(float)(rl));//... and writing it in the second column
			fprintf(pop, "%10.4f ", times[j]); 
			
			
			//calculating  std of the i-th cell heteroplasmy 
			//and writing it in the file after the respective heteroplasmies
			for (i=0;i<n;i++) {	
				hvars[j][i]=hetsquared[j][i]/(float)(rl)-(het[j][i]/(float)(rl))*het[j][i]/(float)(rl);
				hsigmas[j][i]=sqrt(hvars[j][i]);
				fprintf(fp, "%12.6f ", het[j][i]/(float)(rl)); //writing heteroplasmy
				fprintf(fp, "%12.6f ", hsigmas[j][i]/sqrt((float)(rl))); //writing corresponding std
				hvars[j][i]=0; // must be reset to zero
				hsigmas[j][i]=0; // must be reset to zero

			}


			//calculating  std of the i-th cell population
			//and writing it in the file after the respective population size
			for (i=0;i<2*n;i++) {

				popvars[j][i]=popsquared[j][i]/(float)(rl)-(populations[j][i]/(float)(rl))*populations[j][i]/(float)(rl);
				popsigmas[j][i]=sqrt(popvars[j][i]);	
				fprintf(pop, "%8.2f ", populations[j][i]/(float)(rl)); //writing population size (w and m for each cell)
				fprintf(pop, "%8.2f ", popsigmas[j][i]/sqrt((float)(rl))); //writing corresponding heteroplasmy
				popvars[j][i]=0;    // must be reset to zero
				popsigmas[j][i]=0;  // must be reset to zero
			}
			
			fprintf(fp,"\n");
			fprintf(pop,"\n");
		}

 		fclose(fp);
		fclose(pop);
	}

} 
//**********************************************************END OF LOOP OVER REALIZATIONS*******************************************************

return 0;	
} //END OF MAIN PROGRAM


//This function calculates the replication rate of the mtDNA molecules (specific to each cell) 
//implementing the linear feedback control system.
//Input : the number of wildtype (wt) and mutant (mt) of the interested cell, the parameter delta (ratio between control strenghts)
//Output : the replication rate for the interested cell with the current population
double rep_rate(double wt,double mt, double d){
	
	double Nss, c1;

	Nss=455;
	c1=0.001;

	if (0.07+c1*(Nss-(wt+d*mt)>0)) {
    	return 0.07+c1*(Nss-(wt+d*mt));	
	}
	else{
		return 0;
	}

	
}

//This function takes a vector a its length num_elements as inputs and gives back the sum of the 
//elements of the vector
double sum_array(double a[], int num_elements){	
   int i;
   double sum=0;
   for (i=0; i<num_elements; i++){
	 sum = sum + a[i];
   }
   return(sum);
}




//This function takes a vector and an index as inputs and gives the sum of the elements up to 
//that index.
//Inputs : vector tosum, index stop_element
//Ouput : sum of the elements up to the index stop_element
double partsum(double tosum[],int stop_element){
	double psum=0;
	int h=0;
	for (h=0; h<=stop_element; h=h+1){
			psum=psum+tosum[h];
	}
	return (psum);
}

