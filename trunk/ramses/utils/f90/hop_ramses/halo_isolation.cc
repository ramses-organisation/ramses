#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>

using namespace std;

//const int N=73743;  


int main(int argc, char *argv[]){

//If to few arguments are specified in the command line

if(argc!=6){
	cout << endl;
	cout << "To few arguments" << endl;
	return 0;
}

//Bash line arguments

string filename=argv[1];
string lowermass=argv[2];
string uppermass=argv[3];
string distcond=argv[4];
string masscond=argv[5];

stringstream ein1;
stringstream ein2;
stringstream ein3;
stringstream ein4;
stringstream ein5;

//Count lines (halos) in the .pos file

char dateiname[20];

ein1 << filename;
ein1 >> dateiname;


string t;
int lines=0;
ifstream inFile(dateiname);
while(!inFile.eof()){
	getline(inFile,t);
	lines++;
}
cout << endl;
cout<< "The file " << dateiname << " contains " << lines-2 <<" halos" << endl;


// Definition of entities 

int N=lines-2;

double d;
string h;

double * I= new double [N];
double * M= new double [N];
double * X= new double [N];
double * Y= new double [N];
double * Z= new double [N];
double * U= new double [N];
double * V= new double [N];
double * W= new double [N];
double * CONT= new double [N];
double * R200cubed= new double [N];
double * R200= new double [N];
int * Npart= new int [N];
bool * suit = new bool [N];

double low;
double up;
//double low=5.87849186354*pow(10,-5); //this is 5*10^12 Solarmasses/h in code unit
//double up=5.87849186354*pow(10,-4);  //this is 5*10^13 Solarmasses/h in code unit

int start=0;		//number of the first halo in the list in the desired mass range
int sumsuit=0;		//number of halos in the desired mass range
int sumsuitall=0;	//number of halos in the desired mass range which are fulfilling the distance and mass condition

int a,b,a_sq;

double D=0;		//Distance between the COM of two haloes
double Xnew=0;
double Ynew=0;
double Znew=0;

ifstream eingabe;
ofstream ausgabe;


//Read .pos file

eingabe.open(dateiname);

eingabe >> h >> h >> h >> h >> h >> h >> h >> h >> h >> h ;

for(int i=0;i<N;i++){
	eingabe >> I[i] >> Npart[i] >> M[i] >> CONT[i] >>  X[i] >> Y[i] >> Z[i] >> U[i] >> V[i] >> W[i];
}

eingabe.close(); 


//Calculate R200 from the halo mass

for(int i=0;i<N;i++){

	R200cubed[i]=0.75*(1/3.14159265359)*0.005*M[i]; 
	R200[i]=pow(R200cubed[i],(1.0/3.0));

}


//Set bool array

for(int i=0;i<N;i++){
	suit[i]=false;
}

//Find Halos with desired masses M_e[5*10^12 Solarmasses/h ; 5*10^13 Solarmasses/h]

ein2 << lowermass;
ein2 >> low;
ein3 << uppermass;
ein3 >> up;

for(int i=0;i<N;i++){
	if(M[i]>=low){
		if(M[i]<up){
			//cout << i << "    " << M[i] << "    " << setw(10) << setprecision(8) << R200[i] <<  endl;
			sumsuit++;
			suit[i]=true;
			if(sumsuit==1){
				start=i;
			}									
		}
	}
}

//cout << start << endl;

ein4 << distcond; 
ein4 >> a;
ein5 << masscond; 
ein5 >> b;

//a=13; //Distance criterion parameter
//b=3; //Mass criterion parameter

//for(int a=0;a<13;a++){
	a_sq=a*a;
	sumsuitall=sumsuit;

	//for(int i=start;i<(start+sumsuit);i++){
		//suit[i]=true;
	//}


//Loop to find the Halos, among the ones selected above, which are fulfilling the distance and mass isolation criteria

	for(int j=start;j<(start+sumsuit);j++){
	
		//cout << R200[j] << endl;
		for(int i=0;i<N;i++){
		
			Xnew=X[i];
			Ynew=Y[i];
			Znew=Z[i];
	
			if((X[j]-X[i])>0.5){
				Xnew=X[i]+1;
			}
			if((X[i]-X[j])>0.5){
				Xnew=X[i]-1;
			}
			if((Y[j]-Y[i])>0.5){
				Ynew=Y[i]+1;
			}
			if((Y[i]-Y[j])>0.5){
				Ynew=Y[i]-1;
			}
			if((Z[j]-Z[i])>0.5){
				Znew=Z[i]+1;
			}
			if((Z[i]-Z[j])>0.5){
				Znew=Z[i]-1;
			}

			D=((X[j]-Xnew)*(X[j]-Xnew))+((Y[j]-Ynew)*(Y[j]-Ynew))+((Z[j]-Znew)*(Z[j]-Znew));

			if(D<(a_sq*((R200[j]+R200[i])*(R200[j]+R200[i])))&&M[j]<(b*M[i])){
				if(j!=i){
					suit[j]=false;				
				}		
			}		
	
			//Check
			/*
			if(j==925&&i==42015){		
			cout << setw(4) << j << "     " << setw(6) << i << "        " << setw(10) << setprecision(7) << D << "        " << setw(10) << setprecision(7) << 4*((R200[j]+R200[i])*(R200[j]+R200[i])) << "        " << setw(10) << setprecision(7) << D-4*((R200[j]+R200[i])*(R200[j]+R200[i])) << "        " << suit[j] << endl;
			}
			*/
		}

		if(suit[j]==false){
		sumsuitall--;
		}

	}


	//File to check of the bool array
	/*
	ausgabe.open("halo_suit.dat");

	for(int i=start;i<(start+sumsuit);i++){
	
		
		if(suit[i]==true){
		ausgabe << 1 << endl;
		}
		else{
			ausgabe << 0 << endl;
		}
		
		ausgabe << setw(4) << i << setw(4) << suit[i] << endl;
	}
	
	ausgabe.close();
	*/
	cout << endl;
	cout << "      #   npart       mass      contamin.      R200             x          y          z           u           v           w" << endl;
	//ausgabe.open("halo_isolation.dat");

	for(int i=start;i<(start+sumsuit);i++){
		/*
		if(suit[i]==true){
		ausgabe << 1 << endl;
		}
		else{
			ausgabe << 0 << endl;
		}*/
		if(suit[i]==true){	
		cout << setw(7) << i+1 << setw(8) << Npart[i]  << setw(13) << setprecision(7) << M[i] << setw(13) << setprecision(7) << CONT[i] << setw(13) << setprecision(7) << R200[i] << setw(11) << setprecision(7) << X[i] << setw(11) << setprecision(7) << Y[i] << setw(11) << setprecision(7) << Z[i] << setw(12) << setprecision(7) << U[i] << setw(12) << setprecision(7) << V[i] << setw(12) << setprecision(7) << W[i] << endl;
		}
	}

//	ausgabe.close();


	//cout << a << endl;
	cout << endl;
	cout << "Total number of halos with masses in the specified mass range: " << sumsuit << endl;  
	cout << "Number of those halos which are fulfilling the distance ("<< a << ") and mass (" << b << ") condition: " << sumsuitall << endl;
	cout << endl; 

	//}

}
