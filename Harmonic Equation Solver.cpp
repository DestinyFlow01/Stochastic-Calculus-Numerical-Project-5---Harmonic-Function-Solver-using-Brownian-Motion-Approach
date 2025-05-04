#include<iostream>
#include<random>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<omp.h>
#include<Python.h>
using namespace std;

random_device rd; // Seed with a real random value, if available
thread_local mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

class Particle {
	private : 
		//Position and increment for each time step
		vector<double> x; 
		vector<double> dx; 
		vector<double> y; 
		vector<double> dy;
		double t = 0; 	//Current time
		
		double h_val = 0; //Initialize to zero for default
	public : 
		void Initialization(double, double);
		void Brownian_Motion(double dt); // Generating Brownian motion from current time step to next time step
		double Print_x_singletimestep(int i);
		vector<double> Print_x();
		double Print_h();
};


//Defining member function for class Particle
void Particle::Initialization(double x_init, double y_init) { x.push_back(x_init); y.push_back(y_init);}

void Particle::Brownian_Motion(double dt) {
	//Generating Brownian Motion
	double mean = 0, std = sqrt(dt);
	
	//Generating normal distribution
	normal_distribution<> d(mean, std);
	
	while(true) {
		//Sampling only 1 random number for Brownian motion : 
		double dx_current = d(gen);
		dx.push_back(dx_current);
		
		double dy_current = d(gen);
		dy.push_back(dy_current);
		
		// Ensure the vector `x` is not empty
	    if (x.empty()) {
	        throw runtime_error("Vector x is empty. Ensure Initialization() is called before Brownian_Bridge().");
	    }
		
		// Ensure the vector `y` is not empty
	    if (y.empty()) {
	        throw runtime_error("Vector y is empty. Ensure Initialization() is called before Brownian_Bridge().");
	    }
		
		double x_current = x.back() + dx_current;
		x.push_back(x_current);
		
		double y_current = y.back() + dy_current;
		y.push_back(y_current);
		
		double position = x_current*x_current + y_current*y_current; 
		
		if(position >= 1) {
			if(y_current >= 0) h_val = 1;
			else h_val = -1;
			
			break;
		}
	}
	
	
}

double Particle::Print_x_singletimestep(int i) { return x[i]; }

vector<double> Particle::Print_x() { return x; }

double Particle::Print_h() { return h_val; }

//Analytical solution : 
vector<vector<double>> Analytical_Field(int Nx, int Ny, double xmin, double xmax, double ymin, double ymax) {
	vector<vector<double>> h_field(Ny, vector<double>(Nx, 0));
	
	for(int i = 0; i<Nx; i++) {
		double x = xmin + (xmax - xmin)*i/(Nx-1);
		for(int j = 0; j<Ny; j++) {
			double y = ymin + (ymax - ymin)*j/(Ny-1);
			double r = sqrt(x*x + y*y);
			
			if(r >= 1) {
				if(y >= 0) h_field[j][i] = 1;
				else h_field[j][i] = -1;
			}
			else {
				double theta = atan2(y,x);
				
				double sum = 0;
				for(int k = 0; k<100; k++) {
					int n = 2*k+1;
					sum += 4/(M_PI*n)*pow(r,n)*sin(n*theta);
				}
				
				h_field[j][i] = sum;
			}
		}
	}
	return h_field;
}

//Absolute error with respect to analytical solution
vector<vector<double>> Absolute_Error(vector<vector<double>> Analytic, vector<vector<double>> Numeric, int Nx, int Ny) {
	vector<vector<double>> error_field(Ny, vector<double>(Nx, 0));
	
	for(int i = 0; i<Nx; i++) {
		for(int j = 0; j<Ny; j++) error_field[j][i] = abs(Analytic[j][i] - Numeric[j][i]);
	}
	
	return error_field;
}

//Output CSV
void OutputCSV_2D_Field(string name, vector<vector<double>> h_field, int Nx, int Ny, double xmin, double xmax, double ymin, double ymax) {
	ofstream ofs;
	ofs.open(name);
	
	ofs << "x,y,h\n";
	
	for(int i = 0; i<Nx; i++) {
		double x = xmin + (xmax - xmin)*i/(Nx-1);
		for(int j = 0; j<Ny; j++) {	
			double y = ymin + (ymax - ymin)*j/(Ny-1);
			ofs<<x<<","<<y<<","<<h_field[j][i]<<"\n";
			
			if(x == 0 and y == 0.5) {
				cout<<"h field at position ("<<x<<", "<<y<<") = "<<h_field[j][i]<<endl;
			}
		}
	}
	
	ofs.close();
}


//Plotting using Python
void call_Python() {
	// Python code for plotting
    const char* python_code = R"(
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as patches

df1 = pd.read_csv("Harmonic function field.csv")
df2 = pd.read_csv("Harmonic function field (Analytic).csv")
df3 = pd.read_csv("Harmonic function field (Error).csv")

heatmap_data1 = df1.pivot(index="y", columns="x", values="h")
heatmap_data2 = df2.pivot(index="y", columns="x", values="h")
heatmap_data3 = df3.pivot(index="y", columns="x", values="h")

fig, ax = plt.subplots(figsize=(10, 8))

sns.heatmap(heatmap_data1, cmap="coolwarm", annot=False, ax=ax)

xmin, xmax = df1["x"].min(), df1["x"].max()
ymin, ymax = df1["y"].min(), df1["y"].max()

plt.xlabel("X Axis")
plt.ylabel("Y Axis")
plt.title('Harmonic Function Field (Numeric)')
plt.show(block=False)



fig, ax = plt.subplots(figsize=(10, 8))

sns.heatmap(heatmap_data2, cmap="coolwarm", annot=False, ax=ax)

xmin, xmax = df2["x"].min(), df2["x"].max()
ymin, ymax = df2["y"].min(), df2["y"].max()

plt.xlabel("X Axis")
plt.ylabel("Y Axis")
plt.title('Harmonic Function Field (Analytic)')
plt.show(block=False)


fig, ax = plt.subplots(figsize=(10, 8))

sns.heatmap(heatmap_data3, cmap="coolwarm", annot=False, ax=ax)

xmin, xmax = df3["x"].min(), df3["x"].max()
ymin, ymax = df3["y"].min(), df3["y"].max()

plt.xlabel("X Axis")
plt.ylabel("Y Axis")
plt.title('Harmonic Function Field (Error)')
plt.show(block=False)


input("Press Enter to close all plots...")
    )";

    // Execute the Python code
    PyRun_SimpleString(python_code);
}

int main() {	
	system("CLS");
	//Parameters
	double t_initial = 0; 
	double dt = 0.01;
	int t_out = 1000;
	int num_thread = 13;
	
	
	int num_particle_prompt;
	cout<<"Number of particles (each) = "; cin>>num_particle_prompt;
    const int num_particle = num_particle_prompt;
    vector<Particle> particles(num_particle);
    
    //Simulation domain
    double xmin = -1, xmax = 1;
    double ymin = -1, ymax = 1;
    double dx = 0.01, dy = 0.01; 
    int Nx = (xmax - xmin)/dx + 1, Ny = (ymax - ymin)/dy + 1;
    
    //Making field of h : 
    vector<vector<double>> h_field(Ny, vector<double>(Nx, 0));
    
    //Numerical simulation
    cout<<"Simulating for each x and y position : \n";
    
    //Set the number of threads
    omp_set_num_threads(num_thread);
	
	#pragma omp parallel for collapse(2)
    for(int i = 0; i<Nx; i++) {
    	for(int j = 0; j<Ny; j++) {
    		vector<Particle> particles(num_particle);
    		double x_init = xmin + (xmax - xmin)*i/(Nx-1);
			double y_init = ymin + (ymax - ymin)*j/(Ny-1);
			double sum_h = 0;
			
			//Simulating for each particles
    		for(int k = 0; k<num_particle; k++) {
    			particles[k].Initialization(x_init, y_init);
    			particles[k].Brownian_Motion(dt);
    			
    			double h_eachparticle = particles[k].Print_h();
    			sum_h += h_eachparticle;
			}
			
			h_field[j][i] = sum_h/num_particle;
		}
	}
	
	//Analytical solution : 
	vector<vector<double>> h_field_analytic = Analytical_Field(Nx, Ny, xmin, xmax, ymin, ymax);
	vector<vector<double>> h_field_error = Absolute_Error(h_field_analytic, h_field, Nx, Ny);
	
    //Output CSV
	cout<<"Outputing to CSV : \n";
	cout<<"- Values of h (Numeric) : ";
	string name = "Harmonic function field.csv"; OutputCSV_2D_Field(name, h_field, Nx, Ny, xmin, xmax, ymin, ymax);  
	
	cout<<"\n- Values of h (Analytic) : ";
	name = "Harmonic function field (Analytic).csv"; OutputCSV_2D_Field(name, h_field_analytic, Nx, Ny, xmin, xmax, ymin, ymax);  
	
	cout<<"\n- Values of h (Absolute error) : ";
	name = "Harmonic function field (Error).csv"; OutputCSV_2D_Field(name, h_field_error, Nx, Ny, xmin, xmax, ymin, ymax);  
	
	cout<<"Done\n\n";
	
	//Calling Python for pkotting
	cout<<"Plotting using Python\n";
    Py_Initialize();
    call_Python();
    Py_Finalize();
    
    /*
    Type this in CMD
    g++ -fopenmp "Harmonic Equation Solver.cpp" -IC:\Python313\include -LC:\Python313\libs -lpython313 -o "Harmonic Equation Solver" -Wl,--enable-auto-import
	"Harmonic Equation Solver.exe"
    */
}