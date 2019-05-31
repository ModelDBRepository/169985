/*
 *  network.h
 *  stnetwork
 *
 *  Created by Honi Sanders on 8/30/12.
 *  Copyright 2012 Brandeis University. All rights reserved.
 *
 */

#ifndef NEWNETWORK_H
#ifndef STNETWORK_H
#ifndef NETWORK_H
#define NETWORK_H

#include <cmath>
#include <valarray>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <ctime>
#include <map>
// // if not debugging:
// #define NDEBUG
#include <assert.h>
using namespace std;



class Neuron {
private:
	char celltype[4];	// exc or inh
	double V_s;			// Somatic voltage
	double V_d;			// Dendritic voltage
	int spike;			// did it spike? 1 or 0
	double AMPA;		// postsynaptic activation of AMPA receptors
	double NMDA;		// postsynaptic activation of NMDA receptors
	double GABAA;		// postsynaptic activation of GABAA receptors
	double GABAB;		// postsynaptic activation of GABAB receptors
	double T,R,D,G,B;	// [GABA], activated & desensitized GABABR, [G protein] for GIRK
	double g_refr;		// refractory conductance
	bool inact;			// inactivation of Na+ channels
	char input;				// are input synapses onto dendrite 'd' or soma 's'
	char inhibition;		// are inhibitory synapses onto dendrite 'd' or soma 's'
	double gating[4];	// m, h, n for HH equations and x for synaptic gating
	void spiking();
	double vdepend(double g, char type); // type is 'N' for NMDA or 'G' for GABAB
	double compartmentstep(char compartment, double V, valarray<double> g, 
								  double g_in, double V_coupling, double I_ext, double g_PING);
	void EIFstep(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB, 
					 double g_in, double I_ext, double g_PING);
	void HHstep(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB, 
					double g_in, double I_ext, double g_PING);
public:
	Neuron();			// default constructor
	Neuron(const char ty[], const char inputset, const char inhibitionset);		// constructor
	void initialize(const char ty[], const char inputset='s', const char inhibitionset='s');
	void fillwstats(double arr[7]);	// fills arr w/ V_s, V_d, spike, AMPA, NMDA, GABAA, GABAB
	void step(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB,
				 double g_in, double I_ext=0, double g_PING=0);
};



class Network {
private:
	size_t N;							// # of neurons in network
	size_t N_E;							// # of excitatory neurons
	size_t N_I;							// # of inhibitory neurons
	size_t N_in;						// # of incoming axons		
	valarray<double> w_AMPA;	// weights of local excitatory synapses
	valarray<double> w_NMDA;	// weights of local excitatory synapses
	valarray<double> w_GABAA;	// weights of local inhibitory synapses
	valarray<double> w_GABAB;	// weights of local inhibitory synapses
	valarray<double> w_in;		// weights of input synapses
	valarray<double> w_ffAMPA;	// weights of feedforward (buffer-buffer) synapses
	valarray<double> w_ffNMDA;	// weights of feedforward (buffer-buffer) synapses
public:
	Network();				// default constructor
	Network(const size_t N_A, const size_t N_input, const double connectivity,
			  const char inputset, const char inhibitionset);	// constructor
	Network(const size_t N_A, const valarray<double>& w_input, const double connectivity,
			  const char inputset, const char inhibitionset);	// constructor
	void initialize(const size_t N_A, const size_t N_input, const double connectivity,
						 const char inputset, const char inhibitionset);
	void initialize(const size_t N_A, const valarray<double>& w_input, const double connectivity,
						 const char inputset, const char inhibitionset);
	void fillwsize(size_t arr[]);
	Neuron* allneurons;		// array with all of the neurons in the network
	void step(const valarray<double> syn_in, const double I_exts[]);
	void step(const valarray<double>& syn_in, const valarray<double>& syn_ff, 
				 const double I_exts[], const double g_PING=0);
};

// time in ms, voltage in mV

#define Cm		10 			//Membrane Capacitance 		in nF/mm2
#define g_L		1			//Leak Conductance			in uS/mm2		
#define g_refr_max	50	//Refractory Conductance	     "
#define g_c 1 //changed from 10			//Soma-Dendrite Coupling Conductance "
//#define p	0.5			//directionality of Soma-Dendrite Coupling
#define	g_Na	450  // changed from 450		// Sodium Conductance for HH
#define g_K		180	 // changed from 180	// Potassium Conductance for HH
#define Delta	3.5			//Constant for EIF

#define kf	12			// synaptic opening variable in kHz (AMPA and GABAA)
#define kr	1			// synaptic closing variable in kHz (AMPA)
#define beta 0.1		// synaptic closing variable in kHz (GABAA)

#define E_L	-80				// leak reversal potential in mV
#define E_g_E	0			// excitatory  "			"
#define E_g_I	-70			// inhibitory  "			"
#define E_Na	55			// Sodium		"			"	(for HH)
#define E_K		-90			// Potassium	"			"		"  and refractory reversal for EIF
#define Phi_m	2.5			// Temperature Factors				"
#define Phi_h	2.5			
#define Phi_n	5			


#define tau_AMPA	2		// time constant of AMPA channel closing in ms
#define tau_NMDA	100		//	"	of NMDA				"			"
#define tau_GABAA	10		//	"	of GABAA			"			"
#define tau_GABAB	250		//	"	of GABAB			"			"
#define tau_refrI	1		//	"	of refractory current in I-cells
#define tau_refrE	2		//	"		"	      in E-cells

#define thresh	-55			// spiking threshold
#define reset	-70			// reset potential after spike
#define peak	50			// spike peak potential
#define alpha	0.5			// ratio of unactivated receptors activated by spike
#define alpha_GABAB	0.03	//	"			"	GABAB	"			"		"

#define Mg	0.5			//[Mg2+] in mM



#define V_init -70		// initial voltage
#define s_init	0		// initial synaptic activation



#define dt	0.025			//length of time step in ms



#define w_noise	0.5		// variance of noise over a millisecond (uS/mm2)
#define scale	g_L



#define ratio	4				//ratio:1::excitatory:inhibitory neurons



#define eps	0.000000000000001
#define pi	3.141592653589793



double doubrand(double minimum=0, double maximum=1); 	// random double between minimum and maximum
double normrand(double mean=0, double std=1);	// normally distributed random double
double rectify(double x);								// =x if x>0, =0 if x<0
valarray<double> connect(valarray<double>& w, double connection);	
// changes valarray to a connection matrix with connection probability of "connection"
void outputvalarray(valarray<double> arr, size_t Nrows, size_t Ncols,
						  string filename = "./results/array.dsv");
void outputvalarray(valarray<int> arr, size_t Nrows, size_t Ncols,
						  string filename = "./results/array.dsv");
void outputvalarray(valarray<bool> arr, size_t Nrows, size_t Ncols,
						  string filename = "./results/array.dsv");
valarray<double> loadvalarray(string filename, size_t Nrows, size_t Ncols);
valarray<double> valarraysort(valarray<double> &arr);




#endif
#endif
#endif
