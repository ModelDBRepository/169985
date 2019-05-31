/*
 *  functions.cpp
 *  stnetwork
 *
 *  Created by Honi Sanders on 5/25/11.
 *  Copyright 2011 Brandeis University. All rights reserved.
 *
 *
 requires that AMPAscale, NMDAscale, GABAscale, etc. be declared somewhere else
 requires that EIF, HH be declared somewhere else, 
	one as 0, the other as 1 to determine which kind of simulation is being run
 requires that funcoutfiles and funcout be declared somewhere else
 
	all valarray matrices are listed column by column
	time moves forward column by column
	in connectivity valarrays connectivity is from rows to columns
	
	Output files have all data for a given time step on a single line. 
	Time progresses row by row.  The first entry on each row should be the time
	associated with that time step (implemented in main).
 
 
 > only changes necessary for the new (spatial-temporal) implementation are
 > 1. alpha_s = 0.5 in HHstep for 50% activation of NMDARs with single spike
 > 2. change constitutive from 0.25 to 0.05


 */

#include "network.h"
double E_revValues[] = {E_L, E_g_E, E_g_I, 0.0, E_K, E_Na, E_K};
extern double eeAMPAscale, eiAMPAscale, inscale, iniscale, eeNMDAscale, eiNMDAscale;
extern double GABAAscale, GABABscale;
extern double ffAMPAscale, ffNMDAscale;
extern bool EIF, HH;
extern int funcout;
extern ofstream funcoutfile1,funcoutfile2;  
double constitutive=0.05;	// fraction of KIR that is constitutively active





/********** Neuron Functions **********/


Neuron::Neuron() {
	this->initialize("not", 'n', 'n');
}


Neuron::Neuron(const char ty[], const char inputset, const char inhibitionset) {
	this->initialize(ty, inputset, inhibitionset);
}



void Neuron::initialize(const char ty[], const char inputset, const char inhibitionset) {
	strncpy(celltype, ty, 4);
	V_s = V_init + doubrand(-10, 10);
	V_d = V_init + doubrand(-10, 10);
	g_refr = 0;
	AMPA = s_init;
	NMDA = s_init;
	GABAA = s_init;
	GABAB = s_init;
	T = 0;  R = 0;  D = 0;  G = 0;	B = 0;
	inact = 0;
	spike = 0;
	for (int i=0; i<4; i++) {
		gating[i] = 0;
	}
	input = inputset;
	inhibition = inhibitionset;
}



void Neuron::fillwstats(double arr[7]) {
	arr[0] = V_s; arr[1] = V_d; arr[2] = spike; 
	arr[3] = AMPA; arr[4] = NMDA; arr[5] = GABAA; arr[6] = GABAB;
}



void Neuron::step(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB, 
						double g_in, double I_ext, double g_PING) {
	if (EIF) 
		EIFstep(g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in, I_ext, g_PING);
	else if (HH) 
		HHstep(g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in, I_ext, g_PING);
}






void Neuron::EIFstep(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB, 
							double g_in, double I_ext, double g_PING) {
	valarray<double> E_revs(E_revValues,5), g(E_revs), I(g); g=0; I=0;
	double g_E, g_I, V_d0=V_d, V_s0=V_s, refrfact;
	
	
	AMPA *= exp(-dt/tau_AMPA);		// relax synaptic activations, refractory conductance
	NMDA *= exp(-dt/tau_NMDA);
	GABAA *= exp(-dt/tau_GABAA);
	GABAB *= exp(-dt/tau_GABAB);
	if (strcmp(celltype, "exc") == 0) {refrfact = exp(-dt/tau_refrE);}
	else if (strcmp(celltype, "inh") == 0) {refrfact = exp(-dt/tau_refrI);}
	else {cerr << celltype << " "; refrfact=1;}
	g_refr *= refrfact;
	
	
	g_E = g_AMPA + vdepend(g_NMDA,'N') + w_noise*doubrand(-1, 1)/sqrt(dt);
	g_I = g_GABAA + vdepend(g_GABAB,'G') + w_noise*doubrand(-1, 1)/sqrt(dt);
	
	
	
	g[0]=g_L; g[1]=g_E; g[2]=g_I; g[3]=g_c; g[4] = g_refr;	//leak, exc, inh, refr, coupling conductances
	
	
	V_s = compartmentstep('s', V_s0, g, g_in, V_d0, I_ext, g_PING);	// new somatic voltage
	if (V_s > peak){						// spike?
		Neuron::spiking();
	}
	else if (V_s<thresh && V_s0>thresh) {	// Na+ deinactivation?
		inact = false;
	}
	else {
		spike = 0;
	}
	
	
	if (strcmp(celltype, "exc") == 0) {
		V_d = compartmentstep('d', V_d0, g, g_in, V_s0, 0, g_PING);	// new dendritic voltage
	}
	else {V_d = V_s;}										// no dendritic compartment for inh cells
	
	
	
	V_s = abs(V_s)/V_s * min(100.0,abs(V_s));	//restrict V_s to [-100,100]
	V_d = abs(V_d)/V_d * min(100.0,abs(V_d));	//restrict V_d to [-100,100]
}




void Neuron::HHstep(double g_AMPA, double g_NMDA, double g_GABAA, double g_GABAB, 
						  double g_in, double I_ext, double g_PING) {
	double alpha_m, beta_m, alpha_h, beta_h, alpha_n, beta_n;
	double  m=gating[0], h=gating[1], n=gating[2], x=gating[3];
	valarray<double> E_revs(E_revValues,7), g(E_revs); g=0;
	double g_E, g_I, V_d0=V_d, V_s0=V_s;
	
	
	g_E = g_AMPA + vdepend(g_NMDA,'N');// + w_noise*doubrand(-1, 1)/sqrt(dt);
	g_I = g_GABAA;// + w_noise*doubrand(-1, 1)/sqrt(dt);

	
	g[0]=g_L; g[1]=g_E; g[2]=g_I; g[3]=g_c;			//leak, exc, inh, coupling conductances
	if (strcmp(celltype, "exc") == 0) {			// include GABAB in excitatory neurons only
		g[6] = vdepend(g_GABAB,'G');	
	}
	
	if (funcout) {			// writes to file
		funcoutfile1 << g_in*(E_g_E-V_d) << " " << vdepend(g_NMDA,'N')*(E_g_E-V_d) << " ";	// currents
		funcoutfile1 << g_GABAA*(E_g_I-V_d) << " " << vdepend(g_GABAB,'G')*(E_K-V_d) << " ";
		if (funcout == 2) 
			funcoutfile2 << g_AMPA << " " << g_NMDA << " " << g_GABAA << " " << g_GABAB << " ";// conductances
	}	

	
	// Hodgkin-Huxley equations
	if (strcmp(celltype, "exc") == 0) {
		alpha_m = -0.1 * (V_s + 32)/(exp(-0.1*(V_s + 32)) - 1);
		beta_m = 4 * exp(-(V_s + 57)/18);
		alpha_h = 0.07 * exp(-(V_s + 48)/20);
		beta_h = 1 / (exp(-0.1*(V_s + 18)) + 1);
		m += (Phi_m * (alpha_m*(1 - m) - beta_m*m)) * dt;
		h += (Phi_h * (alpha_h*(1 - h) - beta_h*h)) * dt;
		g[5] = g_Na*pow(m,3)*h;						// Na conductance
		
		alpha_n = -0.01 * (V_s + 34)/(exp(-0.1*(V_s + 34))-1);
		beta_n = 0.125 * exp(-(V_s + 44)/80);
		n += (Phi_n * (alpha_n*(1 - n) - beta_n*n)) * dt;
		g[4] = g_K*pow(n,4);						// K conductance
	}
	else if (strcmp(celltype, "inh") == 0) {		// see wang & buzsaki, 1996
		alpha_m = -0.1 * (V_s + 35)/(exp(-0.1*(V_s + 35)) - 1);
		beta_m = 4 * exp(-(V_s + 60)/18);
		alpha_h = 0.07 * exp(-(V_s + 58)/20);
		beta_h = 1 / (exp(-0.1*(V_s + 28)) + 1);
		m = alpha_m/(alpha_m + beta_m);
		h += (Phi_n * (alpha_h*(1 - h) - beta_h*h)) * dt;
		g[5] = g_Na*pow(m,3)*h;				// Na conductance
		
		alpha_n = -0.01 * (V_s + 34)/(exp(-0.1*(V_s + 34))-1);
		beta_n = 0.125 * exp(-(V_s + 44)/80);
		n += (Phi_n * (alpha_n*(1 - n) - beta_n*n)) * dt;
		g[4] = g_K*pow(n,4);						// K conductance  (changed from 90 to 180)
	}
	else {
		cerr << celltype << " ";
	}

	
	V_s = compartmentstep('s', V_s0, g, g_in, V_d0, I_ext, g_PING);	// new somatic voltage
	if (strcmp(celltype, "exc") == 0) {
		V_d = compartmentstep('d', V_d0, g, g_in, V_s0, 0, g_PING);	// new dendritic voltage
	}
	else {V_d = V_s;}										// no dendritic compartment for inh cells
	
	V_s = abs(V_s)/V_s * min(100.0,abs(V_s));	//restrict V_s to [-100,100]
	V_d = abs(V_d)/V_d * min(100.0,abs(V_d));	//restrict V_d to [-100,100]

	
	
	
	if (spike == 0 && V_s > 0 && V_s > V_s0) {spike = -1;}		// note the action potential upswing
	else if (spike == -1 && V_s > 20 && V_s < V_s0) {spike = 1;}	// declare spike if V decreases after AP upswing
	else if (spike == 1) {spike = 2;}						// don't count spike again until repolarization
	else if (spike == 2 && V_s < 0) {spike = 0;}			// allow spike after repolarization


	
	// 1/(1 + exp(-V_s/2.0)*dt integrated over one spike ~= 0.3 for i-cells and 0.4 for e-cells												

	if (strcmp(celltype, "exc") == 0) {				// adjust excitatory synaptic gating variables
		AMPA += (kf * 1/(1 + exp(-V_s/2.0)) * (1 - AMPA) - kr * AMPA) * dt;
		
		double alpha_x = 10, beta_x = 0.5, alpha_s = 0.5, beta_s = 0.01;	// in kHz
		// alpha_s was originally 1, adjust to change amount of NMDA activated w/ one spike
		x += (alpha_x*1/(1 + exp(-V_s/2.0)) * (1 - x) - beta_x*x) * dt;
		NMDA += (alpha_s*x*(1 - NMDA) - beta_s*NMDA) * dt;
	}
	else if (strcmp(celltype, "inh") == 0) {		// adjust inhibitory synaptic gating
		GABAA += (kf * 1/(1 + exp(-V_s/2.0)) * (1 - GABAA) - beta * GABAA) * dt;
		
		bool Destexhe1995=0, Destexhe1999=!Destexhe1995;
		
		if (Destexhe1995) {
			double dT, dR, dD, dG, Kd=100, K1=6.6*100, K2=0.02, K3=0.0053, K4=0.017, K5=0.083, K6=0.0079;
			double Vmax=0.0001, Km=0.000004;
			// Comparison to orignal values in Destexhe & Sejnowski, 1995:  Kd is still in uM b/c G is in uM.
			// K1 is in 1/M*1/ms b/c dt is in ms and T is in M, so it is muliplied by 10^-3
			// K2, K3, K4, and K6 are in 1/ms b/c dt is in ms, so they are multiplied by 10^-3
			// K5 is in uM/ms b/c dt is in ms and G is in uM, so it is multiplied by 10^6*10^-3
			// Vmax is in M*1/ms b/c dt is in ms and T is in M, so it is multiplied by 10^-3
			// Km is in M b/c T is in M, so it is multiplied by 10^-6
			dT = - Vmax * T/(T+Km);  dR = K1*T*(1-R-D) - K2*R + K3*D;  dD = K4*R - K3*D;  dG = K5*R - K6*G;
			T += dT*dt;  R += dR*dt;  D += dD*dt;  G += dG*dt; // [GABA], activated, desensitized GABABR, G protein for GABAB current
			if (spike==1)  T += 0.007;    // GABA release w/ spike
			// doesn't much matter how much T changes w/spike unless there is large diffusion and few spikes
			GABAB = pow(G,4)/(pow(G,4) + Kd);
			// cout << T << " " << R << " " << D << " " << G << " ";
		}
		else if (Destexhe1999) {
			
			double dT, dB, dR, dG, kD=0.1, k1=30, kminus1=0.1, k2=0.02, Bm=1, K1=0.18, K2=0.0096, K3=0.19, K4=0.06, Kd=17.83;
			
			dT = -k1*T*(Bm-B) + kminus1*B - kD*T;	dB = k1*T*(Bm-B) - (kminus1+k2)*B;
			dR = K1*T*(1-R) - K2*R;	dG = K3*R - K4*G;
			T += dT*dt;  R += dR*dt;  B += dB*dt;  G += dG*dt;
			if (spike==1)  T += 1;    // GABA release w/ spike
			GABAB = pow(G,4)/(pow(G,4) + Kd);
			// cout << T << " " << B << " " << R << " " << G << " " << dT << " " << dB << " " << dR << " " << dG;
		}
		
      /*
       GABAB *= exp(-dt/tau_GABAB);
       if (spike==1) GABAB += alpha_GABAB * (1-GABAB);
       */
	}
	gating[0]=m; gating[1]=h; gating[2]=n, gating[3]=x;	// store gating variables

	
}






double Neuron::compartmentstep(char compartment, double V, valarray<double> g, 
										 double g_in, double V_c, double I_ext, double g_PING) {
	valarray<double> E_revs(E_revValues, g.size()), I(g);
	// g is [leak,exc,inh,coupling,K,Na,GABAB]
	double Psi=0, dV;
	bool onecomp = false;
	
	if (strcmp(celltype, "inh") == 0) {		// HH interneurons are messed up, so
		onecomp = true;			// note that there is only one compartment
		E_revs[4] = -90;		// change K reversal potential
	}
	
	E_revs[3] = V_c;			// potential of coupled compartment
	
	bool TTX=1;
	if (compartment == 'd' || (I_ext!=0 && TTX)) {	// in dendrite (or with I_ext in TTX), there is
		g[4] = 0;				// no refractory current if EIF and no HH currents if HH
		if (HH) {g[5] = 0;}
	}
	if (compartment == 's' && !onecomp) {	// in soma (except for HH inh cells), there is
		g[1] = 0;				// no excitatory synaptic current
		g[2] += g_PING;		// PING inhibitory conductance
		g[6] = 0;				// no GABA-B current
	}
	if (EIF && (onecomp || compartment == 's') ) {
		Psi = g_L * Delta * exp((V - thresh)/Delta);	// spike-generating current
	}
	if (compartment != inhibition && !onecomp) {
		g[2] = 0;
	}
	if (compartment == input || onecomp) {
		g[1] += g_in;
	}
	g[1] += w_noise*doubrand(-1, 1)/sqrt(dt);
	g[2] += w_noise*doubrand(-1, 1)/sqrt(dt);
	
	I = g*(E_revs - V);	// + w_noise*doubrand(-1,1)/sqrt(dt);			// current
	dV = (I.sum() + I_ext + Psi*!inact)*dt/Cm;
	V = V + dV;					//new compartment voltage	
	return V;
}



double Neuron::vdepend(double g, char type) {
	valarray<double> E_revs(E_revValues,7);
	double g_new, g_max=1;
	bool weak=0;		// use weakly rectifying GIRK?
	if (type == 'N') {
		g_max = 1/(1 + 0.3*Mg*exp(0.08*(0)));		// conductance at reversal potential
		g_new = g/(1 + 0.3*Mg*exp(0.08*(E_revs[1]-V_d))); // /g_max;
	}
	else if (type == 'G') {					// GABAB uses K+ reversal potential
		if (weak==1) {							// weak rectifier
			g_max = 1/(1 + exp(-0.05*(0)));		// conductance at reversal potential of weak rectifier
			g_new = g/(1 + exp(-0.07*(E_revs[6]-V_d))); // /g_max;		
		}
		else {								// strong rectifier
			g_max = 1/(1 + exp(-0.1*(-10)));	// conductance at reversal potential of strong rectifier
			g_new = g/(1 + exp(-0.1*(E_revs[6]-V_d-10))); // /g_max;	
		}
	}
	else {
		cerr << "Did not specify valid type for rectification." << endl;
	}
	return g_new;
}






void Neuron::spiking() {
	V_s = peak;
	if (strcmp(celltype,"exc") == 0) {
		AMPA += alpha * (1 - AMPA);
		NMDA += alpha * (1 - NMDA);
	}
	else if (strcmp(celltype, "inh") == 0) {
		GABAA += alpha * (1 - GABAA);
		GABAB += alpha_GABAB * (1 - GABAB);
	}
	
	g_refr = g_refr_max;
	spike = 1;
	inact = true;
}	















////////////  Network Functions //////////////



Network::Network() {
	N = 400;					// # of neurons in network
	N_E = (ratio*N/(ratio+1));	// # of excitatory neurons
	N_I = (N - N_E);			// # of inhibitory neurons
	N_in = 100;					// # of incoming axons
	allneurons = new Neuron[N];		// initialize an array of N neurons	
	
	// initialize synaptic weights
	w_AMPA.resize(N_E*N,0.0); w_NMDA.resize(N_E*N,0.0); 
	w_GABAA.resize(N_I*N,0.0); w_GABAB.resize(N_I*N,0.0); 
	w_in.resize(N_in*N,0.0); 
}






Network::Network(const size_t N_A, const valarray<double>& w_input, const double connectivity,
					  const char inputset, const char inhibitionset) {	// make new network of size N_A
	this->initialize(N_A, w_input, connectivity, inputset, inhibitionset);
}



Network::Network(const size_t N_A, const size_t N_input, const double connectivity,
					  const char inputset, const char inhibitionset) {	// make new network of size N_A
	this->initialize(N_A, N_input, connectivity, inputset, inhibitionset);
}



void Network::initialize(const size_t N_A, const size_t N_input, const double connectivity,
								 const char inputset, const char inhibitionset) {	// make new network of size N_A
	slice s;
	valarray<double> temp1, w_input;
	
	N = N_A;					// # of neurons in network
	N_E = (ratio*N/(ratio+1));	// # of excitatory neurons
	N_I = (N - N_E);			// # of inhibitory neurons
	N_in = N_input;				// # of incoming axons
	
	w_input.resize(N_in*N,0.0);
	s = slice(0, N_E-1, N_in+1);
	w_input[s] = inscale;								// 1 input cell to each E-cell
	s = slice(N_in*N_E, N_in*N_I, 1);
	temp1.resize(N_in*N_I, 0.0);
	w_input[s] = connect(temp1, 1.0)*iniscale;			// input cells to all I-cells
	
	this->initialize(N_A, w_input, connectivity, inputset, inhibitionset);
}



void Network::initialize(const size_t N_A, const valarray<double>& w_input, const double connectivity,
								 const char inputset, const char inhibitionset) {	// make new network of size N_A
	size_t i;
	slice s, EEslice, EIslice, IEslice, IIslice;
	valarray<double> temp1, temp2;
	
	
	N = N_A;					// # of neurons in network
	N_E = (ratio*N/(ratio+1));	// # of excitatory neurons
	N_I = (N - N_E);			// # of inhibitory neurons
	N_in = w_input.size()/N;				// # of incoming axons
	
	allneurons = new Neuron[N];		// initialize an array of N neurons
	for (i=0; i<N_E; i++) {
		allneurons[i].initialize("exc", inputset, inhibitionset);
	}
	for (i=N_E; i<N; i++) {
		allneurons[i].initialize("inh", inputset, inhibitionset);
	}	
	
	
	// initialize synaptic weights
	w_AMPA.resize(N_E*N,0.0); w_NMDA.resize(N_E*N,0.0); 
	w_GABAA.resize(N_I*N,0.0); w_GABAB.resize(N_I*N,0.0); 
	w_in.resize(N_in*N,0.0); 
	w_ffAMPA.resize(N_E*N,0.0); w_ffNMDA.resize(N_E*N,0.0);
	
	temp1.resize(N_E*N_E, 0.0);
	EEslice = slice(0,N_E*N_E,1); EIslice = slice(N_E*N_E,N_E*N_I,1); 
	IEslice = slice(0,N_I*N_E,1), IIslice = slice(N_I*N_E,N_I*N_I,1);
	w_AMPA[EEslice] = connect(temp1, connectivity) * eeAMPAscale;	// E-cells to E-cells
	temp2.resize(N_E*N_I,0.0);
	w_AMPA[EIslice] = connect(temp2, connectivity) * eiAMPAscale;	// E-cells to I-cells
	
	w_NMDA[EEslice] = temp1 * eeNMDAscale;							// E-cells to E-cells
	w_NMDA[EIslice] = temp2 * eiNMDAscale;							// E-cells to I-cells
	
	temp1.resize(N_E*N_I, 0.0);
	w_GABAA[IEslice] = connect(temp1, connectivity) * GABAAscale;					// I-cells to E-cells 
	temp2.resize(N_I*N_I, 0.0);
	//w_GABAA[IIslice] = connect(temp2, connectivity) GABAAscale;			// I-cells to I-cells
	w_GABAB[IEslice] = connect(temp1, connectivity) * GABABscale;			// I-cells to E-cells 
	
	
	w_in = w_input;
    
	w_ffAMPA = w_AMPA/eeAMPAscale*ffAMPAscale;							// last buffer to this buffer
	w_ffNMDA = w_NMDA/eeNMDAscale*ffNMDAscale;							// last buffer to this buffer	
	
	s = slice(0, N_E-1, N+1);
	w_AMPA[s] = 0;									// eliminate autapses
	w_NMDA[s] = 0;
}





void Network::fillwsize(size_t arr[]) {
	arr[0] = N; arr[1] = N_E; arr[2] = N_I; arr[3] = N_in;
}



	
double nullarray1[1]={0.0};
void Network::step(const valarray<double> syn_in, const double I_exts[]=nullarray1) {
	int i;
	double g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in=0, I_ext=0;
	double stats[7];
	valarray<double> syn_AMPA(N_E), syn_NMDA(N_E), syn_GABAA(N_I), syn_GABAB(N_I);
	valarray<double> syns_AMPA(N_E), syns_NMDA(N_E), syns_GABAA(N_I), syns_GABAB(N_I), syns_in(N_in);
	valarray<double> w_AMPAi(N_E), w_NMDAi(N_E), w_GABAAi(N_I), w_GABABi(N_I), w_ini(N_in);
	slice iEslice, iIslice, iinslice;

	
	for (i=0; i<N_E; i++) {			// gather local synaptic activations
		allneurons[i].fillwstats(stats);
		syn_AMPA[i] = stats[3];
		syn_NMDA[i] = stats[4];
	}
	for (i=0; i<N_I; i++) {
		allneurons[N_E+i].fillwstats(stats);
		syn_GABAA[i] = stats[5];
		syn_GABAB[i] = stats[6]*(1-constitutive) + constitutive;   // set constitutive=1 for fixed KIR
	}
	

	
	for (i=0; i<N; i++) {					// step over all neurons	
		iEslice =  slice(i*N_E, N_E, 1);
		iIslice =  slice(i*N_I, N_I, 1);
		iinslice = slice(i*N_in, N_in, 1);
		w_AMPAi = w_AMPA[iEslice];			// initialize weights
		w_NMDAi = w_NMDA[iEslice];
		w_ini = w_in[iinslice];
		w_GABAAi = w_GABAA[iIslice];
		w_GABABi = w_GABAB[iIslice];
		
		syns_AMPA = syn_AMPA * w_AMPAi;		// calculate synaptic conductances
		syns_NMDA = syn_NMDA * w_NMDAi;
		syns_GABAA = syn_GABAA * w_GABAAi;
		syns_GABAB = syn_GABAB * w_GABABi;
		
		syns_in = syn_in * w_ini;
		g_AMPA = syns_AMPA.sum();
		g_NMDA = syns_NMDA.sum();
		g_GABAA = syns_GABAA.sum();
		g_GABAB = syns_GABAB.sum();
		g_in = syns_in.sum();
		if (I_exts[0]!=0) {I_ext=I_exts[i];}

		
		allneurons[i].step(g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in, I_ext);	// calculate result of time step for that neuron

	}
	
	
}
	



void Network::step(const valarray<double>& syn_in, const valarray<double>& syn_ff, // N_E*2 long
						 const double I_exts[]=nullarray1, const double g_PING) {
	size_t i;
	double g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in=0, I_ext=0;
	double stats[7];
	valarray<double> syn_AMPA(N_E), syn_NMDA(N_E), syn_GABAA(N_I), syn_GABAB(N_I);
	valarray<double> syns_AMPA(N_E), syns_NMDA(N_E), syns_GABAA(N_I), syns_GABAB(N_I);
	valarray<double> syns_in(N_in), syns_ffAMPA(N_E), syns_ffNMDA(N_E);
	valarray<double> w_AMPAi(N_E), w_NMDAi(N_E), w_GABAAi(N_I), w_GABABi(N_I);
	valarray<double> w_ini(N_in), w_ffAMPAi(N_E), w_ffNMDAi(N_E);
	slice iEslice, iIslice, iinslice, ffAMPAslice(0,N_E,1), ffNMDAslice(N_E,N_E,1);
	
	
	for (i=0; i<N_E; i++) {			// gather local synaptic activations
		allneurons[i].fillwstats(stats);
		syn_AMPA[i] = stats[3];
		syn_NMDA[i] = stats[4];
	}
	for (i=0; i<N_I; i++) {
		allneurons[N_E+i].fillwstats(stats);
		syn_GABAA[i] = stats[5];
		syn_GABAB[i] = stats[6]*(1-constitutive) + constitutive;   // set constitutive=1 for fixed KIR
	}
	
	
	
	for (i=0; i<N; i++) {					// step over all neurons	
		iEslice =  slice(i*N_E, N_E, 1);
		iIslice =  slice(i*N_I, N_I, 1);
		iinslice = slice(i*N_in, N_in, 1);
		w_AMPAi = w_AMPA[iEslice];			// initialize weights
		w_NMDAi = w_NMDA[iEslice];
		w_ini = w_in[iinslice];
		w_ffAMPAi = w_ffAMPA[iEslice];
		w_ffNMDAi = w_ffNMDA[iEslice];
		w_GABAAi = w_GABAA[iIslice];
		w_GABABi = w_GABAB[iIslice];
		
		syns_AMPA = syn_AMPA * w_AMPAi;		// calculate synaptic conductances
		syns_NMDA = syn_NMDA * w_NMDAi;
		syns_GABAA = syn_GABAA * w_GABAAi;
		syns_GABAB = syn_GABAB * w_GABABi;
		syns_in = syn_in * w_ini;
		syns_ffAMPA = syn_ff[ffAMPAslice] * w_ffAMPAi;
		syns_ffNMDA = syn_ff[ffNMDAslice] * w_ffNMDAi;
		
		
		g_AMPA = syns_AMPA.sum() + syns_ffAMPA.sum();
		g_NMDA = syns_NMDA.sum() + syns_ffNMDA.sum();
		g_GABAA = syns_GABAA.sum();
		g_GABAB = syns_GABAB.sum();
		g_in = syns_in.sum();
		
		if (I_exts[0]!=0) {I_ext=I_exts[i];}
		
		// calculate result of time step for that neuron
		allneurons[i].step(g_AMPA, g_NMDA, g_GABAA, g_GABAB, g_in, I_ext, g_PING);	
		
	}
	
}











map<string,bool> arrayfileappend;

void outputvalarray(valarray<int> arr, size_t Nrows, size_t Ncols, string filename) {
	valarray<double> arr_doub(arr.size());
	for (size_t i=0; i<arr.size(); i++) {
		arr_doub[i] = double(arr[i]);
	}
	outputvalarray(arr_doub, Nrows, Ncols, filename);
}

void outputvalarray(valarray<bool> arr, size_t Nrows, size_t Ncols, string filename) {
	valarray<double> arr_doub(arr.size());
	for (size_t i=0; i<arr.size(); i++) {
		arr_doub[i] = double(arr[i]);
	}
	outputvalarray(arr_doub, Nrows, Ncols, filename);
}

void outputvalarray(valarray<double> arr, size_t Nrows, size_t Ncols, string filename) {
#ifdef __APPLE__
	// only write valarray to disk on laptop; not on cluster
	if (arr.size()!=Nrows*Ncols) {	// check that array is expected size
		cerr << "Nrows*Ncols = " << Nrows*Ncols << ", but arr.size() = " << arr.size() << endl;
		assert(arr.size()==Nrows*Ncols);
	}
	ofstream arrayfile;
	size_t i,j;
	if (arrayfileappend.count(filename)==0)		// append if file already opened, open if not
		arrayfile.open(filename.c_str());
	else
		arrayfile.open(filename.c_str(),ios_base::app);
	
	for (i = 0; i<Nrows; i++) {
		for (j = 0; j<Ncols; j++) {
			arrayfile << arr[i + j*Nrows] << " ";
		}
		arrayfile << endl;
	}
	arrayfile.close();
	arrayfileappend[filename]=true;		// next time, append instead of overwriting
#endif
}





valarray<double> loadvalarray(string filename, const size_t Nrows, const size_t Ncols) {
   ifstream infile(filename.c_str());
   istringstream *lines = new istringstream[Nrows];
   string line;
   valarray<double> arr(Nrows*Ncols);
   size_t count=0;
	bool err1=0,err2=0,err3=0,err4=0;
   
   if (infile.is_open()) {
		
      while (getline(infile, line)) {
         if (count == Nrows && !err1) {
            cerr << "Valarray file has more rows than expected." << endl;
				err1 = 1;
			}
         lines[count].str(line);
         count++;
      }		// read "filename" into "lines" (arr of istringstreams)
		
      if (count < Nrows && !err2) {
         cerr << "Valarray file has fewer rows than expected." << endl;
			err2 = 1;
		}
		
      infile.close();
		
      for (size_t col=0; col<Ncols; col++) {
         for (size_t row=0; row<Nrows; row++) {
            if (lines[row].good())
               lines[row] >> arr[Nrows*col + row];
            else if (!err3) {
               cerr << "Valarray file has fewer columns than expected." << endl;
					err3 = 1;
				}
				
				
            if (!err4 && col == Ncols-1 && !lines[row].eof()) {
					double entry;
					lines[row] >> entry;
					if (abs(entry)>eps || !lines[row].eof()) {	// taking care of bug in Mac compiler
						cerr << "Valarray file has more columns than expected. ";
						cerr << "Row " <<row << " has leftover " << entry << endl;
						err4=1;
					}
				}
         }
      }		// transfer "lines" to "arr" (valarray)
	}
	else cerr << "Unable to open valarray input file." << endl;

	
	return arr;
}




valarray<double> valarraysort(valarray<double> &arr) {
	size_t N=arr.size();
	if (N<=1)
		return arr;
	
	vector<double> vec(N);
	for (size_t i=0; i<N; i++) {
		vec[i] = arr[i];
	}
	sort(vec.begin(),vec.end());
	for (size_t i=0; i<N; i++) {
		arr[i] = vec[i];
	}
	return arr;
}




	
valarray<double> connect(valarray<double>& w, 
								 double connection) 
								 { // returns a matrix with non-zero values at a frequency of "connection"
	for (size_t i=0; i<w.size(); i++) {
		w[i] = doubrand();
		w[i] = rectify(w[i] - 1 + connection)/connection;	
		// is in [0,1] at a rate of "connection", =0 at a rate of 1-"connection"
		if (w[i]>0) {w[i] = 1.0;}			// =1 at a rate of "connection", =0 at a rate of 1-"connection"
	}
	return w;
}
	


double doubrand(double minimum, double maximum) {		// random double between minimum and maximum (defaults [0,1])
	double a;
	a = minimum + double(rand())/RAND_MAX * (maximum - minimum);
	return a;
}


double normrand(double mean, double std) {		// normally distributed random double with default mean=0, std dev=1
	double U=doubrand(), V=doubrand(),X;
	X = sqrt(-2*log(U))*cos(2*pi*V);
	return mean + X*std;
}



double rectify(double x) {		// =x if x>0, =0 if x<0
	x = (x + abs(x))/2;
	return x;
}

	
