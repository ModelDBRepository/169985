// maybe change sweep output to accord with newnetwork's sweep output
#include "network.h"
const size_t N_A = 400, N_E = size_t(N_A*ratio/(ratio+1)), N_I = N_A - N_E;
size_t N_in = N_E, N_prf = size_t(N_A*0.25);
const int Nnetworks=3, Nstims=3;
int gamma_prime = 0;			// gamma cycle when external priming (syn_ff0) starts, index starts at 0


bool EIF = 0, HH = !EIF;
int funcout = 2;  // what to funcoutfile? 
// 0--nothing, 1--currents, 2--currents & (non v-dependent) conductances

const bool PING = 0;		// is there PIN gamma?

const double t_max = 160;
const size_t Ntsteps = size_t(t_max/dt);
double t_rest = 15;		// time before first feedforward spike
double t_stims[Nstims], t_ends[Nstims];		// times of external input starts and stops
double t_length = 5;		// length of external input
double t_gamma = 30;			// length of gamma cycle in ms

/*
 const double t_stim = 20;				//time stimulus goes on
 const double t_end = t_stim + 10;				//time stimulus goes off
 const double t_stim2 = t_end + 20;		//time 2nd stimulus goes on
 const double t_end2 = t_stim2 + (t_end-t_stim);			//time 2nd stimulus goes off
 */		// t_stim & t_end for only 2 stims
const double f_in_off = eps;				//input firing rate
const double f_in_on = 500;				//input firing rate for active cells

const double connectivity = 1;
const double act = double(N_E*connectivity);
const char input = 'd', inhibition = 'd';

double ee = 1.0, ff = 1.0;  // ee/ff is the ratio of recurrent to feed-forward input
double eeNMDAscale = 90.0/act * ee/(ee+ff), eiNMDAscale = 3.0/act * ee/(ee+ff);
double eeAMPAscale = eeNMDAscale/4, eiAMPAscale = eeNMDAscale/4;
//double eeAMPAscale = 40.0/act * ee/(ee+ff), eiAMPAscale = 10.0/act * ee/(ee+ff);
double GABAAscale = 2.0/act*ratio, GABABscale = 1300.0/act*ratio;
double inscale = 4.0, iniscale = 2.0/N_E;//inscale/N_E/6;
double ffAMPAscale = eeAMPAscale * ff/ee, ffNMDAscale = eeNMDAscale * ff/ee;
double PING_AMPAscale = eiAMPAscale*4, PING_NMDAscale = eiNMDAscale*4, PING_GABAscale = 20.0;
/*
 // from Fellous: 
 double eeAMPAscale=0.4/act, eiAMPAscale=0.1/act, eeNMDAscale = 30.0/act, eiNMDAscale=3.0/act;			
 double GABAAscale = 1.0/act, GABABscale=0;
 double inscale = 32.0;
 */		// from Fellous


//output files
ofstream AV_dfile, AV_sfile, inspikesfile, Aspikecoordsfile;
ofstream Atransmitfile, sweepfile, outfile, funcoutfile1, funcoutfile2;



valarray<double> initinput(gslice_array<double> ginput, 
													 double f_base, double t_on, double f_on, double t_off);
valarray<double> spikes2act(valarray<double> spikes, double alpha_x = 0.5, double tau = tau_AMPA, bool w=0);

void openfiles();
void record(double stats[7], size_t i, size_t j, double t);
void closefiles();












int main (int argc, char * const argv[]) {
	/******* READ COMMAND LINE INPUT *********/
	if (argc > 1 && *argv[1] != '_') {
		eeNMDAscale = atof(argv[1])/act * ee/(ee+ff);	
		//eiNMDAscale = eeNMDAscale/10;
		eeAMPAscale = eeNMDAscale/4;	// comment if you don't want AMPA to scale with NMDA
		//		inscale = max(inscale, eeAMPAscale*act);
		eiAMPAscale = eeNMDAscale/4;
		ffAMPAscale = eeAMPAscale * ff/ee;
		ffNMDAscale = eeNMDAscale * ff/ee;
	}
	if (argc > 2 && *argv[2] != '_') 
		GABAAscale = atof(argv[2])/act*ratio;
	if (argc > 3 && *argv[3] != '_')
		N_prf = size_t(atof(argv[3])*N_A);
	assert(Nstims*N_prf<=N_E);
	assert(Nstims*N_prf<=N_in);
	assert(Nstims*(t_gamma-1)+t_rest+t_gamma/2<t_max);
	if (argc > 4 && *argv[4] != '_')
		GABABscale = atof(argv[4])/act*ratio;
	if (argc > 5 && *argv[5] != '_')
		 gamma_prime = atoi(argv[5]);
	
	srand ( time(NULL) );					// initialize rand()
	
	
	/******* INITIALIZE VARIABLES *********/
	size_t i, j, cell, Ns[4], tstep, N_oth;   // changed to Ns[4] from Ns[5] I don't know why it was Ns[5]
	double t, stats[7], spiketime, off, PING_AMPA, PING_NMDA, PING_GABA=0, I_exts[N_A]={0.0}, dr;
//	string filenames[] = {"results/array1.dsv","results/array2.dsv","results/array3.dsv"};
	valarray<int> spikes_A(N_A*Nnetworks*Ntsteps);	
	valarray<double> w_in(N_in*N_A), syn_in(N_in), temp1;
	valarray<double> allsyn_ff0(0.0,N_E*2*Ntsteps), syn_ff[Nnetworks];
	for (j=0; j<Nnetworks; j++)
		syn_ff[j].resize(N_E*2,0.0);
	slice tinslice, s, tIslice, tffslice;
	// valarray<int> Ispikes(N_I);
	size_t start=0, alengths[] = {Ntsteps,N_prf}, astrides[] = {N_in,1};
	valarray<size_t> lengths(alengths,2), strides(astrides,2), lengths2(lengths), strides2(strides);
	lengths2[1] = size_t(N_prf*3/2), strides2[0] = N_E*2;
	gslice AMPAsyns(0,lengths2,strides2),NMDAsyns(N_E,lengths2,strides2), othcells;
	valarray<double> in_act(N_in*Ntsteps), taus(N_in), rate(N_in);
	gslice actcells(start,lengths,strides);			//input cells that fire with stimulus
	
	
	Neuron Gamma_Inh("inh",input,inhibition);
	Network allnetworks[Nnetworks];
	//	w_in = loadvalarray("w_in.dsv", N_in,N_A);
	//	outputvalarray(w_in, N_in, N_A,"w_in.dsv");
	for (j=0; j<Nnetworks; j++) 
		allnetworks[j].initialize(N_A, N_in, connectivity, input, inhibition);
	allnetworks[0].fillwsize(Ns); assert(N_E == Ns[1]); assert(N_I == Ns[2]);
	allnetworks[0].allneurons[0].fillwstats(stats);
	/*	
	 Network Buffer1(N_A, w_in, connectivity, input, inhibition);
	 Buffer1.fillwsize(Ns); assert(N_E == Ns[1]); assert(N_I == Ns[2]);
	 Network Buffer2(N_A, w_in, connectivity, input, inhibition);
	 */
#ifdef __APPLE__
	openfiles();	// for writing data. on cluster, __APPLE__ is not defined, so the files are not opened
#endif
	if (!funcoutfile1.is_open()) 
		funcout = 0;		// if funcoutfile isn't open (eg. we're on cluster), don't write to it
	
	
	
	
	/******* INITIALIZE EXTERNAL INPUT *********/
	// external input to all buffers
	const bool spikinginput=0;
	if (spikinginput) {
		for (j=0; j<Nstims; j++) {
			t_stims[j] = t_gamma*(j+1)-t_length + dt;		// changed t_rest to dt
			t_ends[j] = t_gamma*(j+1) + dt;
			actcells = gslice(N_prf*j,lengths,strides);
			in_act[actcells] = initinput(in_act[actcells], f_in_off, t_stims[j], f_in_on, t_ends[j]);
		}
		N_oth = N_in - N_prf*Nstims;
		lengths2[0] = Ntsteps; lengths2[1] = N_oth; strides2[0] = N_in; strides2[1] = 1;
		othcells = gslice(N_prf*Nstims,lengths2,strides2);
		in_act[othcells] = initinput(in_act[othcells], f_in_off, t_stims[0], f_in_off, t_ends[0]);
		/*	i=1; s = slice(N_in,N_in/i,i);
		 in_act[s] = 1;		*/ // get initial interneuron inhibition (now implemented in Neuron::initialize())
		/*
		 lengths2[0] = size_t(t_max/t_gamma); lengths2[1] = size_t(N_in/2); strides2[0] = size_t(t_gamma/dt)*N_in; strides2[1] = 2;
		 gslice gs = gslice(size_t(2/dt)*N_in, lengths2, strides2);
		 in_act[gs] = 1;	// everyone spikes first time step
		 */	// simulate gamma inhibition
		in_act = spikes2act(in_act, 0.9, tau_AMPA,1);
	}
	else {
		bool gammamodulated=1;
		for (cell = 0; cell<N_in; cell++) {
			if (gammamodulated) {
				taus[cell] = floor(cell/N_prf+1)*t_gamma+normrand(0,t_gamma*0.2)*1;
			}
			else {
				taus[cell] = doubrand(2*t_gamma/4,Nstims*t_gamma+2*t_gamma/4);
			}
			taus[cell] = floor(taus[cell]/dt)*dt;
			rate[cell] = normrand(30,5);
		}
		taus = valarraysort(taus);
		outputvalarray(taus, N_in, 1,"taus.dsv");
	}		// rate-based input where taus gives the timing of the sharp event 
	
	
	// feedforward input to first buffer
	if (1 && PING) {
		spiketime = doubrand(0, 5);
		for (i=0; i<N_prf; i++) {
			for (off = t_gamma*gamma_prime + t_rest; off < t_max-7; off += t_gamma) {
				allsyn_ff0[2*N_E*size_t((off+spiketime)/dt) + i]  = 1;				// spike in AMPAsyns
				allsyn_ff0[2*N_E*size_t((off+spiketime)/dt) + i+N_E]  = 1;		// spike in NMDAsyns
			}
		}
	}
	else {	
		for (i=0; i<N_prf; i++) {
			double duty = 15;
			spiketime = doubrand(0, duty);
			for (off = t_gamma*gamma_prime + t_rest; off < t_max-duty; off += 10) {
				allsyn_ff0[2*N_E*size_t((off+spiketime)/dt) + i]  = 1;				// spike in AMPAsyns
				allsyn_ff0[2*N_E*size_t((off+spiketime)/dt) + i+N_E]  = 1;		// spike in NMDAsyns
			}
		}
	}		// make the input to the first buffer of the same form as the expected activity of each buffer
	outputvalarray(allsyn_ff0, N_E*2, Ntsteps,"results/inspikes.dsv");
	allsyn_ff0[AMPAsyns] = spikes2act(allsyn_ff0[AMPAsyns], 0.9, tau_AMPA,0);
	allsyn_ff0[NMDAsyns] = spikes2act(allsyn_ff0[NMDAsyns], 0.5, tau_NMDA,0);
	
	
	
	/******* SIMULATION LOOP *********/
	for (tstep = 0; tstep < Ntsteps; tstep++) {
		t = tstep*dt;
		if (spikinginput) {
			tinslice = slice(tstep*N_in, N_in, 1);
			syn_in = in_act[tinslice];					// external input for this timestep
		}
		else {
			for (cell = 0; cell<N_in; cell++) {
				double fr_sharp=200, tau_sharp=10, baseline=20;
				// 150 Hz is peak firing rate during sharp event, 15 ms is length of sharp event, and 20 is baseline rate
				dr = normrand(0,5.0)/sqrt(dt) + fr_sharp/dt*(t==taus[cell]) - (rate[cell]-baseline)/tau_sharp;
				rate[cell] = rectify(rate[cell] + dr*dt);
				syn_in[cell] = 0.00175*rate[cell];
			}
			outputvalarray(rate, 1, N_in,"results/OBrates.dsv");
		}
		tffslice = slice(tstep*N_E*2, N_E*2, 1);
		syn_ff[0] = allsyn_ff0[tffslice];
		PING_AMPA = 0, PING_NMDA = 0;
		
		
		
		/*		if (floor(t)==1) {for (i=0; i<N_prf; i++) {I_exts[i] = 100;}}
		 else if (floor(t)==150) {for (i=0; i<N_prf; i++) {I_exts[i] = -75;}}
		 else if (floor(t)==50 || floor(t)==200) {for (i=0; i<N_prf; i++) {I_exts[i] = 0;}}
		 */	// I_exts for this timestep	

		if (AV_sfile.is_open()) {						// if we are recording data
			AV_sfile << tstep*dt << " ";				
			AV_dfile << tstep*dt << " ";
			Atransmitfile << tstep*dt << " ";
			if (funcout) {
				funcoutfile1 << tstep*dt << " ";
				if (funcout ==2) 
					funcoutfile2 << tstep*dt << " ";
			}
		}			// start output file lines with time
		
		for (j=0; j<Nnetworks; j++) {
			allnetworks[j].step(syn_in, syn_ff[j], I_exts, PING_GABA*PING_GABAscale);
			for (i=0; i<N_A; i++) {		
				allnetworks[j].allneurons[i].fillwstats(stats);
				spikes_A[tstep*Nnetworks*N_A + j*N_A + i] = (stats[2] == 1);
				if (i<N_E) {
					PING_AMPA += stats[3];	PING_NMDA += stats[4];
					if (j<Nnetworks-1) {
						syn_ff[j+1][i] = stats[3];			
						syn_ff[j+1][i+N_E] = stats[4];	// 0:N_E-1 have AMPA and N_E:2*N_E-1 have NMDA
						//	outputvalarray(syn_ff[j], 1,N_E*2,filenames[j]); 
					}
				}
				if (AV_sfile.is_open()) 				// if we are recording data
					record(stats,i,j,t);
			}	// get stats of all neurons in network
		}		// simulate timestep
		if (PING) {
			Gamma_Inh.step(PING_AMPA*PING_AMPAscale, PING_NMDA*PING_NMDAscale, 0, 0, 0, 0, 0);
			Gamma_Inh.fillwstats(stats);
			PING_GABA = stats[5];
			if (AV_sfile.is_open()) 				// if we are recording data
				record(stats,25,Nnetworks,t);
		}
		
		
		if (AV_sfile.is_open()) {					
			AV_sfile << endl;
			AV_dfile << endl;
			Atransmitfile << endl;
			if (funcout) {
				funcoutfile1 << endl;
				if (funcout == 2) 
					funcoutfile2 << endl;
			}			
		}			// if we are recording data, end lines
		
	}			// end of timestep loop
	
	
	
	if (AV_sfile.is_open()) {		
		closefiles();
	}		// close all files
	
	
	
	/******* WRITE TO SWEEPFILE *********/
	if (1 || argc >4) {
		valarray<int> idealact(0,N_E*Nnetworks), trueact(0,N_E*Nnetworks), actoverlap(0,N_E*Nnetworks);
		valarray<int> celllast50spikes(size_t(50/dt));
		slice celllast50slice;
		size_t net;
		
		for (i=0; i<taus.size(); i++) {
			net = size_t( (taus[i] - t_gamma/2)/t_gamma );	// which buffer should store this sharp event
			idealact[net*N_E+i] = 1;					// which cell should store this sharp event
		}		// determine ideal activity
		for (j=0; j<Nnetworks; j++) {
			for (i=0; i<N_E; i++) {
				allnetworks[j].allneurons[i].fillwstats(stats);
				if (!isnan(stats[0])) {
					celllast50slice = slice(size_t((t_max-50)/dt)*Nnetworks*N_A + j*N_A + i, size_t(50/dt)-2, Nnetworks*N_A);
					celllast50spikes = spikes_A[celllast50slice];
					if (celllast50spikes.sum() > 1) {
						trueact[j*N_E + i] = 1;
					}
					else {
						trueact[j*N_E + i] = 0;
					}
				}
				else {
					trueact[j*N_E + i] = -1;
				}// trueact only catalogs activity of excitatory cells
			}
		}		// determine true activity
		actoverlap[idealact==trueact] = 1;
		
		
		// calculate the gamma modulation index by comparing the phase distribution of sharp event onsets with
		// a cosine wave with peak at the center of the cycle
		valarray<double> phase = (taus - t_gamma/2), moving(taus.size()), osc(taus.size());
		valarray<bool> bool1(taus.size()), bool2(taus.size()), bool3(taus.size());
		valarray<int> bincount(taus.size());
		size_t Nbins = phase.size();
		double phasebin = 1.0/20, modindex, low, high;

		for (i = 0; i<phase.size(); i++) {
			phase[i] = (phase[i] - floor(phase[i]/t_gamma)*t_gamma )/t_gamma;	// phase is in [0,1]
			osc[i] = double(i)/osc.size()*2*pi;					// osc is 0:2*pi/osc.size:2*pi
		}
		phase = valarraysort(phase);
		for (i=0; i<Nbins; i++) {
			low = double(i)/Nbins-phasebin; high = double(i)/Nbins+phasebin;
			bincount = 0;
			bool1 = phase>low;
			bool2 = phase<high;
			if (low<0)
				bool2 = bool2 + ((phase - 1.0)>low);
			else if (high>1)
				bool1 = bool1 + ((phase + 1.0)<high);
			bincount[bool1 && bool2] = 1;
			moving[i] = bincount.sum();
		}
		moving = moving/(2*phasebin*Nbins);
		moving = (moving - moving.sum()/moving.size())*-cos(osc)*2/Nbins;
		modindex = moving.sum();
		
		valarray<double> phasediff(phase);
		double modindex1, modindex2;
		phasediff = pow(phase - 0.5, 2);	
		modindex1 = sqrt(phasediff.sum()/phasediff.size());	// deviation of phase from t_gamma/2 (i.e. phase=0.5)
		phasediff = pow(phase - phase.sum()/phase.size(), 2);
		modindex2 = sqrt(phasediff.sum()/phasediff.size());	// standard deviation of phase



		
		outputvalarray(idealact,actoverlap.size(),1,"results/activities1.dsv");
		outputvalarray(trueact,actoverlap.size(),1,"results/activities2.dsv");
		outputvalarray(actoverlap,actoverlap.size(),1,"results/activities3.dsv");

		sweepfile.open("./results/sweepresults.dsv", ios::app);
		sweepfile << (eeNMDAscale+ffNMDAscale)*act << "\t" << GABAAscale*act/ratio << "\t";
		sweepfile << double(N_prf)/N_A << "\t" << GABABscale*act/ratio << "\t";
		
		allnetworks[0].allneurons[size_t((N_A + N_E)/2)].fillwstats(stats);
		sweepfile << setprecision(3) << stats[6] << "\t";
		allnetworks[Nnetworks-1].allneurons[size_t((N_A + N_E)/2)].fillwstats(stats);
		sweepfile << setprecision(3) << stats[6] << "\t";
		
		sweepfile << setprecision(3) << idealact.sum() << "\t" ;
		sweepfile << setprecision(3) << trueact.sum() << "\t";
		sweepfile << setprecision(3) << actoverlap.sum() << "\t";
		sweepfile << setprecision(3) << 1 - ( int(N_E*Nnetworks) - actoverlap.sum() )/double(idealact.sum()) << "\t";
		sweepfile << setprecision(3) << modindex << "\t" << modindex1 << "\t" << modindex2;
		sweepfile << endl;
		
		sweepfile.close();
		
	}		// end of sweepfile write

	
	return 0;
}		// end of main()










/******* HELPER FUNCTIONS *********/


valarray<double> initinput(gslice_array<double> ginput, 
									double f_base, double t_on, double f_on, double t_off
									) {	// fills a valarray with spikes given firing parameters for those cells
	valarray<double> tmp(ginput);
	
	size_t i, N = tmp.size()/Ntsteps;
	
	if (t_on>=t_max)
		t_on = t_max-dt;
	slice oslice(0, N, 1), baseslice(0,N*size_t(t_on/dt),1); 
	slice onslice(N*size_t(t_on/dt),N*size_t((t_off-t_on)/dt),1);
	slice offslice(N*size_t(t_off/dt),N*size_t((t_max-t_off)/dt),1);
	valarray<double> basein(tmp[baseslice]), onin(tmp[onslice]), offin(tmp[offslice]);
	valarray<bool> pos(tmp.size());							// spike times
	tmp[baseslice] = connect(basein, f_base * dt/1000);		// baseline spiking pattern
	tmp[onslice] = connect(onin, f_on * dt/1000);			// stimulus on spiking pattern
	tmp[offslice] = connect(offin, f_base * dt/1000);		// stimulus off spiking pattern
	for (i=0; i<pos.size(); i++) {pos[i] = (tmp[i]>0);}
	tmp[pos] = 1;						// spike times
	tmp[oslice] = 0;					// no spikes in first time step
	
	
	return tmp;
}


valarray<double> spikes2act(valarray<double> spikes, double alpha_x, double tau, bool w) {
	// bool w tells whether to write the spikes of this input to inspikesfile
	size_t i, t, N = spikes.size()/Ntsteps, spike;
	for (i=0; i<N; i++){				//turn spikes into synaptic activation
		for (t=1; t<Ntsteps; t++){
			spike = int(spikes[i+t*N]);
			//if (inspikesfile.is_open() && w) {inspikesfile << spike << " ";}
			spikes[i+t*N] = spikes[i+(t-1)*N] * exp(-dt/tau);
			if (spike == 1) {
				spikes[i+t*N] = spikes[i+t*N]+alpha_x*(1-spikes[i+t*N]);
				if (inspikesfile.is_open() && w) {inspikesfile << t*dt << "\t" << i << endl;}
			}
		}
		if (inspikesfile.is_open() && w) {inspikesfile << endl;}
	}
	return spikes;
}


void openfiles() {
	AV_sfile.open("./results/AV_s.dsv");
	AV_dfile.open("./results/AV_d.dsv");
	Atransmitfile.open("./results/Atransmit.dsv");
	Aspikecoordsfile.open("./results/Aspikecoords.dsv");
	if (funcout) {	// if we are writing to funcoutfile
		funcoutfile1.open("./results/funcoutfile1.dsv");
		funcoutfile1 << "Currents\nTime (ms)\tAMPA\tNMDA\tGABAA\tGABAB\n";
		if(funcout==2){
			funcoutfile2.open("./results/funcoutfile2.dsv");
			funcoutfile2 << "Conductances\nTime (ms)\tAMPA\tNMDA\tGABAA\tGABAB\n";
		}			
	}
	outfile.open("./results/outfile.dsv");
	inspikesfile.open("./results/inspikes.dsv");	//initialize input
}	



void record(double stats[7], size_t i, size_t j, double t) {
	if (stats[2] == 1)
		Aspikecoordsfile << t << "\t" << i+j*N_A << endl;
	AV_sfile << stats[0] << " ";
	AV_dfile << stats[1] << " ";
	if (i<N_E && j<Nnetworks) {Atransmitfile << stats[3] << " " << stats[4] << " ";}
	else {Atransmitfile << stats[5] << " " << stats[6] << " ";}
}


void closefiles() {
	AV_sfile.close();
	AV_dfile.close();
	Atransmitfile.close();
	if (funcout) {
		funcoutfile1.close();
		if (funcout==2) 
			funcoutfile2.close();
	}
	Aspikecoordsfile.close();
	outfile.close();
	inspikesfile.close();
}	







