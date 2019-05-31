#include "network.h"
//#include "newnetwork.h"
const size_t N_A = 400, N_E = size_t(N_A*ratio/(ratio+1)), N_I = N_A - N_E;
size_t N_in = N_E, N_prf = size_t(N_A*0.4);

bool EIF = 0, HH = !EIF;
int funcout = 2;  // what to funcoutfile? 
// 0--nothing, 1--currents, 2--currents & (non v-dependent) conductances

const double t_max = 250;
const size_t Ntsteps = size_t(t_max/dt);
const double t_stim = 0;				//time stimulus goes on
const double t_end = 100;				//time stimulus goes off
const double f_in = 0;				//input firing rate
const double f_in_on = 200;				//input firing rate for active cells

const double connectivity = 1;
double act = double(N_E*connectivity);
const char input = 'd', inhibition = 'd';


double eeNMDAscale = 70.0/act, eiNMDAscale = 3.0/act;
double eeAMPAscale= eeNMDAscale/2, eiAMPAscale=eeAMPAscale/8;//5.0/act;
double GABAAscale=7.0/act*ratio, GABABscale=500.0/act*ratio;
double inscale = 20.0, iniscale = eiAMPAscale/2.0;
double ffAMPAscale = 0, ffNMDAscale = 0;
/*
// from Fellous: 
double eeAMPAscale=0.4/act, eiAMPAscale=0.1/act, eeNMDAscale = 30.0/act, eiNMDAscale=3.0/act;			
double GABAAscale = 1.0/act, GABABscale=0;
double inscale = 32.0;
*/

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
		eeNMDAscale = atof(argv[1])/act;	
		//eiNMDAscale = eeNMDAscale/10;
		eeAMPAscale = eeNMDAscale/2;	// comment if you don't want AMPA to scale with NMDA
//		inscale = max(inscale, eeAMPAscale*act);
		eiAMPAscale = eeAMPAscale/6; //min(eeAMPAscale*4,inscale/act);
	}
	if (argc > 2 && *argv[2] != '_') 
		GABAAscale = atof(argv[2])/act*ratio;
	if (argc > 3 && *argv[3] != '_')
		N_prf = size_t(atof(argv[3])*N_A);
	if (argc > 4 && *argv[4] != '_')
		GABABscale = atof(argv[4])/act*ratio;
	
	srand ( time(NULL) );					// initialize rand()
	
	
	/******* INITIALIZE VARIABLES *********/
	size_t i, Ns[4], tstep;
	double t, stats[7], I_exts[N_A]={0.0};					// I_exts[N_A] is hardcoded here
	valarray<int> spikes_A(N_A*Ntsteps);	
	valarray<double> in_act(N_in*Ntsteps), syn_in(N_in);
	slice tinslice, s, tIslice;
	// valarray<int> Ispikes(N_I);
	size_t start=0, alengths[] = {Ntsteps,N_prf}, astrides[] = {N_in,1};
	valarray<size_t> lengths(alengths,2), strides(astrides,2);
	valarray<size_t> lengths2(lengths), strides2(strides);
	gslice actcells(start,lengths,strides);			//input cells that fire with stimulus
	start++; start=start-1;
	
	start = N_prf; lengths2[1]= N_in-N_prf;
	gslice othcells(start,lengths2,strides);		 //input cells that don't fire with stimulus
	
	Network WorkingMemory(N_A, N_in, connectivity, input, inhibition);
	WorkingMemory.fillwsize(Ns); assert(N_E == Ns[1]); assert(N_I == Ns[2]);
	
	

	
#ifdef __APPLE__
	openfiles();	// for writing data. on cluster, __APPLE__ is not defined, so the files are not opened
#endif
	if (!funcoutfile1.is_open()) 
		funcout = 0;		// if funcoutfile isn't open (eg. we're on cluster), don't write to it
	
	
	
	
	/******* INITIALIZE EXTERNAL INPUT *********/	
	in_act[actcells] = initinput(in_act[actcells], f_in, t_stim, f_in_on, t_end);
	in_act[othcells] = initinput(in_act[othcells], f_in, t_stim, f_in, t_end);
	in_act = spikes2act(in_act, alpha, tau_AMPA,0);
	
	
	/******* SIMULATION LOOP *********/
	for (tstep = 0; tstep < Ntsteps; tstep++) {
		t = tstep*dt;
		tinslice = slice(tstep*N_in, N_in, 1);
		syn_in = in_act[tinslice];					// external input for this timestep
		

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
		
		
		// cout << t << " ";
		WorkingMemory.step(syn_in, I_exts);					// timestep
		// cout << endl;
		

		
		for (i=0; i<N_A; i++) {		
			WorkingMemory.allneurons[i].fillwstats(stats);
			spikes_A[tstep*N_A+i] = (stats[2] == 1);
			
			if (AV_sfile.is_open()) {						// if we are recording data
				record(stats,i,0,t);
			}
			
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
	}
	

	/******* WRITE TO SWEEPFILE *********/
	if (argc >4) {
		gslice preslice, onprfslice, onnprfslice, postprfslice, postnprfslice, postinhslice;
		int tpre=int(t_stim/dt), ton=int((t_end-t_stim-10)/dt), tpost=int(min(t_max-t_end-10,50.0)/dt);
		valarray<int> prespikes(N_E*tpre);
		valarray<int> onprfspikes(N_prf*ton), onnprfspikes((N_E-N_prf)*ton);
		valarray<int> postprfspikes(N_prf*tpost), postnprfspikes((N_E-N_prf)*tpost), postinhspikes(N_I*tpost);
		double pref, onprff, onnprff, postprff, postnprff, postinh;
		
		
/*		lengths2 = lengths;		//lengths = {Ntsteps,N_prf}
		lengths2[0] = tpre;
		lengths2[1] = N_E;
		strides2 = strides;		// strides = {N_in,1}
		strides2[0] = N_A;
		preslice = gslice(0,lengths2,strides2);
		prespikes = spikes_A[preslice];
		pref = prespikes.sum()*1000/t_stim/N_E;
*/ 
		
						// from 10 ms after the stimulus starts until it ends
		lengths2 = lengths;		//lengths = {Ntsteps,N_prf}
		lengths2[0] = ton;
		lengths2[1] = N_E;
		strides2 = strides;		// strides = {N_in,1}
		strides2[0] = N_A;
		onprfslice = gslice(int((t_stim+10)/dt*N_A),lengths2,strides2);
		onprfspikes = spikes_A[onprfslice];
		lengths2[1] = N_E - N_prf;
		onnprfslice = gslice(int((t_stim+10)/dt*N_A) + N_prf,lengths2,strides2);
		onnprfspikes = spikes_A[onnprfslice];
		onprff = onprfspikes.sum()*1000/(ton*dt)/N_prf;
		onnprff = onnprfspikes.sum()*1000/(ton*dt)/(N_E-N_prf);
		
		lengths2 = lengths;						// last 50 ms of the simulation
		lengths2[0] = tpost;
		postprfslice = gslice(int(t_max/dt-tpost)*N_A,lengths2,strides2);
		postprfspikes = spikes_A[postprfslice];
		lengths2[1] = N_E - N_prf;
		postnprfslice = gslice(int(t_max/dt-tpost)*N_A + N_prf,lengths2,strides2);
		postnprfspikes = spikes_A[postnprfslice];
		postprff = postprfspikes.sum()*1000/(tpost*dt)/N_prf;
		postnprff = postnprfspikes.sum()*1000/(tpost*dt)/(N_E-N_prf);

		lengths2[1]=N_I;
		postinhslice = gslice(int(t_max/dt-tpost)*N_A+N_E,lengths2,strides2);
		postinhspikes = spikes_A[postinhslice];
		postinh = postinhspikes.sum()*1000/(tpost*dt)/(N_I);
		
		sweepfile.open("./results/sweepresults.dsv", ios::app);
		sweepfile << atof(argv[1]) << "\t" << atof(argv[2]) << "\t";
		sweepfile << atof(argv[3]) << "\t" << atof(argv[4]) << "\t";
		
		sweepfile << setprecision(3) << onprff << "\t" << onnprff << "\t";
		sweepfile << setprecision(3) << postprff << "\t" << postnprff << "\t";
		sweepfile << setprecision(3) << postprff/(postnprff+1) << "\t" << postinh;
		WorkingMemory.allneurons[350].fillwstats(stats);
		sweepfile << "\t" << stats[6];
		
		//output by column (all output is on one line):
		//						g_NMDA	g_GABAA		pattern_size	g_GABAB 
		// then stimulus on:	FR_targeted_cells	FR_nontargeted_cells
		// then end of simul:	FR_targeted_cells	FR_nontargeted_cells	Ratio	FR_inh_cells	GABAB_activation
		
		sweepfile << endl;
		
		sweepfile.close();
		
	}
	
    return 0;
}




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
	// this function is copied from stnetwork, where there are multiple subnetworks indexed by j
	// references to j in this function have therefore been removed
	if (stats[2] == 1)
		Aspikecoordsfile << t << "\t" << i << endl;
	AV_sfile << stats[0] << " ";
	AV_dfile << stats[1] << " ";
	if (i<N_E) {Atransmitfile << stats[3] << " " << stats[4] << " ";}
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







