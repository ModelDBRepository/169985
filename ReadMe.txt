These files go with two papers:

Sanders H, Berends M, Major G, Goldman MS, Lisman JE (2013) NMDA and
GABAB (KIR) Conductances: The "Perfect Couple" for Bistability. J
Neurosci 33:424-429.  doi: 10.1523/JNEUROSCI.1854-12.2013
http://www.jneurosci.org/cgi/doi/10.1523/JNEUROSCI.1854-12.2013

Sanders H, Kolterman BE, Shusterman R, Rinberg D, Koulakov A, Lisman J
(2014) A network that performs brute-force conversion of a temporal
sequence to a spatial pattern: relevance to odor recognition. Front
Comput Neurosci 8:1-11.  doi: 10.3389/fncom.2014.00108
http://journal.frontiersin.org/Journal/10.3389/fncom.2014.00108/


The simulations are written in C++, requiring no external libraries.
The plotting is written in gnuplot (http://gnuplot.info).  The wrapper
functions for compiling, running, submitting jobs to cluster, and
plotting are written in BASH.  My condolences to Windows users, you
are on your own.


**** NMDA and GABAB (KIR) Conductances: The "Perfect Couple" for
     Bistability ***

The files necessary for producing the figures from the 2013 paper are
in the folder "NMDA-GABAB_2013".  network.h is a header file
containing constants and function declarations.  It is identical for
both papers.  functions.cpp contains function definitions.  It has two
points of difference between the two papers, pointed out at the top of
the file.  main.cpp contains the simulation loop.  It is obviously
fairly different between the two papers, as they are different
simulations.


To compile and run the simulation with the default parameter values,
type
./crnewnetwork
at the command line.  This will save the results of the simulation to
results/.  It will additionally plot the results if you have gnuplot
installed at /usr/local/bin/gnuplot.

If the crnewnetwork script fails due to platform differences/package
dependencies you can likely compile and run the model with

g++ -g main.cpp functions.cpp -o newnetwork
./newnetwork

You can add up to four command line arguments after ./crnewnetwork.
These command line arguments are the eeNMDAscale (NMDA conductance at
excitatory-excitatory synapses), GABAAscale (the GABAA conductance at
inhibitory-excitatory synapses), N_prf (the fraction of excitatory
cells activated by external input, expressed as a decimal).  The
conductances at each synapse are scaled by number of active cells so
that the total conductance seen by each cell is approximately the same
as the conductance passed as a command line argument.  Also of note is
that the conductances in the model are expressed in uS/mm^2, whereas
in the paper conductances are expressed in mS/cm^2, which leads to the
numbers in the simulation being 10 times larger than the numbers in
the paper (1000/10^2=10).

clusterscripts/ contains scripts for submitting jobs to the Brandeis
HPCC cluster.  For questions about this or about anything else,
contact Honi Sanders at honi@brandeis.edu.
