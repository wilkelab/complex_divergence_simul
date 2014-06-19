complex_divergence_simul
========================

The PPI_Marcotte_Setup.tar.bz2 file is an intermittent backup for ease of email.  It can mostly be ignored.

The main setup folder for running the evolutionary simulation is the PPI_Marcotte_Setup directory.  To start a simulation, the whole directory can simply be copied as it stands.  Once the directory is completely copied, a simulation can be started with the command 'python main.py'.  The default number of mutations to test is 1000.

The code depends on many things.  First, biopython must be installed and importtable with the 'import Bio' command.  Next the location of foldx is currently hardcoded as /home/austin/local/bin/foldx64Linux.  Of course, that will need to be changed for your location.

Also of note, the BuildModel command in Foldx appears to randomly throw segmentation faults.  The overcome this fact, I catch the faults until the program returns a successful exit.  Obviously, I am assuming a successful exit means the program work correct; this assumption might not be true.
