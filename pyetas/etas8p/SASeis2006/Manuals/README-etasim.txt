Program ETASIM

These programs produce simulated data files for given sets of parameters in the point process model used in ETAS. See Ogata (1981) or Ogata (1998) for theoretical basis. It is noted that the intensity defined by a combination of parameter values should be non-negative and well-defined. Due to some combinations of parameter values, the simulated data can be explosive. The simulated events are given in the file [work.etasim] in the same format as [work.etas].

Control file [etasim.open]. 

The parameters of the ETAS model and other necessary parameters should be given in the control file [etasim.open]. There are two options; either simulating magnitude by Gutenberg-Richter's Law (ic = 0) or using magnitudes from [work.etas] (ic = any number except 0, say 1). For the first option, you have to provide b-value of G-R law and number of events to be simulated. The second option simulates the same number of events that are not less than threshold magnitude in the data [work.etas], and simulation starts after a precursory period depending on the same history of events in [work.etas] in the period. Therefore the occurrence times of simulated events in [work.etasim] is the same as those in [work.etas] during the precursory period.

The first row consists of two numbers, the first number is ic and the second number is b-value. If ic is not 0 second number can be any integer (dummy number). 

The second row consists of two numbers, the first number is tstart where [0, tstart] is precursory period. The second number is the number of the simulated events if ic = 0, otherwise dummy number.
 
The third row consists of two numbers, the first number is cutoff (threshold) magnitude of the simulated data, and the second number is the reference magnitude to calculate the intensity rate of the ETAS model. This should be the same reference magnitude given in [etas.open]. For example, this could be the magnitude of the mainshock in case of aftershock simulation.

The fourth row consists of five numbers of the ETAS parameters, which are \mu, K, c, \alpha and p in the order from the first to the fifth number.

Calculated record of the program [etasim] is stored by the name [etasim.print] in the directory of [Calculations] for your initial check of the program.

The FORTRAN programs are originally designed and programmed (1985) by Yosihiko Ogata, Institute of Statistical Mathematics, Tokyo, Japan. The R-module was designed and programmed by Yosihiko Ogata (2003).


References

Ogata, Y, On Lewis' simulation method for point processes, IEEE Information Theory, IT-27, 23-31, 1981. 

Ogata, Y. (1998) Space-time point-process models for earthquake occurrences, Ann. Inst. Statist. Math., 50, 379-402.
