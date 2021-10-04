Statistical Analysis of Seismicity - updated version (SASeis2006) 

Ogata, Yosihiko
The Institute of Statistical Mathematics

This program package consists of stand-alone programs, AFTPOI, EPTREN, ETAS, ETASIM, PGRAPH, RAFTPOI and RETAS. The source files are written in FORTRAN77, and can be complied the execution files here were compiled by GNU FORTRAN Ågg77Åh on any platforms where GNU FORTRAN is available. the cygwin, and should be executable at least on Windows XP. These programs only output numerical files for graphics. Especially, these can be used in free statistical software R, and wWe also provide the R source modules and for plotting figures in postscript files format in the same directory, written in R, a widely-used free statistical programming language. All the names are corresponding to the names of the FORTRAN sources and the output files. 
Some test data files for these programs are included. They are either simulated data or actual earthquake (aftershocks) data. It is recommended that you first try these data files to understand how each program works.

AFTPOI: Maximum likelihood estimates (MLEs) of parameters in the Omori-Utsu (modified Omori) formula representing for the decay of occurrence rate of aftershocks with time, which is assumed to distributed according toformulated as a the non-stationary Poisson process.

EPTREN: MLEs of parameters of in a non-stationary Poisson process with a rate function being represented by thean exponential polynomial or an exponential Fourier series. The optimal order of the polynomial and series can beare  determined by the minimized AIC value.

ETAS: MLEs of parameters of in the ETAS model for applying the general seismicity and an aftershock sequences. 

ETASIM: Simulation of earthquake data set to test based on the ETAS program or to examine forecasting variabilitymodel.

PGRAPH: To display some calculate some graphical statistics to grasp the elementary statistical features of your dataset as of a point process or its residual point process. . Since these graphical statistics also show the significance lines of the variability to test the stationary Poisson process, it is useful to make the residual analysis in detail for the dataset produced by the programs RAFTPOI and RETAS below. 

RAFTPOI: Occurrence times of earthquakes are transformed using the Omori-Utsu formula with the MLEs estimated by AFTPOI. This is called as the residual data and can be d used for the diagnostic analysis of the model or the detection of anomalies in the aftershock sequence.

RETAS: Occurrence times of earthquakes are transformed using the ETAS model with the MLEs. This is called as the residual data and can be used for the diagnostic analysis or the detection of anomalies in the aftershock sequence.This is called as the residual data and used for the diagnostic analysis of the model or anomalies of the seismicity.

Some technical notes are as follows:
1. If the programs crash with an error message "cannot find cygwin1.dll", copy cygwin1.dll to c:\windows\system32.
2. To install updated version of R for windows, MacOS or some other UNIX based systems, please go to "http://www.r-project.org/".
3. To install most recent version of cygwin, please go to "http://www.cygwin.com/".

Acknowledgement.
I am grateful Jancang Zhuang for his advices on R commands, and also for carefully reading the manuals for SASeis2006.

References. 

Ogata, Y. (1988) Statistical models for earthquake occurrences and residual analysis for point processes, Journal of American Statistical Association, Application, Vol. 83, No. 401, pp. 9-27

Ogata, Y. and Katsura, K. (1985) Fortran programs for Statistical Analysis of Series of Events (SASE) consisting of the programs EPTREN, LINLIN, SIMBVP, LINSIM and PGRAPH included in Time Series and Control Program Package, TIMSAC-84 (joint with H. Akaike, T. Ozaki, M. Ishiguro, G. Kitagawa, Y. Tamura, E. Arahata, K. Katsura and R. Tamura), Computer Science Monograph, No. 22/23, The Institute of Statistical Mathematics, Tokyo, Japan. 

Utsu, T. and Ogata, Y. (1997) Computer program package: Statistical Analysis of point processes for Seismicity, SASeis, IASPEI Software Library for personal computers, the International Association of Seismology and Physics of Earth's Interior in collaboration with the American Seismological Society, Vol. 6, pp. 13-94. 

