Programs AFTPOI and RAFTPOI

Estimation of parameter values in the Omori-Utsu Law for the decay rate of aftershocks, which has been called Ågthe modified Omori LawÅh by Tokuji Utsu. 

1. The Poisson process of the Omori-Utsu law and residual Poisson process

1.1 The MLE computation

This FORTRAN program computes the maximum likelihood estimates (MLEs) of the parameter values in the Omori-Utsu formula (the modified Omori formula; Omori, 1894; Utsu, 1961) which represent the decay law of aftershock activity in time (see also Utsu et al., 1995). In these equations, f(t) represents the rate of aftershock occurrence at time t, where t is the time measured from the origin time of the main shock. \mu, K, c and p are non-negative constants. \mu represents constant-rate background seismicity which may be included in the aftershock data.

f(t) = \mu + K / (t + c)^p

The FORTRAN program was written by Ogata (1983), where the negative log-likelihood function is minimized by the Davidon-Fletcher-Powell algorithm. Namely, starting from a given set of initial guess of the parameters that is given in the control file [aftpoi.open], the program [aftpoi.f] repeats calculations of the function values and its gradients at each step of parameter vector. At each cycle of iteration, the linearly searched step (lambda), negative log-likelihood value (-LL), and two estimates of square sum of gradients are shown. At each linear search of the likelihood computation, -LL and the 4 parameter values are shown. The value -LL decreases and tends to a final value. one of the square sum of gradients becomes nearly zero and the iteration is terminated. 

In the present case, the following results are displayed, where AIC = -2 LL + 2 x (number of searched variables), and the number = 4 in this case. The calculated record of the program [aftpoi] is stored by the name [aftpoi.print] in the directory of [Calculations] for your initial check of the program. The program [aftpoi] always call the dataset named [work.etas] with the format as described below.


1.2 The residual point process

The FORTRAN program [raftpoi.f] computes the following output for displaying the goodness-of-fit of the Omori-Utsu model to the data. The cumulative number of earthquakes at time t since t0 is given by the integration of f(t) in Func6 with respect to the time t,

F(t) = \mu (t - t0) + K {c^(1-p) - (t - t_i + c)^(1-p)} / (p - 1)

where the summation of i is taken for all data event. The output of the program is given in the file [work.res] which includes the column of {F(t_i), i = 1, 2, Åc, N}. This file is used for displaying the cumulative curve and magnitude v.s. transformed time F(t_i). For the users who have the free graphic software R, I set a module [r.raftpoi] to display them, in addition to the module [r.seis] for displaying the cumulative curve and magnitude v.s. ordinary time t_i. If the observed rate of occurrence is compared with the calculated one from the model, period of decreased or increased seismic activity (relative quiescence or activation) can be recognized. The calculated record of the program [raftpoi] is stored by the name [raftpoi.print] in the directory of [Calculations] for you to check the program. The program [raftpoi] always call the dataset named [work.res] with the format as described below.



2. Data [work.etas]

The data file for this program is a sequential file with a 9 columns-format of the sequential number, longitude, latitude, magnitude, time from the mainshock in days, depth, year, month, and day, as is given by the prototype data file [work.etas]. The first row of the data is the title. The second row corresponds to the main shock (not used in this program, but is used in other programs such as [etas.f]). The time (usually in days) is measured from the main shock (t = 0). If aftershock in an early stage of the sequence (e.g., from t = 0 to 0.01 day) are not included because of incomplete observation, it is better to set Tstart at the beginning of the aftershock data (e.g., Tstart = 0.01). 0 < Tstart < t_1 < . . . < t_i < t_{i+1} < . . . < t_n < Tend in 5th column of [work.etas]. Magnitudes in 4th column are used only in selecting the aftershocks by giving a threshold magnitude. 


3. Control file [aftpoi.open]

The file [aftpoi.open] provides input variables such as the restricting range of data and initial values of parameters. 

The first row of [aftpoi.open] consists of two numbers; the first number is always set 6 for the selection of the subroutine [func6] for the Omori-Utsu model, and then for the second number, any dummy number can be set.

The second row consists of three numbers; zts, zte and tarst in the above illustrated variables when the file is applied to the program [aftpoi.f]. When the file is applied to the program [raftpoi.f] for the diagnostic outputs, the last number can be ztend (>= zte) that is the last time of the period of available data. The example shows the first 6.5days of aftershock period since the mainshock due to [work.etas], where the target interval is [tarst, zte] days.

The third row consists of two numbers; threshold magnitude Mth and reference magnitude Mz. Here I have taken the magnitude of the mainshock as the reference magnitude, but you may take Mz = Mth for general seismicity. Note that the MLE of the parameter K differs depending on the reference parameter. 

The fourth row provides the initial estimates of the five parameters \mu, K, c, \alpha and p when this is applied the program [aftpoi.f], but \alpha (= 0.0) is dummy parameter to make same format as [etas.open]. If you set \mu = 0.0 in the initial value, this variable remains 0.0 throughout the calculations. Also, if you set p = 1.0, this variable remains 1.0 (the original Omori formula) throughout the calculations. When this control file is used for the program [retas.f] to make the diagnostic outputs, these values in this row have to be the MLEs estimated by the program [aftpoi.f]. 

The fifth row consists of three numbers; zts, zte and tarst only used when the file is applied to the program [raftpoi.f] for the diagnostic outputs, the last number must be zte that is the end of the target period. 

The sixth row is the same as the fourth row, only used for the program [raftpoi.f], and the seventh row of the value of the negative of the maximum likelihood corresponding to the MLEs, used just for the user( note to compare goodness-of-fit with other cases or models.

4. Remarks

4.1 Initial values for parameters in control file [aftpoi.open]

The program works for initial values given fairly arbitrary. A suggested set of the initial values could be set in the control file [aftpoi.open]. Type an appropriate value if you wish to change it. Unreasonable initial values may cause an overflow error. If -LL value and gradients do not seem to converge, we recommend you to use the last parameter estimates as the initial values in the control file [etas.open] and restart the execution.

4.2 Abnormal termination

If the data are very different from the proposed model, an overflow error may occur or the solution may not converge within the fixed number of iterations (30N, N is the number of parameters). The results in this case are often incorrect and not recommended to adopt as the estimates. 

4.3 Standard errors 

The standard errors (SD) for estimates are computed by taking the square root of the corresponding diagonal element of the Inverse Fisher matrix. 

5. Graphs

The R language modules [r.seis], [r.seispoi] and [r.raftpoi] are for the graphs of cumulative number and magnitude of earthquakes against the ordinary time and transformed time, respectively: These modules show the files [work.etas] and [work.res. The last graph includes printed number with certain time values for vertical dotted lines showing target interval and MLEs in control file [aftpoi.open]. The outputs of the graphs in postscript format, named as seis.ps, seispoi.ps, and raftpoi.ps are given in the same directory. 

The R language modules [r.raftpoi] and [r.seisaftpoi] shows the theoretical curve of cumulative frequency (red color) which is a straight line by definition of the transformed time, and the cumulative frequency curves of the data (red color). 

The graphical free software R is also installed in the other directory based on the Microsoft Windows XP (RGui). To command in the R console, write 

 > source('r.*') 

where * correspond above R-module, then the graphical window should appear and draw appropriate figures.

The FORTRAN programs are originally designed and programmed (1983) and reprogrammed (December 2005) by Yosihiko Ogata, Institute of Statistical Mathematics, Tokyo, Japan. The R-module was designed and programmed by Yosihiko Ogata (2003).

References                

Akaike, H. (1974). A new look at the statistical model identification, IEEE Trans. Automat. Control, AC-19, 716-723.

Ogata, Y., Estimation of parameters in the modified Omori formula for aftershock frequencies by the maximum likelihood procedure, J. Phys. Earth, 31, 115-24, 1983. 

Omori, F., On the after-shocks of earthquakes, J. Coll. Sci. Imp. Univ. Tokyo, 7, 111-200, 1894.

Utsu, T., A statistical study on the occurrence of aftershocks, Geophys. Mag., 30, 521-605, 1961.

Utsu, T., Y. Ogata, and R. S. Matsu'ura, The centenary of the Omori formula for a decay law of after-shock activity, J. Phys. Earth, 43, 1-33, 1995.
