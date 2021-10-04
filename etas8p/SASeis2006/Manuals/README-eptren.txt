Overview of the program [eptren]

This program computes the maximum likelihood estimates (MLEs) of the coefficients A1, A2, …, An in an exponential polynomial

       f(t) = exp{A_1 + A_2 t + A_3 t^2 + … }                                             (1)

or A_1, A_2, B_2, … , A_n, B_n in  a Poisson process model with an intensity taking the form of an exponential Fourier series


f(t) = exp{A_1 +A_2 cos(2πt/P) +B_2 sin(2πt/P) +A_3 cos(4πt/P) +B_3 sin(4πt/P)+ … }   (2) 
  
which represents the time varying rate of occurrence (intensity function) of earthquakes in a region. A non-stationary Poisson process is assumed for the model of earthquake occurrences. The optimal order n can be determined by minimize the value of the Akaike Information Criterion (AIC) values (cf., Akaike, 1974). 

The maximum number of coefficients is limited to 21 due to the dimension set in the FORTRAN source. The coefficients A_1, A_2, ... , A_{21} can be determined for equation (1) and A_1, A_2, B_2, … , A_{10}, B_{10} can be determined for equation (2) at the most.

Since a non-stationary Poisson process is assumed in the data, we need careful interpretation of the results when the model is applied to data containing earthquake clusters such as aftershocks, foreshocks, and earthquake swarms.

R-modules are also attached in order display the outputs of numerical data for plotting graphs of the estimated functions.

Structure of the program is the following.

          [eptren]
             |------[inputs]
             |------[reduc1]
             |------[reduc2]
             |------[davidn]----------------[funct]
             |         |--------[hesian]
             |         |--------[linear]----[funct]
             |
             |------[fincal]
             |------[output]---[printr]-----[trenfn]
                                  |---------[cyclfn]

1. The name of the analyzing data and its format is either [work.etas] or [work.res] and their format. The choice of either data is controlled by the following control file. The range of x-coordinate is limited to the non-negative valued variable. When dataset [work.etas] includes remarkable clusters, the interpretation of the trend and cycles has to be careful. Rather try [work.res] due to the time transformation by the ETAS model.

2. Control file including the selection of either trend or cycle fitting is [eptren.open] in which necessary variables to read are given and explained.

3. Output files of [eptern] are [out.eptren1] and [out.eptren2]. To make figure R-module [r. eptrend] or [r. epcycle] is used to make R-display and also produces postscript output [eptrend.ps] and [epcycle.ps].

4. Calculated record of the program [eptren] is stored by the name [eptren??.print] in the directory of [Calculations] for your initial check of the program, where each "?”is either 0, 1 or 2, depending on the same number given in [eptren.open].

This program was originally designed (January 1985) and revised (December 2005) by Yosihiko Ogata, and programmed and also reprogrammed by Koichi Katsura, Institute of Statistical Mathematics, Tokyo, Japan. The R-module was designed and programmed by Jiancang Zhuang and Ogata, the Institute of Statistical Mathematics, Tokyo, Japan (December 2005).

References

Akaike, H. (1974). A new look at the statistical model identification, IEEE Trans. Automat. Control, AC-19, 716-723.

Lewis, P.A.W. and G.S. Shedler (1976). Statistical analysis of non-stationary series in a database system, IBM. J. Res. Develop., 20, 465-481.

Maclean, C.J. (1974). Estimation and testing of an exponential polynomial rate function within the non-stationary Poisson process, Biometrika, 61, 81-86.   
           

