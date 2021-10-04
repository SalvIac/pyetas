Overview of the program [pgraph]

This provides the following several graphical outputs for the point process data set.

1. Cumulative numbers of events versus time, and positions of spikes (subroutine cumlat), and output file is [out.pgCumMT]. The R-module [r.pgCumMT] produces these figures in both console and postscript [pgCumMT.ps].

2. Series of counted events in the moving interval with a fixed length. This is a kind of moving average. Dotted lines indicate i x (standard error), i =1, 2, 3, assuming the stationary poisson process ( subroutine count1 or count2 ), and output file is [out.pgPTnum]. The R-module [r.pgPTnum] produces these figures in both console and postscript [pgPTnum.ps].

3. The log survivor curve with i x (standard error), i =1, 2, 3, assuming the stationary poisson process, and the similar graph in which (x, y) plots are rotated and shifted in such a way that the standard error lines and expectation lines are parallel (subroutine surviv), and output file is [out.pgSurviv] and [out.pgSurDev]. The R-module [r.pgSurviv] produces these figures in both console and postscript [pgSurviv.ps].

4. Empirical distribution of u(i) = exp{-mx(i) } where m is mean occurrence rate and x(i) is the i-th interval length between consecutive points, and lines of .95 and .99 significance bands of the two-sided Kolmogorov-Smirnov test assuming the uniform distribution are given. The related graph of ( u(i), u(i+1) ) plots are also carried out (subroutine unifrm), and output file is [out.pgInter1] and [out.pgInter2]. The R-module [r.pgInterP] produces these figures in both console and postscript [pgInterP.ps].

5. Estimated intensity m_f(t) under the palm probability. This is related to the covariance density c(t) by the relation of c(t) = \mu {m_f(t) - \mu}, where \mu is the mean intensity of the point process. The 0.95 and 0.99 error bands are provided assuming the stationary Poisson process (subroutine palmxy, palmpr), and output file is [out.pgPalm]. To make figure R-module [r.pgPalm] is used and also produces postscript output [pgPalm.ps].
.
6. Estimation of variance time curve with the 0.95 and 0.99 error lines assuming the stationary Poisson process (subroutine vtcxyp, vtcprt) , and output file is [out.pgVTC]. The R-module [r.pgVTC] produces these figures in both console and postscript [pgVTC.ps].
.
7. Subroutine structure of the FORTRAN program [pgraph.f]
        [pgraph]
           |------[input]
           |------[count1]-----[shimiz]
           |------[surviv]-----[errbr2]-----[plsinv]
           |         |---------[unifrm]-----[unitsq]
           |
           |------[palmpr]
           |------[vtc]
           |------[vtcprt]

8. The name of the analyzing data and its format is either [work.etas] or [work.res] and their format. The choice of either data is controlled by the following control file. The range of x-coordinate is limited to the non-negative valued variable.

9. Control file is [pgraph.open] in which necessary variables to read are given and explained.

10. Calculated record of program [pgraph] is stored by the name [pgraph.print] and [pgraphRes.print] corresponding to the dataset [work.etas] and [work.res], respectively, in the directory of [Calculations] for you to check of the execution program.

This program was originally designed (January 1985) and revised (December 2005) by Yosihiko Ogata, and programmed and also reprogrammed by Koichi Katsura. The R-module was designed and programmed by Jiancang Zhuang and Ogata (December 2005).


References

Cox, D.R. and Lewis, P.A.W. (1966). The Statistical Analysis of Series of Events, Methuen & Co. Ltd, London.

