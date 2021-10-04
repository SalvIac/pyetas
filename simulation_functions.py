# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import numpy as np
from scipy.stats import uniform


''' magnitude pdf cdf and inverse cdf '''

def pdf_magnitude_trunc(m, b, mmin, mmax):
    return b*np.log(10) * 10.**(-b*(m-mmin)) / (1.-10.**(-b*(mmax-mmin)))
    
def cdf_magnitude_trunc(m, b, mmin, mmax):
    return (1.-10.**(-b*(m-mmin))) / (1.-10.**(-b*(mmax-mmin)))

def inv_cdf_magnitude_trunc(u, b, mmin, mmax):
    return mmin + (np.log10(1.-u*(1.-10.**(-b*(mmax-mmin)))))/(-b)


''' time (Omori) pdf cdf and inverse cdf '''

def pdf_time(t, c, p):
    return (p - 1)/c * np.power(1. + t/c, -p)

def cdf_time(t, c, p):
    return 1 - np.power(1+t/c, 1-p)
    
def inv_cdf_time(u, c, p):
    return c*(np.power(1-u, 1/(1-p))-1)


''' truncated time (Omori) pdf cdf and inverse cdf '''

def pdf_time_trunc(t, c, p, ta):
    if p != 1:
        pdf = (1-p)/(np.power(c+ta, 1-p) - np.power(c, 1-p)) * np.power(c+t, -p)
    else:
        pdf = np.power(np.log(c+ta)-np.log(c), -1) * np.power(c+t, -1)
    pdf[t > ta] = 0.
    return pdf
    
def cdf_time_trunc(t, c, p, ta):
    if p != 1:
        cdf = (np.power(c+t, 1-p) - np.power(c, 1-p)) / (np.power(c+ta, 1-p) - np.power(c, 1-p))
    else:
        cdf = (np.log(c+t)-np.log(c)) / (np.log(c+ta)-np.log(c))
    cdf[t > ta] = 1.
    return cdf
    
def inv_cdf_time_trunc(u, c, p, ta):
    if p != 1:
        t = np.power(u * (np.power(c+ta, 1-p) - np.power(c, 1-p)) + np.power(c, 1-p), 1/(1-p)) - c
    else:
        t = np.exp(u * (np.log(c+ta)-np.log(c)) + np.log(c)) - c
    t[t > ta] = 0.
    return t


''' space (model 5 Zhuang et al. 2011) pdf cdf and inverse cdf '''

def pdf_r_space5(r, dm, q, D, gamma):
    pdf = (q-1) / (D*np.exp(gamma*dm)*np.pi) * np.power(1+r**2/(D*np.exp(gamma*dm)), -q)
    return pdf

def pdf_space5(x, y, dm, q, D, gamma):
    r = np.sqrt(x**2+y**2)
    pdf = pdf_r_space5(r, dm, q, D, gamma)
    return pdf

def cdf_r_space5(r, dm, q, D, gamma):
    cdf = 1 - np.power(1+r**2/(D*np.exp(gamma*dm)), 1-q)
    return cdf

def cdf_space5(x, y, dm, q, D, gamma):
    r = np.sqrt(x**2+y**2)
    cdf = cdf_r_space5(r, dm, q, D, gamma)
    return cdf

def inv_cdf_space5(u_r, u_theta, dm, q, D, gamma):
    # inverse cdf of the relative locations (f)
    r = np.sqrt( D * np.exp(gamma*dm) * ((1-u_r)**(-1/(q-1))-1) )
    theta = 2*np.pi*u_theta
    return r*np.cos(theta), r*np.sin(theta)


''' truncated space (truncated model 5 Zhuang et al. 2011) pdf cdf and inverse cdf '''

# def fr(r, dm, q, D, gamma):
#     sig = D * np.exp(gamma * dm)
#     return (1 - (1 + r**2 / sig)**(1 - q)) / (2 * np.pi)

def pdf_r_space5_trunc(r, dm, q, D, gamma, r_trunc):
    trunc_factor = 1 - np.power(1+r_trunc**2/(D*np.exp(gamma*dm)), 1-q)
    pdf = (q-1) / (D*np.exp(gamma*dm)*np.pi) * np.power(1+r**2/(D*np.exp(gamma*dm)), -q) / trunc_factor
    pdf[r > r_trunc] = 0.
    return pdf

# def pdf_space5_trunc(x, y, dm, q, D, gamma, r_trunc):
#     r = np.sqrt(x**2+y**2)
#     pdf = pdf_r_space5_trunc(r, dm, q, D, gamma, r_trunc)
#     return pdf

# def cdf_r_space5_trunc(r, dm, q, D, gamma, r_trunc):
#     trunc_factor = fr(r_trunc, dm, q, D, gamma)
#     cdf = (1 - np.power(1+r**2/(D*np.exp(gamma*dm)), 1-q) ) / trunc_factor
#     # cdf[r > r_trunc] = 1.
#     return cdf

# def cdf_space5_trunc(x, y, dm, q, D, gamma, r_trunc):
#     r = np.sqrt(x**2+y**2)
#     cdf = cdf_r_space5_trunc(r, dm, q, D, gamma, r_trunc)
#     return cdf

def inv_cdf_r_space5_trunc(u_r, dm, q, D, gamma, r_trunc):
    # inverse cdf of the relative locations (f)
    r = np.sqrt( D * np.exp(gamma*dm) * ((1-u_r)**(-1/(q-1))-1) )
    if np.any(r > r_trunc):
        r[r > r_trunc] = inv_cdf_r_space5_trunc(uniform.rvs(size=np.sum(r > r_trunc)), 
                                                dm, q, D, gamma, r_trunc)
    return r

def inv_cdf_space5_trunc(u_r, u_theta, dm, q, D, gamma, r_trunc):
    # inverse cdf of the relative locations (f)
    r = inv_cdf_r_space5_trunc(u_r, dm, q, D, gamma, r_trunc)
    theta = 2*np.pi*u_theta
    return r*np.cos(theta), r*np.sin(theta)


''' space (model 3 Zhuang et al. 2011) pdf cdf and inverse cdf '''

def inv_cdf_space3(u_r, u_theta, dm, q, D):
    # inverse cdf of the relative locations (f) - model3
    r = np.sqrt( D * ((1-u_r)**(-1/(q-1))-1) )
    theta = 2*np.pi*u_theta
    return r*np.cos(theta), r*np.sin(theta)




''' completeness vs time after event in catalog '''

def compl_vs_time_hs06(mag, time):
    # estimated incompleteness function for California (Helmstetter et al., 2006)
    # where time is the time (in days) after an earthquake with magnitude mag
    mc = mag-4.5-0.75*np.log10(time)
    return mc

    
def inv_compl_vs_time_hs06(mag, mc):
    # inverse incompleteness function for California (Helmstetter et al., 2006)
    # where time is the time (in days) after an earthquake with magnitude mag
    time = 10**(-(mc + 4.5 - mag)/0.75)
    return time


def compl_vs_time_p16(mag, time):
    # global estimated incompleteness function (Page et al. 2016)
    # where time is the time (in days) after an earthquake with magnitude mag
    mc = mag/2-0.25-np.log10(time)
    return mc

    
def inv_compl_vs_time_p16(mag, mc):
    # inverse incompleteness function for California (Helmstetter and Shaw, 2006)
    # where time is the time (in days) after an earthquake with magnitude mag
    time = 10**(-(mc + 0.25 - mag/2))
    return time    


def compl_vs_time_general(mag, time, *args):
    # global estimated incompleteness function (Page et al. 2016)
    # where time is the time (in days) after an earthquake with magnitude mag
    mc = args[0]*mag-args[1]-args[2]*np.log10(time)
    return mc


def inv_compl_vs_time_general(mag, mc, *args):
    # global estimated incompleteness function (Page et al. 2016)
    # where time is the time (in days) after an earthquake with magnitude mag
    time = 10**(-(mc + args[1] - args[0]*mag)/args[2])
    return time




#%%
    

if __name__ == "__main__":
    
    
    from scipy.stats import uniform
    from scipy.interpolate import interp1d
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
    


    # ''' incompleteness function '''

    # ttt = np.arange(0.0001, 10, 0.0001)
    # mc1 = compl_vs_time_hs06(7, ttt)
    # mc2 = compl_vs_time_p16(7, ttt)
    # fig = plt.figure()
    # plt.plot(ttt, mc1)
    # plt.plot(ttt, mc2)
    # plt.show()


    
    # ''' magnitude distributions '''
    
    # b = 1.
    # mmin = 3.5
    # mmax = 8.
    # dm = 0.1
    # # t = np.arange(dt/2, tmax, dt) # center of bins
    # m = np.arange(mmin, mmax+dm, dm) # bins
    # # u = uniform.rvs(size=100000)
    # u = np.arange(0,1,1e-6)
    # # for discrete integral



    # ###################### untruncated GR distribution #####################

    # pdfm = pdf_magnitude_trunc(m, b, mmin, mmax)
    # cdfm = cdf_magnitude_trunc(m, b, mmin, mmax)
    # invcdfm = inv_cdf_magnitude_trunc(u, b, mmin, mmax)


    # # general plot
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # ax[0].plot(m, pdfm, 'r-', lw=2, alpha=0.6,
    #            label='truncated GR pdf b='+str(b)+' mmax='+str(mmax))
    # # histogram simulations
    # hist = np.histogram(invcdfm, bins=m)
    # dth = np.min(np.diff(hist[1]))
    # ax[0].bar(hist[1][:-1]+dth/2, hist[0]/u.shape[0]/dth,
    #           width=dth, alpha=0.2)
    # ax[0].set_ylabel('f(m)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # ax[1].plot(m, cdfm, 'r-', lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('m')
    # ax[1].set_ylabel('F(m)')
    # plt.show()
    
    
    # # plot for effect of b
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for bb in [0.8, 1, 1.2]:
    #     ax[0].plot(m, pdf_magnitude_trunc(m, bb, mmin, mmax), lw=2, alpha=0.6,
    #                label='truncated GR pdf b='+str(bb)+' mmax='+str(mmax))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for bb in [0.8, 1, 1.2]:
    #     ax[1].plot(m, cdf_magnitude_trunc(m, bb, mmin, mmax), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()

    

    # # plot for effect of mmax
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for mmmax in [4, 6, 8]:
    #     ax[0].plot(m, pdf_magnitude_trunc(m, 1., mmin, mmmax), lw=2, alpha=0.6,
    #                label='truncated GR pdf b='+str(bb)+' mmin='+str(mmin)+' mmax='+str(mmmax))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf
    # for mmmax in [4, 6, 8]:
    #     ax[1].plot(m, cdf_magnitude_trunc(m, 1., mmin, mmmax), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    
    # # plot for effect of mmin
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for mmmin in [3.5, 4., 4.5]:
    #     ax[0].plot(m, pdf_magnitude_trunc(m, 1., mmmin, mmax), lw=2, alpha=0.6,
    #                label='truncated GR pdf b='+str(bb)+' mmin='+str(mmmin)+' mmax='+str(mmax))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for mmmin in [3.5, 4., 4.5]:
    #     ax[1].plot(m, cdf_magnitude_trunc(m, 1., mmmin, mmax), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    
    # ''' time distributions '''
    
    # c = 0.005
    # p = 1.1
    # tmax = 10*365
    # dt = 1
    # # t = np.arange(dt/2, tmax, dt) # center of bins
    # t = np.arange(0, tmax, dt) # bins
    # # u = uniform.rvs(size=100000)
    # u = np.arange(0,1,1e-6)
    # # for discrete integral
    # div = 10000
    # t2 = np.arange(dt/2/div, tmax, dt/div) # center of bins
    # ta = 10*365 # 10yrs

    # # from scipy.integrate import trapz, simps
    # # trapz(g, ttt)
    # # simps(g, ttt)
    

    # ###################### untruncated time distribution #####################
    
    # pdft = pdf_time(t, c, p)
    # cdft = cdf_time(t, c, p)
    # invcdft = inv_cdf_time(u, c, p)

    # # check integral pdf equal to cdf (discrete integral)
    # t2 = np.arange(dt/2/div, tmax, dt/div) # center of bins
    # pdft2 = pdf_time(t2, c, p)
    # cdft2 = cdf_time(t2, c, p)
    # np.testing.assert_almost_equal(cdft2, np.cumsum(pdft2*dt/div), decimal=3)
    # # fig, ax = plt.subplots(1, 1)
    # # plt.plot(t2, cdft2)
    # # plt.plot(t2, np.cumsum(pdft2*dt/div))
    # # plt.show()

    
    # # general plot
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # ax[0].plot(t, pdft, 'r-', lw=2, alpha=0.6,
    #            label='untruncated time pdf c='+str(c)+' p='+str(p))
    # # histogram simulations
    # hist = np.histogram(invcdft, bins=np.arange(0, tmax, 10*dt))
    # dth = np.min(np.diff(hist[1]))
    # ax[0].bar(hist[1][:-1]+dth/2, hist[0]/u.shape[0]/dth,
    #           width=dth, alpha=0.2)
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # ax[1].plot(t, cdft, 'r-', lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    # # plot for effect of p
    # cc = c
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for pp in [1.01, 1.05, 1.1, 1.2, 1.3]:
    #     ax[0].plot(t, pdf_time(t, cc, pp), lw=2, alpha=0.6,
    #                label='untruncated time pdf c='+str(cc)+' p='+str(pp))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for pp in [1.01, 1.05, 1.1, 1.2, 1.3]:
    #     ax[1].plot(t, cdf_time(t, cc, pp), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    # # plot for effect of c
    # pp = p
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for cc in [0.001, 0.005, 0.01, 0.05, 0.1]:
    #     ax[0].plot(t, pdf_time(t, cc, pp), lw=2, alpha=0.6,
    #                label='untruncated time pdf c='+str(cc)+' p='+str(pp))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for cc in [1.01, 1.05, 1.1, 1.2, 1.3]:
    #     ax[1].plot(t, cdf_time(t, cc, pp), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    


    # ####################### truncated time distribution ######################

    # pdft_trunc = pdf_time_trunc(t, c, p, ta)
    # cdft_trunc = cdf_time_trunc(t, c, p, ta)
    # invcdft_trunc = inv_cdf_time_trunc(u, c, p, ta)
    
    # # check integral pdf equal to cdf (discrete integral)
    # pdft2_trunc = pdf_time_trunc(t2, c, p, ta)
    # cdft2_trunc = cdf_time_trunc(t2, c, p, ta)
    # np.testing.assert_almost_equal(cdft2_trunc, np.cumsum(pdft2_trunc*dt/div), decimal=3)
    # # fig, ax = plt.subplots(1, 1)
    # # plt.plot(t2, cdft2_trunc)
    # # plt.plot(t2, np.cumsum(pdft2_trunc*dt/div))
    # # plt.show()


    # # general plot
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # ax[0].plot(t, pdft_trunc, 'b-', lw=2, alpha=0.6,
    #            label='truncated time pdf ta=10yrs c='+str(c)+' p='+str(p))
    # # histogram simulations
    # hist = np.histogram(invcdft_trunc, bins=np.arange(0, tmax, 10*dt))
    # dth = np.min(np.diff(hist[1]))
    # ax[0].bar(hist[1][:-1]+dth/2, hist[0]/u.shape[0]/dth,
    #           width=dth, alpha=0.2)
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # ax[1].plot(t, cdft_trunc, 'b-', lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    # # plot for effect of p
    # cc = c
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for pp in [0.7, 0.9, 1, 1.1, 1.3]:
    #     ax[0].plot(t, pdf_time_trunc(t, cc, pp, ta), lw=2, alpha=0.6,
    #                label='truncated time pdf ta=10yrs c='+str(cc)+' p='+str(pp))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for pp in [1.01, 1.05, 1.1, 1.2, 1.3]:
    #     ax[1].plot(t, cdf_time_trunc(t, cc, pp, ta), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()
    
    
    # # plot for effect of c
    # pp = p
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # for cc in [0.001, 0.005, 0.01, 0.05, 0.1]:
    #     ax[0].plot(t, pdf_time_trunc(t, cc, pp, ta), lw=2, alpha=0.6,
    #                label='truncated time pdf ta=10yrs c='+str(cc)+' p='+str(pp))
    # ax[0].set_ylabel('f(t)')
    # ax[0].set_yscale("log")
    # ax[0].legend()
    # # cdf    
    # for cc in [1.01, 1.05, 1.1, 1.2, 1.3]:
    #     ax[1].plot(t, cdf_time_trunc(t, cc, pp, ta), lw=2, alpha=0.6)
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()



    # ################# comparison time distributions plots ####################
    
    # fig, ax = plt.subplots(2, 1, sharex=True)
    # # pdf
    # ax[0].plot(t, pdft, 'r-', lw=2, alpha=0.6, label='untruncated time pdf')
    # ax[0].plot(t, pdft_trunc, 'b-', lw=2, alpha=0.6, label='truncated time pdf ta=10yrs')
    # ax[0].set_yscale("log")
    # ax[0].set_ylabel('f(t)')
    # ax[0].legend()
    # # cdf    
    # ax[1].plot(t, cdft, 'r-', lw=2, alpha=0.6, label='untruncated time pdf')
    # ax[1].plot(t, cdft_trunc, 'b-', lw=2, alpha=0.6, label='truncated time pdf ta=10yrs')
    # ax[1].set_ylim([0,1])
    # ax[1].set_xlabel('t (days)')
    # ax[1].set_ylabel('F(t)')
    # plt.show()




    # # import matplotlib.pyplot as plt
    # # fig = plt.figure()
    # # ax = fig.gca()
    # # ss = ax.scatter(t_offspr, m_offspr)
    # # ax.set_xlim(0,3000)
    # # plt.show()    





    ''' space distributions '''
    
    q = 1.5
    D = 0.005
    gamma = 0.8
    dm = 6.5-3.5
    
    dr = 0.02
    xymax = 2 # this is in degree


    ######################### 2D spatial distribution ########################

    x, y = np.mgrid[-xymax:xymax+dr:dr, -xymax:xymax+dr:dr]
    pdfs = pdf_space5(x, y, dm, q, D, gamma)
    cdfs = cdf_space5(x, y, dm, q, D, gamma)

    # general plot
    y0 = int(x.shape[0]/2)
    
    # fig = plt.figure()
    
    # ax = fig.add_subplot(221, projection='3d')
    # ax.scatter(x, y, pdfs, s=0.1, alpha=0.5)
    # ax.plot(x[:,y0], y[:,y0], pdfs[:,y0], color='r')
    # ax.set_xlabel('x (deg)')
    # ax.set_ylabel('y (deg')

    # ax = fig.add_subplot(222)
    # ax.plot(x[:,y0], pdfs[:,y0], color='r')
    # ax.set_xlim([0,xymax])
    # ax.set_xlabel('x/r (deg)')
    # ax.set_ylabel('f(r)')
    
    # ax = fig.add_subplot(223, projection='3d')
    # ax.scatter(x, y, cdfs, s=0.1, alpha=0.5)
    # ax.plot(x[:,y0], y[:,y0], cdfs[:,y0], color='r')
    # ax.set_zlim([0,1])
    # ax.set_xlabel('x (deg)')
    # ax.set_ylabel('y (deg')

    # ax = fig.add_subplot(224)
    # ax.plot(x[:,y0], cdfs[:,y0], color='r')
    # ax.set_xlim([0,xymax])
    # ax.set_ylim([0,1])
    # ax.set_xlabel('x/r (deg)')
    # ax.set_ylabel('F(r)')

    # plt.show()
    
    
    # # check one value of the cdf (pdf integral)
    # x2, y2 = np.mgrid[-xymax+dr/2:xymax:dr, -xymax+dr/2:xymax+dr:dr]
    # pdfs2 = pdf_space5(x2, y2, dm, q, D, gamma)
    # cdfs2 = cdf_space5(x2, y2, dm, q, D, gamma)
    # # radial integral pdf   
    # r = np.sqrt(x2**2+y2**2)
    # rlim = np.arange(0.1,2,0.1) # to check the cdf value
    # check = list()
    # for rl in rlim:
    #     ind = r < rl
    #     check.append( np.sum(pdfs2[ind]*dr*dr) ) # integral
    # # cdf values
    # f = interp1d(r[int(r.shape[0]/2):, int(r.shape[1]/2)],
    #              cdfs2[int(r.shape[0]/2):, int(r.shape[1]/2)])
    # check_cdf = f(rlim)
    # np.testing.assert_almost_equal(check_cdf, check, decimal=2)
    # # fig, ax = plt.subplots(1, 1)
    # # plt.plot(rlim, check_cdf)
    # # plt.plot(rlim, check)
    # # plt.show()
    
    
    # # simulations
    # u_r = np.arange(0+0.5e-2,1,1e-2)
    # u_theta = np.arange(0+1e-2,1,1e-2)
    # u_r, u_theta = np.meshgrid(u_r, u_theta)
    # simx, simy = inv_cdf_space5(u_r, u_theta, dm, q, D, gamma)
    # # u = uniform.rvs(size=10000)
    # # u = np.linspace(0,1,10000)
    # # r = np.sqrt( D*((1-u)**(-1/(q-1))-1)/np.exp(-gamma*dm) ) # inverse cdf of gaussian kernel function
    # # theta = uniform.rvs(loc=0., scale=np.pi, size=10000) # ni[par])
    # # simx = r*np.cos(theta)
    # # simy = r*np.sin(theta)    

    # fig = plt.figure()
    # ax = fig.add_subplot(111, projection='3d')
    # hist, xedges, yedges = np.histogram2d(simx.flatten(), simy.flatten(),
    #                                       bins=202, range=[[-xymax, xymax], [-xymax, xymax]])
    # xpos, ypos = np.meshgrid(xedges[:-1] + 0.01, yedges[:-1] + 0.01, indexing="ij")
    # xpos = xpos.ravel()
    # ypos = ypos.ravel()
    # zpos = 0
    # dx = dy = 0.01 * np.ones_like(zpos)
    # dz = hist.ravel()/(np.sum(hist)*np.min(np.diff(xedges))*np.min(np.diff(yedges)))
    # ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
    # # ax.scatter(x, y, pdfs, s=0.1)
    # plt.show()
    
    
    
    # # check 2d
    # u_r = np.arange(0,1,1e-3)
    # u_theta = 0.
    # simx, simy = inv_cdf_space5(u_r, u_theta, dm, q, D, gamma)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, gamma))
    # ax.plot(simx.flatten(), u_r.flatten(), "--")
    # ax.set_xlim([0,2])
    # ax.set_ylim([0,1])
    # plt.show()
    
    
    
    
    # # check 2d - spatial pdf (effect of q)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(0,2,1e-3), pdf_space5(0, np.arange(0,2,1e-3), dm, 2.42, 0.00191, 0.771),
    #         label="q=1.2")
    
    
    
    # ax.set_xlim([0,2])
    # # ax.set_ylim([0,1])
    # ax.legend()
    # plt.show()


    # cdf_r_space5(1., dm, 2.431, 0.001872, 0.7775)
    # inv_cdf_space5(1., 0., dm, 2.431, 0.001872, 0.7775)


    



























    
    # # check 2d - spatial pdf (effect of q)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, 1.2, D, gamma),
    #         label="q=1.2")
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, gamma),
    #         label="q="+str(q))
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, 2., D, gamma),
    #         label="q=2.0")
    # ax.set_xlim([0,2])
    # ax.legend()
    # plt.show()  
    
    
    # # check 2d - spatial pdf (effect of D)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D*10, gamma),
    #         label="D="+str(D*10))
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, gamma),
    #         label="D="+str(D))
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D/10, gamma),
    #         label="D="+str(D/10))
    # ax.set_xlim([0,2])
    # ax.legend()
    # plt.show()  
    
    
    # # check 2d - spatial pdf (effect of gamma)
    # fig = plt.figure()
    # ax = fig.add_subplot(111)
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, 0.4),
    #         label="gamma=0.8")
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, gamma),
    #         label="gamma="+str(gamma))
    # ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, q, D, 1.5),
    #         label="gamma=1.5")
    # ax.set_xlim([0,2])
    # ax.legend()
    # plt.show()
    
    import sys
    sys.path.append('C:\\Users\\Salvatore\\Dropbox\\SalvIac')
    from myutils.latex_plot_mode import latex_plot_mode
    latex_plot_mode()
    
    # check 2d - spatial pdf (effect of gamma)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, 1.5, 1.8e-5, 0.79),
            label=r"$q=1.5$  $D=1.8e-4$  $\gamma=0.79$")
    ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, 1.84, 2.7e-4, 0.8),
            label=r"$q=1.84$  $D=2.7e-4$  $\gamma=0.80$")
    ax.plot(np.arange(0,2,1e-3), cdf_r_space5(np.arange(0,2,1e-3), dm, 2.0670190147517316, 0.00037982731849683953, 1.0327436857148076),
            label=r"$q=1.84$  $D=2.7e-4$  $\gamma=0.80$")
    ax.set_xlim([0,1])
    ax.set_xlabel("r (deg)")
    ax.set_ylabel(r"$F(r)$")
    ax.legend(loc='lower right')
    plt.show()    
    stop
    
    
    
    # simulations
    u_r = np.random.rand(100000)
    u_theta = np.random.rand(100000)
    simx, simy = inv_cdf_space5(u_r, u_theta, dm, q, D, gamma)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    hist, xedges, yedges = np.histogram2d(simx.flatten(), simy.flatten(),
                                          bins=50, range=[[-xymax, xymax], [-xymax, xymax]])
    xpos, ypos = np.meshgrid(xedges[:-1] + 0.02, yedges[:-1] + 0.02, indexing="ij")
    xpos = xpos.ravel()
    ypos = ypos.ravel()
    zpos = 0
    dx = dy = 0.01 * np.ones_like(zpos)
    dz = hist.ravel()/(np.sum(hist)*np.min(np.diff(xedges))*np.min(np.diff(yedges)))
    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, zsort='average')
    ax.scatter(x, y, pdfs, s=0.5)
    plt.show()
    
    
    
    
    # #################### truncated 2D spatial distribution ####################
    
    rr = np.arange(0,2,1e-3)
    
    # check 2d - spatial pdf (effect of D)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(rr, pdf_r_space5(rr, dm, q, D, gamma),
            label="untruc")
    ax.plot(rr, pdf_r_space5_trunc(rr, dm, q, D, gamma, 1.),
            label="trunc")
    ax.set_xlim([0,2])
    ax.legend()
    plt.show()  










    
    
    
    