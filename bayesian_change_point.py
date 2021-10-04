# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import datetime
from copy import deepcopy
import numpy as np
from scipy.special import gammaln
from shapely.geometry import Point
from openquake.hmtk.seismicity.utils import haversine
from pyrisk.etas.utils import get_region



'''
This file performs the change point calculations as described in the
in-review paper - Gupta and Baker, Estimating Spatially Varying Event 
Rates with a Change Point using Bayesian Statistics: Application to Induced Seismicity

% This script divides the state into equally spaced grid points and
% implements change point at each grid point. 
% There is a stop date defined. The stop date is for the end of training
% catalog. The end date is for the end of test catalog. 
% It calculates the probability gain 
% using the methodology specified in Werner et al, 2011, with reference to 
% a uniform rate model.
'''

class BayesianChangePoint():

    magMin = 3.5
    magMax = 8
    
    gridStep = 0.1
    radkm = 25 # Radius of circular region in km
    
    # Parameters for no change model
    k0 = 0.5
    theta0 = float("inf")
    # Parameters for single change model
    k1 = 0.5
    k2 = 0.5
    theta1 = float("inf")
    theta2 = float("inf")
    
    # Critical bayesFactor
    critBF = 0.001
    
    # Critical number of events required to detect change
    critEQ = 4 # If a region has less than critEQ events, no processing is done
    
    bayesFactors = None
    
    
    
    def __init__(self, catalog, region, startDate, stopDate, **kwargs):
        
        # override default attributes of the class
        for key, value in kwargs.items():
            if key in dir(self):
                setattr(self, key, value)

        self.catalog = catalog
        
        # Process inputs in function input format
        self.paramsCP = np.array([[self.k1, self.theta1], 
                                  [self.k2, self.theta2]]) # Parameters for change point model
        self.paramsNC = np.array([self.k0, self.theta0]) # Parameters for no change model used in Bayes Factor calculation
        
        # Create grid of points
        self.region = region
        self.latGrid, self.longGrid = self._create_point_grid(region)
        self.inpoints = self.inpolygon(self.longGrid, self.latGrid, region)
        
        self.startDate = np.datetime64(startDate)
        self.stopDate = np.datetime64(stopDate)
        
        # Decluster training catalog separately, since we are assuming that we have
        # not seen any events post the training catalog
        # catalogTraining = catalogFilter(catalog, startDate, stopDate, -2, 10) # Avoid magnitude filtering at this point
        # catalogTraining = declusterCatalog(catalogTraining) # Decluster training catalog
        self.catalogTraining = self.catalogFilter(catalog,
                                                  self.startDate, self.stopDate,
                                                  self.magMin, self.magMax) # magnitude filtering
        self.totTrainEvents = self.catalogTraining.data["magnitude"].shape[0]
        self.numTrainDays = (stopDate-startDate).days
        
        return
    
    
    
    def get_bayes_factors(self):
        if self.bayesFactors is None:
            # Run change point detection on all grid points
            bayesFactors = np.zeros(self.latGrid.shape)
            bayesFactors[:] = np.nan
            for i in range(0, self.latGrid.shape[0]):
                for j in range(0, self.latGrid.shape[1]):
                    # print(i,j)
                    # Continue only if this point is within state boundary
                    if self.inpoints[i, j]:
                        bayesFactors[i, j] = self.run_one_grid_point(i, j)    
                    else:
                        bayesFactors[i, j] = np.nan
            self.bayesFactors = bayesFactors
        return self.bayesFactors




    
    
    
    def _create_point_grid(self, region):
        maxlon = np.ceil(np.max(region["lon"])*10)/10
        minlon = np.floor(np.min(region["lon"])*10)/10
        maxlat = np.ceil(np.max(region["lat"])*10)/10
        minlat = np.floor(np.min(region["lat"])*10)/10
        lons = np.arange(minlon, maxlon+0.01, self.gridStep)
        lats = np.arange(minlat, maxlat+0.01, self.gridStep)
        ym, xm = np.meshgrid(lats, lons)
        return ym, xm



    @classmethod
    def inpolygon(cls, lons, lats, region):
        region_pol = get_region(region["lon"], region["lat"])
        grid2d = np.array([lons.flatten(), lats.flatten()]).transpose()
        flag = np.array([Point(xxx, yyy).within(region_pol) for xxx, yyy in
                         zip(grid2d[:,0], grid2d[:,1])])
        flag = flag.reshape(lons.shape)
        return flag
    
    
    
    '''
    This function filters the catalog for events between startDate and
    endDate, magnitudes between magMin and MagMax, and in a local region
    defined by localVect - [siteLat, siteLong, radkm]. localVect is optional.
    The function also returns the last date recorded in the catalog, prior to
    filtering
    catalog is an instance of the Openquake class for catalog
    startDate and endDate: datetime.datetime objects
    magMin, magMax: floats
    localVect: [lat, lon, radkm] # Circular regior
    '''
    @classmethod
    def catalogFilter(cls, catalog, startDate, endDate, magMin, magMax,
                      localVect=None):
        catalog = deepcopy(catalog)
        
        if localVect is None:
            localFlag = 0
        else:
            localFlag = 1
            siteLat = localVect[0]
            siteLong = localVect[1]
            radkm = localVect[2]
        
        if len(catalog.data["magnitude"]) < 1:
            return catalog
        
        # Determine the time range over which to count events
        catalogDates = catalog.data["datetime"]
        lastDate = catalogDates[-1]
        indx = (catalogDates >= startDate) & (catalogDates <= endDate)
        catalog.purge_catalogue(indx) # Reduce catalog to catalog of interest
        
        # Determine the magnitude range over which to count events
        catalogMags = catalog.data["magnitude"]
        indx = (catalogMags >= magMin) & (catalogMags <= magMax)
        catalog.purge_catalogue(indx)
        
        # Determine the local region over which to count events
        if localFlag:
            latVect = catalog.data["latitude"]
            lonVect = catalog.data["longitude"]
            distanceVect = haversine(lonVect, latVect, siteLong, siteLat)
            indx = distanceVect <= radkm
            catalog.purge_catalogue(indx)
        return catalog
        
    
    
    @classmethod
    def getTimeBetweenEvents(cls, catalog):
        catalogDates = catalog.data["datetime"]
        return np.diff(catalogDates)/np.timedelta64(1, "D")
    
    
    
    def run_one_grid_point(self, i, j):
        # Get earthquake catalog within circular region
        localVect = [self.latGrid[i, j], self.longGrid[i, j], self.radkm] # Circular region
        catalogTrainingLocal = self.catalogFilter(self.catalogTraining, 
                                                  self.startDate,
                                                  self.stopDate,
                                                  self.magMin,
                                                  self.magMax,
                                                  localVect)
        # Obtain timevect for inter-event times
        if catalogTrainingLocal.data["magnitude"].shape[0] != 0:
            dataVect = self.getTimeBetweenEvents(catalogTrainingLocal)
            self.save = dataVect
            startDateCP = np.min(catalogTrainingLocal.data["datetime"])
            # Change data vector to include information that no events occured between
            # startdate and startdateCP (date of first event in local calaog)
            numDaysNoEvent = (self.startDate - startDateCP)/np.timedelta64(1, "D")
            if numDaysNoEvent > 0:
                dataVect = np.append(numDaysNoEvent, dataVect)
            # Perform change-point analysis only if events in local region exceed
            # critEQ
            if catalogTrainingLocal.data["magnitude"].shape[0] >= self.critEQ:
                # Calculate bayes factor and change point for stop date data set
                bFCP = self.changePointBayesFactorGeneralized(dataVect, 
                                                              self.paramsCP,
                                                              self.paramsNC,
                                                              self.startDate)
                return bFCP
            else:
                return np.nan
        else:
            return np.nan


            #     if bFCP <= self.critBF: # Implement change point model
            #         [datesVect, probVect] = changePointProbability(dataVect, paramsCP, startDate)
            #         [maxProb, maxIndx] = max(probVect)
            #         dateCP = datesVect(maxIndx)

            #         changeDatesVect(j) = dateCP





    '''
    getDatesVect returns a dates vector containing every day in the time range
    specified in dataVect and starting on startDate to be used in changePoint
    calculations
    '''
    @classmethod
    def getDatesVect(cls, dataVect, startDate):
        timeVect = np.cumsum(dataVect)
        totDays = np.floor(timeVect[-1]) + 1    # 1 is added to account for startDate
        endDate = startDate + np.timedelta64(int(totDays), 'D')
        datesVect = np.arange(startDate+np.timedelta64(1, 'D'), endDate, dtype='datetime64[D]')
        return datesVect, totDays



    '''
    returns a number of events vector containing number of events
    that have occured till every day in the time range
    specified in dataVect and starting on startDate to be used in changePoint
    calculations
    '''
    @classmethod
    def getNumEventsVect(cls, dataVect, totDays):
        timeVect = np.cumsum(dataVect)
        days = np.arange(1, totDays)
        totEvents = len(dataVect) + 1 # 1 is added to account for first event on day 0
        numEventsVect = list()
        for t in days:
            numEventsVect.append(np.sum(timeVect <= t) + 1)
        numEventsVect[-1] = numEventsVect[-1]+1
        return np.array(numEventsVect), totEvents


    '''
    # scaleExponentials scales the vector logVect such that when each term in
    # the vector is exponentiated and summed, the sum is less than infinity i.e.
    # the sum is within the range of floating point calculations in MATLAB.
    # The unscaled value can be obtained by expVect*exp(scale) 
    '''
    @classmethod
    def scaleExponentials(cls, logVect):
        maxScale = 1e250 # The number to which the maximum of the exp(logVect - mean) is scaled.
        scaleRatio = 1e10 # Ratio by which the maxScale is reduced if overflow is not eliminated
        # Subtract mean of logs to get first scale
        scaleMean = np.mean(logVect)
        logVect = logVect - scaleMean
        # Scale the largest value such that the sum is not infinity
        maxLog = np.max(logVect)
        scale = 0
        while True:
            expVect = np.exp(logVect - scale)
            if not np.isinf(np.sum(expVect)):
                break
            else:
                scale = maxLog - np.log(maxScale)
                maxScale = maxScale/scaleRatio
        scale = scale + scaleMean
        return expVect, scale


    '''
    changePointBayesFactor calculates the bayes factor for no change against a
    change point using Raftery Akman (1986) approach.
    changePointBayesFactorGeneralized calculates the Bayes Factor of the
    change point over the time range specified in dataVect.
    dataVect is the vector with the time between each successive event. 
    paramsCM is a 2x2 matrix for parameters of the change model prior rates 
    where each row defines the shape and scale of the event rate
    before and after the change point in a gamma distribution. 
    paramsNC is 1x2 matrix for parameters of no change model. The event
    rates are assumed with gamma distribtuion priors and the change point
    is assumed with a continuous uniform distribtuion over the total
    duration of data. 
    Smaller than 1 Bayes Factor imply a shift towards change model.
    The file uses a boundary condition that a single event mid-way in the
    observation period should yield a Bayes factor of 1.
    '''
    def changePointBayesFactorGeneralized(self, dataVect, paramsCM,
                                          paramsNC, startDate):
        k1 = paramsCM[0, 0]
        theta1 = paramsCM[0, 1]
        k2 = paramsCM[1, 0]
        theta2 = paramsCM[1, 1]
        
        k0 = paramsNC[0]
        theta0 = paramsNC[1]

        if np.sum(dataVect) < 1:
            return np.nan
        # Obtain date and event vectors from data
        datesVect, totDays = self.getDatesVect(dataVect, startDate)
        numEventsVect, totEvents = self.getNumEventsVect(dataVect, totDays)

        # Calculate numerater of Bayes Factor
        if np.isinf(theta0):
            numLog = gammaln(totEvents + k0) - gammaln(k0) - (totEvents + k0)*np.log(1/theta0 + totDays)
        else:
            numLog = gammaln(totEvents + k0) - gammaln(k0) - k0*np.log(theta0) - (totEvents + k0)*np.log(1/theta0 + totDays)
        
        # Calculate denominator of Bayes Factor
        tau = (datesVect - startDate)/np.timedelta64(1, "D")
        numEvents = numEventsVect
        r1 = numEvents + k1
        r2 = totEvents - numEvents + k2
        s1 = tau + 1/theta1
        s2 = totDays - tau + 1/theta2
        denVect = -np.log(totDays) + gammaln(r1) + gammaln(r2) - r1*np.log(s1) - r2*np.log(s2)
        
        denVect, scale = self.scaleExponentials(denVect) # Exponentiate probability vector to account for very large or very small values.
        if np.isinf(theta1):
            if np.isinf(theta2):
                denLog = np.log(np.sum(denVect)) + scale - gammaln(k1) - gammaln(k2)
            else:
                denLog = np.log(np.sum(denVect)) + scale - k2*np.log(theta2) - gammaln(k1) - gammaln(k2)
        elif np.isinf(theta2):
            denLog = np.log(np.sum(denVect)) + scale - k1*np.log(theta1) - gammaln(k1) - gammaln(k2)
        else:
            denLog = np.log(np.sum(denVect)) + scale - k1*np.log(theta1) - k2*np.log(theta2) - gammaln(k1) - gammaln(k2)
            
        bayesFactor = np.exp(numLog - denLog) # This value does not include the proportionality factor based on boundary condition.
        
        # Calculate the proportioanlity factor based on description by Raftery
        # Akman
        if k1 != k2:
            print('The constant in Bayes factor calculation requires the', 
                  ' parameters k1 and k2 to be equal.', 'Currently, k1 = ', str(k1), 
                  ' and k2 = ', str(k2), 'This Bayes factor may only be used for ',
                    'comparison between datasets with equal days of observation.')
            return
        else:
            # Calculate appropriate constant factor by equating the Bayes factor
            # for the boundary condition to 1. The boundary condition is one event
            # happening at the middle of the time range.
            totEvents = 1
            if np.isinf(theta0):
                numLogBoundary = gammaln(totEvents + k0) - gammaln(k0) - (totEvents + k0)*np.log(1/theta0 + totDays)
            else:
                numLogBoundary = gammaln(totEvents + k0) - gammaln(k0) - k0*np.log(theta0) - (totEvents + k0)*np.log(1/theta0 + totDays)
            tau = np.arange(1, totDays)#'
            numEvents = np.zeros(len(tau))
            numEvents[int(np.ceil(totDays/2))-1:] = 1
            r1 = numEvents + k1
            r2 = totEvents - numEvents + k2
            s1 = tau + 1/theta1
            s2 = totDays - tau + 1/theta2
            denVectBoundary = -np.log(totDays) + gammaln(r1) + gammaln(r2) - r1*np.log(s1) - r2*np.log(s2)
            denVectBoundary, scale = self.scaleExponentials(denVectBoundary)   # Exponentiate probability vector to account for very large or very small values.
            if np.isinf(theta1):
                if np.isinf(theta2):
                    denLogBoundary = np.log(np.sum(denVectBoundary)) + scale - gammaln(k1) - gammaln(k2)
                else:
                    denLogBoundary = np.log(np.sum(denVectBoundary)) + scale - k2*np.log(theta2) - gammaln(k1) - gammaln(k2)
            elif np.isinf(theta2):
                denLogBoundary = np.log(np.sum(denVectBoundary)) + scale - k1*np.log(theta1) - gammaln(k1) - gammaln(k2)
            else:
                denLogBoundary = np.log(np.sum(denVectBoundary)) + scale - k1*np.log(theta1) - k2*np.log(theta2) - gammaln(k1) - gammaln(k2)
            bayesFactorBoundary = np.exp(numLogBoundary - denLogBoundary)
            constFactor = 1/bayesFactorBoundary
            bayesFactor = constFactor*bayesFactor
        
        return bayesFactor



