# -*- coding: utf-8 -*-
# pyetas
# Copyright (C) 2021-2022 Salvatore Iacoletti
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
# 
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
"""


import pandas as pd
import numpy as np
from scipy.stats import uniform
from shapely.geometry import Point, Polygon, box, polygon
from openquake.hmtk.seismicity.catalogue import Catalogue




class CatalogueEtas(Catalogue):
    
    '''
    creates an object of class "catalog" representing an earthquake catalog
    dataset. An earthquake catalog is a chronologically ordered list of time,
    epicenter and magnitude of all recorded earthquakes in geographical 
    region during a specific time period.
    
    data data.frame date, time, latitude, longitude and magnitude of earthquakes
    time.begin The beginning of time span of the catalog. A character string
    or an object that can be converted to date-time (calendar dates plus time
    to the nearest second). The default NULL sets it to the date-time of the
    first event.
    study.start The start of the study period. A character string or an object
    that can be converted to date-time. If not specified (NULL), then time.begin
    is used
    study.end The end of the study period. A character string or an object that
    can be convertedto date-time. The default NULL sets it to the date-time of
    the last event
    study.length A single numeric value specifying the length of the study
    period in decimal days. Incompatible with study.end
    lat.range The latitude range of a rectangular study region. A numeric
    vector of size 2 giving (latmin, latmax). By default (NULL) the range of 
    the latitudes of events is used
    long.range The longitude range of a rectangular study region. A numeric
    vector of size 2 giving (longmin, longmax). By default (NULL) the range of
    the longitudes of events is used
    region.poly Polygonal boundary of a non-rectangular study region. A list
    with components lat and long of equal length specifying the coordinates of
    the vertices of a polygonal study region. The vertices must be listed 
    in anticlockwise order
    mag.threshold The magnitude threshold of the catalog. A positive numeric 
    value. The default (NULL) sets it to the minimum magnitude of all events.
    flatmapLogical flag indicating whether to transform the spherical 
    coordinates (long,lat) on the earth surface to flat map (planar) 
    coordinates (x,y) in order to approximate  the great-circle distance on the
    sphere by the corresponding Euclideandistance on the flat map.
    dist.unitA character string specifying the unit of geographical coordinates and spatialdistances between events. Options are"degree"(the default case) and"km".tzA character string specifying the time zone to be used for the date-time conver-sion inas.POSIXlt.  The default"GMT"is the UTC (Universal Time, Coordi-nated).DetailsThedatais required to have at least 5 columns with namesdate,time,lat,longandmagcontain-ing, respectively, the date, time, latitude, longitude and magnitude of each event in the catalog.The geographical study region can be rectangular or polygonal:•rectangular study regioncan be specified bylat.rangeandlong.rangewhich must benumeric vectors of length 2.•polygonal study regioncan be specified byregion.polywhich contains coordinates of thevertices of the polygon. It must be either alistwith componentslatandlongof equal lengthor adata.framewith columnslatandlong. The vertices must be listed inanticlockwiseorderand no vertex should be repeated (i.e. do not repeat the first vertex).The functioninside.owinin thespatstatis used to indicate whether events lie inside the studyregion.  Only events inside the study region and the study period (study.start,study.end) areconsidered astargetevents. Other events are assumed to becomplementaryevents.If the events indataare not chronologically sorted, then a warning will be produced and the eventswill be sorted in ascending order with respect to time of occurrence.Ifflatmap=TRUE, longitude-latitude coordinates convert to flat map coordinates:•  ifdist.unit="degree", then the Equirectangular projectionx= cos(cnt.lat/180π)(long−cnt.long
    
    4catalogandy=lat−cnt.latis used to obtain the flat map coordinates(x,y)in degrees, wherecnt.latandcnt.longare, respectively, the latitude and longitude of the centroid of the geographicalregion.•  ifdist.unit="km", then the projectionx= 111.32 cos(lat/180π)longandy= 110.547latis used wherexandyare in (approximate) kilometers.ValueAn object of class"catalog"containing an earthquake catalog dataset.Author(s)Abdollah Jalilian<jalilian@razi.ac.ir>ReferencesZhuang J (2012).   Long-term Earthquake Forecasts Based on the Epidemic-type Aftershock Se-quence (ETAS) Model for Short-term Clustering.Research in Geophysics,2(1), 52–57. doi:10.4081/rg.2012.e8.See Alsoetas.Examplessummary(iran.quakes)# creating a catalog with rectangular study regioniran.cat = catalog(iran.quakes, time.begin="1973/01/01",study.start="1985/01/01", study.end="2016/01/01",lat.range=c(25, 42), long.range=c(42, 63),mag.threshold=4.5)print(iran.cat)## Not run:plot(iran.cat)## End(Not run)# equivalently, specifying the length of the study periodiran.cat2 = catalog(iran.quakes, time.begin="1973/01/01",study.start="1985/01/01", study.length=11322,lat.range=c(25, 42), long.range=c(42, 63),mag.threshold=4.5)print(iran.cat2)# specifying a polygonal geographical region
    
    '''
    
    def __init__(self, data, time_begin=None, study_start=None, study_end=None,
                 study_length=None, lat_range=None, long_range=None,
                 region_poly=None, mag_threshold=None, flatmap=True,
                 dist_unit="degree", roundoff=True, tz="GMT", map_faults=None):

        # accounting for round-off error in coordinates of epicenters
        if roundoff:
            # two random numbers uniformly distributed in (-0.005, 0.005)
            print("roundoff activated")
            np.random.seed(seed=42)
            # ind = data[['longitude', 'latitude']].duplicated()
            # size = np.sum(ind)
            size = data['longitude'].shape[0]
            loc = -0.005
            scale = 0.005-loc
            r = uniform.rvs(loc=loc, scale=scale, size=(size, 2))
            # data.loc[ind, 'longitude'] += r[:,0]
            # data.loc[ind, 'latitude'] += r[:,1]
            data['longitude'] += r[:,0]
            data['latitude'] += r[:,1]
        
        super().__init__()
        self.data = data
        self.update_end_year()
        self.update_start_year()
        self.number_earthquakes = self.get_number_events()
        
        data = pd.DataFrame(data)
        
        # dnames = tolower(names(data))
        # names(data) = dnames
        # vnames = c("date", "time", "long", "lat", "mag")
        # if (!all(vnames %in% dnames))
        # stop(paste("argument", sQuote("data"),
        #             "must be a data frame with column names ",
        #              toString(sQuote(vnames))))
        # if (any(is.na(data[, vnames])))
        #     stop(paste(sQuote(vnames), "must not contain NA values"))
        # if (!is.numeric(data$lat) || !is.numeric(data$long) || !is.numeric(data$mag))
        #     stop("lat, long and mag columns must be numeric vectors")


        # extract spatial coordinates and magnitude
        xx = data['longitude']  # longitude of epicenter: coordinates
        yy = data['latitude']   # latitude of epicenter : coordinates
        mm = data['magnitude']  # magnitude
        

        # extract date and time of events
        dt = pd.to_datetime(data[['year', 'month', 'day', 'hour', 'minute', 'second']])
        if dt.duplicated().any():
            raise Exception("no more than one event can occur simultaneously! check events: \n" + \
                            str(data[dt.duplicated()]))

        if not dt.is_monotonic:
            print("events were not chronologically sorted: they have been sorted in ascending order")
            data['datetime'] = dt
            data.sort_values(by=['datetime'], inplace=True)
            dt.sort_values(inplace=True)

        if time_begin is None:
            time_begin = dt.min()
        else:
            time_begin = time_begin # datetime.datetime(time_begin, tz=tz) #todo
            if (dt < time_begin).all():
              raise Exception("change time.begin: no event has occurred after: " + \
                              str(time_begin))
        
        if study_start is None:
            study_start = time_begin
        else:
            # study_start = datetime.datetime(study_start, tz=tz) #todo
            if study_start < time_begin:
             raise Exception("study.start "+str(study_start)+\
                             " cannot be set before time.begin = "+\
                             str(time_begin))

        if study_length is not None:
            if study_end is not None:
                raise Exception("either study.end or study.length needs to be specified, not both")
            else:
                if (not isinstance(study_length,float)) and \
                   (not isinstance(study_length,int)):
                    raise Exception("study.length must be single numeric value: in deciaml days")
            study_end = study_start + study_length * 24 * 60 * 60
          
        if study_end is None:
            study_end = dt.max()
        else:
            # study_end = datetime.datetime(study_end, tz=tz) #todo
            if study_end < study_start:
                raise Exception("study.end"+ str(study_end)+\
                                "can not be set before study.start"+\
                                str(study_start))
        tt = self.date2day(dt, time_begin, tz=tz)
        
        # spatial region
        if lat_range is None:
            dif = yy.max() - yy.min()
            lat_range = [yy.min()-0.01*dif, yy.max()+0.01*dif]
        elif len(lat_range) != 2 or lat_range[1] <= lat_range[0]:
            raise Exception("lat.range must be a vector of length 2 giving (lat.min, lat.max)")
        if long_range is None:
            dif = xx.max() - xx.min()
            long_range = [xx.min()-0.01*dif, xx.max()+0.01*dif]
        elif len(long_range) != 2 or long_range[1] <= long_range[0]:
            raise Exception("long.range must be a vector of length 2 giving (long.min, long.max)")
        
        if region_poly is None:
            lon = [long_range[0],long_range[1],long_range[1],long_range[0]]
            lat = [lat_range[0],lat_range[0],lat_range[1],lat_range[1]]
            coords  = [(x, y) for x,y in zip(lon, lat)]
            region_poly = dict(lon=lon, lat=lat)
            # region_win = dict(xrange=long_range, yrange=lat_range) # spatstat::owin
            region_win = box(long_range[0], lat_range[0],
                             long_range[1], lat_range[1])
        else:
            if not isinstance(region_poly, Polygon):
                if isinstance(region_poly, dict):
                    if all([k in region_poly.keys() for k in ["lat", "lon"]]) and \
                       len(region_poly['lat']) == len(region_poly['lon']):
                        coords  = [(x, y) for x,y in zip(region_poly['lon'], region_poly['lat'])]
                        # region_poly = dict(lon=region_poly['lon'], lat=region_poly['lat'])
                else:
                    raise Exception("region.poly must be a dict with components lat and lon")
            region_win = Polygon(coords)
            # adjust counter-clockwise
            region_win = polygon.orient(region_win)
            region_poly = dict(lon=region_win.exterior.xy[0][:-1],
                               lat=region_win.exterior.xy[1][:-1])
            region_area = region_win.area
            if region_area < 0:
                raise Exception("Area of polygon is negative - maybe traversed in wrong direction?")
          
        # magnitude threshold
        if mag_threshold is None:
            mag_threshold = mm.min()
        elif np.isnan(mag_threshold):
            raise Exception("mag.threshold must be a single numeric value")
        
        # project long-lat coordinates to flat map coordinates
        longlat_coord = pd.DataFrame(zip(xx,yy), columns=['long','lat'])
        if flatmap:
            proj = self.longlat2xy(long=xx, lat=yy, region_poly=region_poly,
                                   dist_unit=dist_unit)
            xx = proj['x']
            yy = proj['y']
            region_win = proj['region_win']

        ok = (dt <= study_end) & (dt >= time_begin) & (mm >= mag_threshold)
        xx = xx.loc[ok]
        yy = yy.loc[ok]
        tt = tt.loc[ok]
        mm = mm.loc[ok] - mag_threshold
        
        flag = np.array([int(Point(xxx, yyy).within(region_win)) for xxx, yyy in zip(xx,yy)])
        flag[dt[ok] < study_start] = -2

        revents = pd.DataFrame(dict(tt=tt, xx=xx, yy=yy, mm=mm, flag=flag,
                                    bkgd=0., prob=1., lambd=0.))
        longlat_coord = longlat_coord[ok]
        longlat_coord['flag'] = flag
        longlat_coord['dt'] = dt[ok]
        
        
        # fault map
        temp = np.array([None]*revents.shape[0])
        if map_faults is not None:
            for pubid in map_faults:
                ind = np.where(self.data['publicid'] == str(pubid))[0]
                if ind.shape[0] != 0:
                    proj = self.longlat2xy(long=map_faults[pubid].lons,
                                            lat=map_faults[pubid].lats,
                                            region_poly=region_poly,
                                            dist_unit=dist_unit)
                    temp[ind[0]] = {'xx': proj['x'],
                                    'yy': proj['y'],
                                    'depth': map_faults[pubid].depths}        
        revents['fault'] = temp

        
        
        # # fault map
        # temp = np.array([None]*revents.shape[0])
        # if map_faults is not None:
        #     for pubid in map_faults:
        #         print(pubid)
        #         ind = np.where(self.data['publicid'] == str(pubid))[0]
        #         if ind.shape[0] != 0:
        #             proj = self.longlat2xy(long=map_faults[pubid].lons,
        #                                     lat=map_faults[pubid].lats,
        #                                     region_poly=region_poly,
        #                                     dist_unit=dist_unit)
        #             temp[ind[0]] = {'xx': proj['x'],
        #                             'yy': proj['y'],
        #                             'depth': map_faults[pubid].depths}        
        # revents['fault'] = temp


        X = pd.DataFrame(dict(t=tt, x=xx, y=yy, m=mm)) # coord.type=c("t", "s", "s", "m"))
        
        px = list(region_win.exterior.coords.xy[0]) # ) [region_win.bounds[0], region_win.bounds[2],
              # region_win.bounds[2], region_win.bounds[0]]
        py = list(region_win.exterior.coords.xy[1]) # [region_win.bounds[1], region_win.bounds[1],
              # region_win.bounds[3], region_win.bounds[3]]
        
        # repeat the first vertex
        # px.append(px[0])
        # py.append(py[0])
        rpoly = dict(px=px, py=py)

        rtperiod = dict(study_start = self.date2day(study_start, time_begin, tz=tz),
                        study_end = self.date2day(study_end, time_begin, tz=tz))
        

        
        self.revents = revents
        self.map_faults = map_faults
        self.rpoly = rpoly
        self.rtperiod = rtperiod
        self.X = X
        self.region_poly = region_poly
        self.region_win = region_win
        self.time_begin = time_begin
        self.study_start = study_start
        self.study_end = study_end
        self.study_length=study_length
        self.mag_threshold = mag_threshold
        self.longlat_coord = longlat_coord
        self.dist_unit = dist_unit
        return


    @classmethod
    # number of decimal places
    def decimalplaces(cls, x):
        f = [str(xx) for xx in x]
        out = [ff[::-1].find('.') for ff in f]
        return out


    # @classmethod
    # def roundoffErr(x):
    #     x + [10**dx for dx in decimalplaces(x)] # runif(length(x), -1/2, 1/2) * 10**(-decimalplaces(x))

    
    @classmethod
    def date2day(cls, dates, start=None, tz=""):
        # dates is a datetime vector
        if start is None:
            start = dates.min()
        out = (dates-start) / pd.to_timedelta(1, unit='D')
        return out


    @classmethod
    def longlat2xy(cls, long, lat, region_poly, dist_unit="degree"):
        # Equirectangular projection
        if len(long) != len(lat):
            raise Exception("long and lat must be of equal length.")
        region_poly_lon = np.array(region_poly['lon'])
        region_poly_lat = np.array(region_poly['lat'])#.exterior.coords.xy[1])
        if dist_unit == "degree":
            coords = [(x, y) for x,y in zip(region_poly_lon, region_poly_lat)]
            region_win = Polygon(coords)
            region_cnt = region_win.centroid.xy
            x = np.cos(region_cnt[1][0] * np.pi / 180) * (long - region_cnt[0][0])
            y = lat - region_cnt[1][0]
            px = np.cos(region_cnt[1][0]/180 * np.pi) * (region_poly_lon - region_cnt[0][0])
            py = region_poly_lat - region_cnt[1][0]
        elif dist_unit == "km":
            x = 111.320 * np.cos(lat/180 * np.pi) * long
            y = 110.574 * lat
            px = 111.320 * np.cos(region_poly_lat/180 * np.pi) * region_poly_lon
            py = 110.574 * region_poly_lat
        else:
            raise Exception("dist.unit argument must be either degree or km.")
        coords  = [(x, y) for x,y in zip(px, py)]
        region_win = Polygon(coords) #  box(px.min(), py.min(), px.max(), py.max()) # region.win = spatstat::owin(poly=list(x=px, y=py))
        return dict(x=x, y=y, region_win=region_win)
    
    
    @classmethod
    def xy2longlat(cls, x, y, region_poly, dist_unit="degree"):
        if len(x) != len(y):
            raise Exception("x and y must be of equal length.")
        region_poly_lon = np.array(region_poly['lon'])
        region_poly_lat = np.array(region_poly['lat'])
        if dist_unit == "degree":
            coords = [(x, y) for x,y in zip(region_poly_lon, region_poly_lat)]
            region_win = Polygon(coords)
            region_cnt = region_win.centroid.xy
            long = x / np.cos(region_cnt[1][0] * np.pi / 180) + region_cnt[0][0]
            lat = y + region_cnt[1][0]
        elif dist_unit == "km":
            lat = y / 110.574
            long = x / (111.320 * np.cos(lat/180 * np.pi))
        else:
            raise Exception("dist.unit argument must be either degree or km.")
        return dict(long=long, lat=lat)



    def __str__(self):
        string = "earthquake catalog:\n  time begin "+str(self.time_begin)+ \
                  "\n  study period: "+str(self.study_start)+ \
                  " to "+str(self.study_end)+" (T = "+ \
                  str(self.rtperiod['study_end']-self.rtperiod['study_start'])+"days)"+ \
                  "\ngeographical region:\n "+str(self.region_win)+ \
                  "\n polygonal with vertices:\n "+ str(self.region_poly)+ \
                  "\nthreshold magnitude: " +str(self.mag_threshold)+ \
                  "\nfault geometry: "+str([self.data["publicid"][f] for f, fa in 
                                            enumerate(self.revents['fault']) if fa is not None])+ \
                  "\nnumber of events:\n  total events "+ str(self.revents.shape[0]) + \
                  ": "+str((self.revents['flag']==1).sum())+" target events, " + \
                  str((self.revents['flag']!=1).sum())+" complementary events\n  (" + \
                  str((self.revents['flag']==0).sum())+" events outside geographical region, " + \
                  str((self.revents['flag']==-2).sum())+" events outside study period)\n"
        return string

    def __repr__(self):
        return self.__str__()


# plot.catalog = function(x, ...)
# {
#   oldpar = par(no.readonly = True)
#   lymat = matrix(c(1, 1, 2, 1, 1, 3, 4, 5, 6), 3, 3)
#   layout(lymat)
#   par(mar=c(4, 4.25, 1, 1))
#   plot(x$longlat.coord$long, x$longlat.coord$lat, xlab="long", ylab="lat",
#        col=8, cex=2 * (x$revents[, 4] + 0.1)/max(x$revents[, 4]),
#        asp=True, axes=False)
#   maps::map('world', add=True, col="grey50")
#   axis(1); axis(2)
#   ok = x$revents[, 5] == 1
#   points(x$longlat.coord$long[ok], x$longlat.coord$lat[ok], col=4,
#          cex=2 * (x$revents[ok, 4] + 0.1)/max(x$revents[ok, 4]))
#   polygon(x$region.poly$long, x$region.poly$lat, border=2)
#   #
#   mbk = seq(0, max(x$revents[, 4]), 0.1) + x$mag.threshold
#   mct = cut(x$revents[, 4] + x$mag.threshold, mbk)
#   mcc = log10(rev(cumsum(rev(table(mct)))))
#   plot(mbk[-length(mbk)], mcc, type="b",
#        xlab="mag", ylab=expression(log[10]*N[mag]), axes=False)
#   graphics::abline(stats::lm(mcc ~ mbk[-length(mbk)]), col=4, lty=4)
#   graphics::axis(1); graphics::axis(2)
#   #
#   tbk = seq(0, max(x$revents[, 1]), l=100)
#   tct = cut(x$revents[, 1], tbk)
#   tcc = (cumsum(table(tct)))
#   plot(tbk[-length(tbk)], tcc, type="l",
#        xlab="time", ylab=expression(N[t]), axes=False)
#   tok = (tbk[-length(tbk)] >= x$rtperiod[1]) & (tbk[-length(tbk)] <= x$rtperiod[2])
#   graphics::abline(stats::lm(tcc[tok] ~ tbk[-length(tbk)][tok]), col=4, lty=4)
#   graphics::axis(1); graphics::axis(2)
#   abline(v=x$rtperiod[1], col=2, lty=2)
#   abline(v=x$rtperiod[2], col=2, lty=2)
#   #
#   plot(x$revents[, 1], x$revents[, 3], xlab="time", ylab="lat",
#        cex=2 * (x$revents[, 4] + 0.1)/max(x$revents[, 4]), col=8, axes=False)
#   points(x$revents[ok, 1], x$revents[ok, 3], col=4,
#          cex=2 * (x$revents[ok, 4] + 0.1)/max(x$revents[ok, 4]))
#   axis(1); axis(2)
#   #
#   plot(x$revents[, 1], x$longlat.coord$long, xlab="time", ylab="long",
#        cex=2 * (x$revents[, 4] + 0.1)/max(x$revents[, 4]), col=8, axes=False)
#   points(x$revents[ok, 1], x$longlat.coord$long[ok], col=4,
#          cex=2 * (x$revents[ok, 4] + 0.1)/max(x$revents[ok, 4]))
#   axis(1); axis(2)
#   #
#   plot(x$revents[, 1], x$revents[, 4] + x$mag.threshold, xlab="time",
#        ylab="mag", pch=20, cex=0.5, col=8, axes=False)
#   points(x$revents[ok, 1], x$revents[ok, 4] + x$mag.threshold, col=4, pch=20, cex=0.5)
#   axis(1); axis(2)
#   #
#   layout(1)
#   par(oldpar)
# }


#%%

if __name__ == "__main__":
    
    df = pd.read_csv('japan_quakes.csv')
    df['eventID'] = df['Unnamed: 0'].apply(lambda x: str(x))
    df['date'] = pd.to_datetime(df['date'])
    df['year'] = df['date'].apply(lambda x: x.year)
    df['month'] = df['date'].apply(lambda x: x.month)
    df['day'] = df['date'].apply(lambda x: x.day)
    df['time'] = pd.to_datetime(df['time'])
    df['hour'] = df['time'].apply(lambda x: x.hour)
    df['minute'] = df['time'].apply(lambda x: x.minute)
    df['second'] = df['time'].apply(lambda x: x.second)
    df['longitude'] = df['long']
    df['latitude'] = df['lat']
    df['magnitude'] = df['mag']
    df['depth'] = -df['depth']
    
    # df = df.append(df.iloc[1])
    # df['year'].iat[-1] += 1
    
    keys = ['eventID','year', 'month', 'day', 'hour', 'minute', 'second',
            'longitude', 'latitude', 'depth', 'magnitude']
    data = dict()
    for key in keys:
        data[key] = df[key].to_numpy()
    
    cat = CatalogueEtas(data)
    print(cat)

