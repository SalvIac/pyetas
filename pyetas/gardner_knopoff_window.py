# -*- coding: utf-8 -*-
"""
@author: Salvatore
"""

import numpy as np
from scipy import interpolate
from openquake.hmtk.seismicity.declusterer.distance_time_windows import (
    GardnerKnopoffWindow, BaseDistanceTimeWindow, time_window_cutoff,
    GruenthalWindow, UhrhammerWindow)


DAYS = 364.75

class GardnerKnopoffWindowOrig(BaseDistanceTimeWindow):
    """
    Gardner Knopoff method for calculating distance and time windows
    from Table 1 of the original paper
    https://watermark.silverchair.com/BSSA0640051363.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAoswggKHBgkqhkiG9w0BBwagggJ4MIICdAIBADCCAm0GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMmu64YM1zl6-4TZluAgEQgIICPp7S20FDHua2N40JGi8AOnBAAC0HTEGmbc1sDqtPkYa6ibAOQ-8SGH9WY_MZ1MZbCMaQ1-DsV1vsJnOY6lc-CNm3J125rb8jTG3iu73QqdHZtUR-YKn8N0zSz-4rpJaY8SEGrOIO6f15adMAQVj-Rc00zvnpdw7VfkbhcRv6mgJiMVuruDJEcrVHsRtBeIYDRuUNFb-MyFS0yYMjpkxhxQKzhSRcL1ujju1JhGS0kADJ2t8NsBN6m0yPm0LAXo3oblO_ayXsB6WVa4GVxOCLQZjqpCnkEOwuLLT0A4EKe6Xm80ACbqhB4Xq9jLlZuwagIvV78uGXCXcxH9DN_UndPpHDMXVnoWnkZ1DJ_H8V25aUFomNrZUF3vKNUWGtpWGF5U3SL7vYT7u8ehwFRg1KUdV2t9-hqOc0rseE9E_y9jVC5SllOIg1dUrQHjCTi0KykuUak1V73VYczWB9nUIVqVjCLy6K7K_nRpC90LR1tesj4LRvPIh_5Hn-NqvQgEUIVObKXNfZiezxJxccBBlWCGkFdM9UNHo6xVTp65bgVTGWuOP5LJmIqCKaF4fRwQ-LBvYIfiwLAqjR-3ZEo0IDcWu1yZYJUDpti-841ybJg1pGdZxFcIIjQTag_EUWRTh3SmSpTBqHFbo3QBllPNspdBSQ0DyLTre-lug2kuh-md5uHR8TEY02-3lR1SkmY1pe37SBNUTtJbHS4vM7d9pxgWm1d9_yKKRM0OU87BfULhmdfCWMwou1LFz5IOSbZ1M
    """
    
    m = np.arange(1, 9.5, 0.5)
    l = [19.5, 19.5, 19.5, 19.5, 22.5, 16, 30, 35, 40, 47, 54, 61, 70, 81, 94, 94, 94] # km
    t = [6, 6, 6, 6, 11.5, 22, 42, 83, 155, 290, 510, 790, 915, 960, 985, 985, 985] # days
    
    interpolant_time = interpolate.interp1d(m, t, kind="cubic")
    interpolant_space = interpolate.interp1d(m, l, kind="cubic")
    
    def calc(self, magnitude, time_cutoff=None):
        sw_space = self.interpolant_space(magnitude)
        sw_time = self.interpolant_time(magnitude) / DAYS
        if time_cutoff:
            sw_time = time_window_cutoff(sw_time, time_cutoff)
        return sw_space, sw_time


if __name__ == "__main__":
    
    
    import matplotlib.pyplot as plt
    
    gw = GruenthalWindow()
    oq_gw = gw.calc(np.linspace(2, 8, 100))[1]*364.75
    uw = UhrhammerWindow()
    oq_uw = uw.calc(np.linspace(2, 8, 100))[1]*364.75
    oq = GardnerKnopoffWindow()
    oq_gk = oq.calc(np.linspace(2, 8, 100))[1]*364.75

    orig = GardnerKnopoffWindowOrig()
    or_gk = orig.calc(np.linspace(2, 8, 100))[1]*364.75

    fig = plt.figure()
    # oq routines
    plt.plot(np.linspace(2, 8, 100), oq_gk, label="GardnerKnopoffWindow")
    plt.plot(np.linspace(2, 8, 100), oq_gw, label="GruenthalWindow")
    plt.plot(np.linspace(2, 8, 100), oq_uw, label="UhrhammerWindow")
    # modified routine
    plt.plot(np.linspace(2, 8, 100), or_gk, label="Modified GardnerKnopoffWindow")
    # original
    plt.scatter(np.arange(2.5,8.5,0.5), [6,11.5,22,42,83,155,290,510,790,915,960,985],
                label="Original Gardner Knopoff table", color="k")
    plt.legend()
    plt.show()
    

