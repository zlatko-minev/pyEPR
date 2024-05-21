# -*- coding: utf-8 -*-
"""

Debug purpose only
@author: Zlatko Minev
"""

import win32com.client
import win32com.client.combrowse

# win32com.client.combrowse.main()


# A tree heading for registered type libraries"
c = win32com.client.combrowse.HLIHeadingRegisterdTypeLibs()

for s in c.GetSubList():
    # print(s)
    name = s.GetText()
    if ("ansys" in name.lower()) or ("hfss" in name.lower()):
        print(name)
        # HFSSAppDLL 1.0 Type Library
        # C:\Program Files\AnsysEM\AnsysEM17.0\Win64\HfssDesktop.tlb
