#!/usr/bin/python

import sys
import xml.dom.minidom as minidom

#-----------------------------------------------------------------
# Routine to read TMMC meta data from the specified XML file
def Parse_Standard_XML_data(xml_file):
    dom = minidom.parse( xml_file )
    temperature = float(dom.getElementsByTagName('temperature')[0].firstChild.data)
    log_activity = float(dom.getElementsByTagName('log_activity')[0].firstChild.data)
    volume = float(dom.getElementsByTagName('volume')[0].firstChild.data)
    fileprefix = dom.getElementsByTagName('file_prefix')[0].firstChild.data
    units_type = dom.getElementsByTagName('units_type')[0].firstChild.data

    if units_type == 'Reduced':
        return temperature, log_activity, volume, fileprefix, units_type
    elif units_type == 'Absolute':
        LJ_sigma = float(dom.getElementsByTagName('LJ_sigma')[0].firstChild.data)
        LJ_epsilon_kB = float(dom.getElementsByTagName('LJ_epsilon_kB')[0].firstChild.data)
        return temperature, log_activity, volume, fileprefix, units_type, LJ_sigma, LJ_epsilon_kB
    else:
        print 'Error: metadata.xml specifies units type: '+units_type
        print 'Unknown units type'
        sys.exit()

#-----------------------------------------------------------------
