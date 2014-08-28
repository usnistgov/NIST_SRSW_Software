#!/usr/bin/python

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

    return temperature, log_activity, volume, fileprefix, units_type
#-----------------------------------------------------------------
