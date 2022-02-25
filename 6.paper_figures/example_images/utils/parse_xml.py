"""
Functions to navigate Perkin Elmer HarmonyV4 XML files
"""

import os
import re
import pathlib
import xml.etree.ElementTree as ET


def parse_well(well):
    """
    Convert traditional well id nomenclature (e.g. A01) to XML expectations

    Parameters:
    -----------
    well - str
        A length three string indicating traditional well nomenclature

    Returns:
    --------
    row, col
        A tuple of Perkin Elmer row and column id (e.g. A = 1, B = 2, ...)
    """
    col = str(int(well[1:]))
    row = str(ord(well[0].lower()) - 96)
    return row, col


class PerkinElmerXML:
    """Parse Perkin Elmer Harmony XML file"""

    def __init__(self, index_xml_file):
        """Initialize class

        Arguments:
        ----------
        index_xml_file - str or pathlib.Path
            The relative or absolute file path of the xml file of interest

        Attributes:
        -----------
        tree - xml.etree.ElementTree
            XML data of the plate of interest
        treeroot - xml.etree.ElementTree.root
            An object to interact with XML tree elements
        image_element - xml.etree.ElementTree.object
            The XML location pointing to all image metadata info
        """
        self.index_xml_file = pathlib.Path(index_xml_file)

        # Parse xml
        self.tree = ET.parse(self.index_xml_file)
        self.treeroot = self.tree.getroot()

        # Get the image element object
        self.image_element = self.treeroot.find(
            "{http://www.perkinelmer.com/PEHH/HarmonyV4}Images"
        )

    def get_well_image_ids(self, well):
        """Pull identifiers for all images taken of the specific well

        Arguments:
        ----------
        well - str
            A length three string indicating traditional well nomenclature

        Returns:
        --------
        The ids of all images per well
        """
        # Parse well into PerkinElmer format
        well_row, well_col = parse_well(well)

        # Extract well elements from the XML
        well_elements = self.treeroot.find(
            "{http://www.perkinelmer.com/PEHH/HarmonyV4}Wells"
        ).findall("{http://www.perkinelmer.com/PEHH/HarmonyV4}Well")
        for well_element in well_elements:
            well_elem_id = well_element.find(
                "{http://www.perkinelmer.com/PEHH/HarmonyV4}id"
            ).text
            well_elem_row = well_element.find(
                "{http://www.perkinelmer.com/PEHH/HarmonyV4}Row"
            ).text
            well_elem_col = well_element.find(
                "{http://www.perkinelmer.com/PEHH/HarmonyV4}Col"
            ).text

            # Select the well image ids that align with the given well row and col
            if (well_elem_row == well_row) & (well_elem_col == well_col):
                well_elem_image_ids = well_element.findall(
                    "{http://www.perkinelmer.com/PEHH/HarmonyV4}Image"
                )
                well_elem_image_ids = [x.get("id") for x in well_elem_image_ids]
                break

        return well_elem_image_ids

    def _parse_image(self, im_elem):
        """Internal function to extract individual image ids.
        Note these ids are different than well image ids.

        Arguments:
        ----------
        im_elem - Individual image element in the "Images" treeroot

        Returns:
        --------
        im_id, im_url, im_channel
            The image ID, image URL (the image url is what we need to download), and channel
        """
        im_id = im_elem.find("{http://www.perkinelmer.com/PEHH/HarmonyV4}id").text
        im_url = im_elem.find("{http://www.perkinelmer.com/PEHH/HarmonyV4}URL").text
        im_channel = im_elem.find(
            "{http://www.perkinelmer.com/PEHH/HarmonyV4}ChannelName"
        ).text
        return im_id, im_url, im_channel

    def extract_image_urls(self, well_target_ids):
        """Map target ids (derived from get_well_image_ids()) to image urls for download

        Arguments:
        ----------
        well_target_ids - list
            list of image identifiers

        Returns:
        --------
        target_urls
            dictionary mapping image ids to urls
        """
        target_urls = {}
        counter = 0
        for im_element in self.image_element:
            im_id, im_url, im_channel = self._parse_image(im_element)
            if im_id in well_target_ids:
                counter += 1
                target_urls[im_id] = [im_url, im_channel]
                if counter == len(well_target_ids):
                    break
        return target_urls


def get_resolution(xml_path, tiff_file):
    """Pull the image resolution from the XML file

    Arguments:
    ----------
    xml_path - str
        Absolute or relative path to the index xml file
    tiff_file - str
        Absolute or relative path to an example tiff file

    Returns:
    --------
    resolution
        A float indicating the resolution the image was captured using
    """
    tree = ET.parse(xml_path)
    root = tree.getroot()
    for element in root:
        for sub_element in element:
            for sub_sub_element in sub_element:
                tag = re.sub("[\{].*?[\}]", "", sub_sub_element.tag)
                if tag == "URL":
                    tiff_file_tag = sub_sub_element.text
                    if tiff_file_tag == os.path.basename(tiff_file):
                        tiff_attributes = sub_element

    for attrib in tiff_attributes:
        attrib_strip = re.sub("[\{].*?[\}]", "", attrib.tag)
        if attrib_strip == "ImageResolutionX":
            resolution = float(attrib.text)
            unit = attrib.get("Unit")
            if unit == "m":
                resolution = resolution * 1000000
    return resolution
