"""
Helper functions to retrieve image plotting constants
"""

import numpy as np


def get_channel_map():
    """Helper function to load channel map

    Returns:
    channel_map - dict
        Dictionary mapping channels to cell compartments
    """

    channel_map = {
        'HOECHST 33342': "DNA",
        'Alexa 488': "ER",
        '488 long': "RNA",
        'Alexa 568': "AGP",
        'Alexa 647': "Mito"
    }

    return channel_map


def get_channel_colors():
    """Helper function to load channel colors for colorization

    Returns:
    channel_colors - dict
        Dictionary of compartments mapping to RGB channel arrays
    """

    channel_colors = {
        "DNA": ["blue", np.array([0, 0, 255], dtype=np.uint8)],
        "ER": ["green", np.array([0, 255, 0], dtype=np.uint8)],
        "RNA": ["yellow", np.array([255, 255, 0], dtype=np.uint8)],
        "AGP": ["orange", np.array([255, 150, 0], dtype=np.uint8)],
        "Mito": ["red", np.array([255, 0, 0], dtype=np.uint8)],
    }

    return channel_colors
