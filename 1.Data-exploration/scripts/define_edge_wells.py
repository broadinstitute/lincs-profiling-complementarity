"""
Determine edge wells
"""

import string


def get_edge_wells():
    first_row = "A"
    num_rows = 16

    first_col = 1
    num_cols = 24

    first_col_edge = [f"{x}{first_col:02d}" for x in string.ascii_uppercase[0:num_rows]]
    last_col_edge = [f"{x}{num_cols:02d}" for x in string.ascii_uppercase[0:num_rows]]

    first_row_edge = [f"{first_row}{x:02d}" for x in range(first_col, num_cols + 1)]
    last_row_edge = [
        f"{string.ascii_uppercase[num_rows-1]}{x:02d}"
        for x in range(first_col, num_cols + 1)
    ]

    edge_wells = list(
        set(first_col_edge + last_col_edge + first_row_edge + last_row_edge)
    )

    return edge_wells
