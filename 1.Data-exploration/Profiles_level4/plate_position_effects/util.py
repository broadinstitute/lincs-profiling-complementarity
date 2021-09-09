import argparse
import pathlib
import pandas as pd
import numpy as np


def increment(pos_id, diffusion, row=True, reverse=True):
    if row:
        _id = ord(pos_id)
    else:
        _id = int(pos_id)

    if reverse:
        update = _id - diffusion
    else:
        update = _id + diffusion

    if row:
        return chr(update)
    else:
        return f"{update:02}"


def compile_positions(pos_id, diffusion, row=True, reverse=False):
    add_pos_id = str(pos_id)
    pos_ids = [add_pos_id]
    for _ in range(1, diffusion + 1):
        add_pos_id = increment(add_pos_id, 1, row=row, reverse=reverse)
        pos_ids.append(add_pos_id)

    return pos_ids


def mirror_data(row_or_col, all_row_or_col_list):
    row_or_col_id_idx = all_row_or_col_list.index(row_or_col)
    return all_row_or_col_list[::-1][row_or_col_id_idx]


def get_mirror(position, all_positions):
    row = position[0]
    col = position[1:]

    allrows = [x[0] for x in all_positions]
    allcols = [x[1:] for x in all_positions]

    mirror_row = mirror_data(row, allrows)
    mirror_col = mirror_data(col, allcols)

    return f"{mirror_row}{mirror_col}"


def get_compatible_positions(
    positions, pos_id, diffusion, min_row, max_row, min_col, max_col
):
    positions = list(positions)

    row = pos_id[0]
    col = pos_id[1:]

    compatible_rows = compile_positions(row, diffusion=diffusion)
    compatible_rows += compile_positions(row, diffusion=diffusion, reverse=True)
    compatible_cols = compile_positions(col, diffusion=diffusion, row=False)
    compatible_cols += compile_positions(
        col, diffusion=diffusion, row=False, reverse=True
    )

    compatible_rows = [x for x in compatible_rows if min_row <= x <= max_row]
    compatible_cols = [x for x in compatible_cols if min_col <= x <= max_col]

    compatible_ids = list(
        set(
            [
                f"{row_id}{col_id}"
                for row_id in compatible_rows
                for col_id in compatible_cols
            ]
        )
    )

    compatible_ids = [x for x in compatible_ids if x in positions]

    return compatible_ids


def corners(positions, diffusion=1, mirror=False, summarize=False):
    row = [x[0] for x in positions]
    col = [x[1:] for x in positions]

    top_row = min(row)
    bot_row = max(row)

    left_col = min(col)
    right_col = max(col)

    top_rows = compile_positions(top_row, diffusion=diffusion)
    bot_rows = compile_positions(bot_row, diffusion=diffusion, reverse=True)
    left_cols = compile_positions(left_col, diffusion=diffusion, row=False)
    right_cols = compile_positions(
        right_col, diffusion=diffusion, row=False, reverse=True
    )

    top_left_corner = [
        f"{row_id}{col_id}" for row_id in top_rows for col_id in left_cols
    ]
    bot_left_corner = [
        f"{row_id}{col_id}" for row_id in bot_rows for col_id in left_cols
    ]
    top_right_corner = [
        f"{row_id}{col_id}" for row_id in top_rows for col_id in right_cols
    ]
    bot_right_corner = [
        f"{row_id}{col_id}" for row_id in bot_rows for col_id in right_cols
    ]

    what_is_left_internal = list(
        set(positions).difference(
            set(top_left_corner + bot_left_corner + top_right_corner + bot_right_corner)
        )
    )

    if mirror:
        return_group = {
            "set1": top_left_corner + bot_right_corner,
            "set2": top_right_corner + bot_left_corner,
        }

    else:
        return_group = {
            "top_left_corner": top_left_corner,
            "bot_left_corner": bot_left_corner,
            "top_right_corner": top_right_corner,
            "bot_right_corner": bot_right_corner,
        }

    if summarize:
        set_df = pd.DataFrame(return_group).melt(
            var_name="group_set", value_name="well_position"
        )
        full_well_df = pd.DataFrame(positions, columns=["well_position"])
        summary_df = full_well_df.merge(set_df, on="well_position", how="outer").fillna(
            "interior"
        )
        summary_df = summary_df.assign(
            row=[x[0] for x in summary_df.well_position],
            col=[x[1:] for x in summary_df.well_position],
        )

        summary_df.row = pd.Categorical(
            summary_df.row, categories=summary_df.row.unique()[::-1], ordered=True
        )
        summary_df.col = pd.Categorical(summary_df.col, ordered=True)

    return_group["interior"] = what_is_left_internal

    if summarize:
        return return_group, summary_df
    else:
        return return_group


def diffuse_wells(positions, diffusion, mirror=False, keep_same_position=False):
    row = [x[0] for x in positions]
    col = [x[1:] for x in positions]

    min_row = min(row)
    max_row = max(row)

    min_col = min(col)
    max_col = max(col)

    diffuse_set = {}
    for pos_id in positions:
        compatible_ids = get_compatible_positions(
            positions=positions,
            pos_id=pos_id,
            diffusion=diffusion,
            min_row=min_row,
            max_row=max_row,
            min_col=min_col,
            max_col=max_col,
        )

        if mirror:
            mirror_id = get_mirror(pos_id, positions)
            compatible_ids += get_compatible_positions(
                positions=positions,
                pos_id=mirror_id,
                diffusion=diffusion,
                min_row=min_row,
                max_row=max_row,
                min_col=min_col,
                max_col=max_col,
            )

        if keep_same_position:
            return_compatible_ids = list(set([x for x in compatible_ids]))
        else:
            return_compatible_ids = list(
                set([x for x in compatible_ids if x != pos_id])
            )

        diffuse_set[pos_id] = return_compatible_ids

    return diffuse_set


def load_args():
    parser = argparse.ArgumentParser(description="Parse arguments")
    parser.add_argument(
        "--data_dir", type=str, help="directory that contains the profile data"
    )
    parser.add_argument(
        "--output_dir",
        type=str,
        help="directory to save output files",
        default="results/diffusion",
    )
    parser.add_argument(
        "--profile_file", type=str, help="the file name of the profile file"
    )
    parser.add_argument(
        "--diffusion",
        type=int,
        help="how many neighboring wells to consider as non-replicates",
    )
    parser.add_argument(
        "--mirror",
        action="store_true",
        help="Whether to consider mirrored wells as non-replicates (if flipped plate)",
    )
    parser.add_argument(
        "--drop_same_position",
        action="store_true",
        help="Whether to drop the same well position from comparison",
    )
    parser.add_argument(
        "--l1000",
        action="store_true",
        help="Boolean if the data are L1000 (different column names)",
    )
    return parser.parse_args()
