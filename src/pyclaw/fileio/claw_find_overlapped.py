import numpy as np


def set_overlapped_status(sol):
    """
    return a list, overlapped_states,
    whose entries denote overlapped status for each patch.

    @type sol:  pyclaw.Solution
    @param sol: Solution obejct of pyclaw that contains all information
                of this time step.
    @rtype:     list
    @return:    add a component to the solution q,
                which contains overlapped status of each patch

    """
    levels = [state.patch.level-1 for state in sol.states]
    # shift base level to 0
    level_count = {}
    level_spacing = {}  # spacing of each level
    for i, level in enumerate(levels):
        if level in level_count.keys():
            level_count[level] = level_count[level] + 1
        else:
            level_count[level] = 1
            spacing = sol.states[i].patch.delta
            spacing.append(spacing[0])  # dz = dx
            spacing = np.array(spacing)
            level_spacing[level] = spacing
    # num_levels = len(level_count.keys())

    # a list of num of patches at each level
    box_per_level = [item[1] for item in
                     sorted(level_count.items(),
                            key=lambda a: a[0])]
    box_per_level = np.array(box_per_level)

    """
    @type level_count:       dictionary
    @variable level_count:   a dictionary that maps levels
                             to number of patches of certain levels
                             e.g. {0:1, 1:2, 2:12}
    @type num_levels:        int
    @variable num_levels:    number of levels in total
    @type box_per_level:     list
    @variable box_per_level: [number of patches on level 0,
                              number of patches on level1, ...]
    """
    for state in sol.states:
        level = state.patch.level-1
        xlower_coarse = state.patch.dimensions[0].lower
        # xupper_coarse = state.patch.dimensions[0].upper
        ylower_coarse = state.patch.dimensions[1].lower
        # yupper_coarse = state.patch.dimensions[1].upper
        dx = state.patch.delta[0]
        dy = state.patch.delta[1]
        nx = state.patch.num_cells_global[0]
        ny = state.patch.num_cells_global[1]
        # in overlapped_status, entry with value 0 denotes
        # that the cell is not overlapped
        # entry with value 8 denotes that the cell is overlapped
        overlapped_status = np.zeros((1, state.q.shape[1], state.q.shape[2]))
        # convert from Fortran-Style to C-style
        # overlapped_status = overlapped_status.transpose(0, 2, 1)
        # In the future, efficiency of this part can be improved
        # by mapping grid levels to
        # a list of states of corresponding levels.
        # Otherwise, we need to scan each states in each outer loop as below
        for state_fine in sol.states:
            # find states with grid level of one higher
            if ((state_fine.patch.level-1) == level + 1):
                xlower_fine = state_fine.patch.dimensions[0].lower
                xupper_fine = state_fine.patch.dimensions[0].upper
                ylower_fine = state_fine.patch.dimensions[1].lower
                yupper_fine = state_fine.patch.dimensions[1].upper
                x_idx_lower = \
                    max(int(round((xlower_fine - xlower_coarse) /
                            float(dx))), 0)
                x_idx_upper = \
                    min(int(round((xupper_fine-xlower_fine)/dx)) +
                        x_idx_lower, nx)
                y_idx_lower = \
                    max(int(round((ylower_fine - ylower_coarse) /
                            float(dy))), 0)
                y_idx_upper = \
                    min(int(round((yupper_fine-ylower_fine)/dy)) +
                        y_idx_lower, ny)
                # set these cells to 8
                overlapped_status[0, x_idx_lower:x_idx_upper,
                                  y_idx_lower:y_idx_upper].fill(8)

            else:
                continue
        # state.q is in fortran style
        # state.q = np.vstack((state.q, overlapped_status.transpose(0, 2, 1)))
        state.q = np.vstack((state.q, overlapped_status))
