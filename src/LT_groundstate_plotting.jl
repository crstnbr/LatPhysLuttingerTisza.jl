################################################################################
#
#   plotting of LTGroundstate objects
#
################################################################################

function plotLTGroundstateKSpace(
    gs :: LTGroundstate,
    bz :: BZ,
    ;
    plot_color_valid = "b",
    plot_color_invalid = "r",
    plot_invalid :: Bool = true, 
    figsize :: Tuple = (6,6),
    showPlot :: Bool = true,
    save_filename :: String = "NONE",
    constraint :: Float64 = 1e-5
    ) where {R, BZ<:AbstractBrillouinZone{R}}

    ###########################
    #   INITIAL SETTINGS
    ###########################

    # configure plot environment
    rc("font", family="serif")

    # create a new figure
    fig = figure(figsize=figsize)

    ###########################
    #   PLOT BRILLOUIN ZONE
    ###########################

    plot(bz)

    ###########################
    #   PLOT POINTS
    ###########################

    # scatter the points
    cv              = gs.constraint_values
    ks              = gs.k_vectors
    N_points        = length(cv)
    d_spatial       = size(ks)[2]
    invalid_indices = collect(1 : N_points)[cv .>= constraint]
    valid_indices   = collect(1 : N_points)[cv  .< constraint]

    # if option plot_invalid == true plot invalid points
    if plot_invalid
        if length(invalid_indices) > 0
            if d_spatial == 2
                scatter(
                [ks[i, 1] for i in invalid_indices],
                [ks[i, 2] for i in invalid_indices],
                color = plot_color_invalid,
                label = "constraint NOT fullfilled"
                )
            elseif d_spatial == 3
            scatter3D(
                [ks[i, 1] for i in invalid_indices],
                [ks[i, 2] for i in invalid_indices],
                [ks[i, 3] for i in invalid_indices],
                color = plot_color_invalid,
                label = "constraint NOT fullfilled"
                )
            end
        end
    end

    # always plot valid points
    if length(valid_indices) > 0
        if d_spatial == 2
            scatter(
            [ks[i, 1] for i in valid_indices],
            [ks[i, 2] for i in valid_indices],
            color = plot_color_valid,
            label = "constraint fullfilled"
            )
        elseif d_spatial == 3
        scatter3D(
            [ks[i, 1] for i in valid_indices],
            [ks[i, 2] for i in valid_indices],
            [ks[i, 3] for i in valid_indices],
            color = plot_color_valid,
            label = "constraint fullfilled"
            )
        end
    end

    # add a legend
    legend()

    ###########################
    #   CONFIGURE AXIS & TITLE
    ###########################

    # set the title
    title("E_min = $(gs.E_min)")
    
    # equal axis
    gca()[:set_aspect]("equal")

    ###########################
    #   FINISH THE PLOT
    ###########################

    # tighten the layout
    tight_layout()

    # save the plot
    if save_filename != "NONE"
    # make sure the directory exists
    if contains(save_filename, "/")
        # get the containing folder
        folder = save_filename[1:findlast(save_filename, '/')]
        # build the path to that folder
        mkpath(folder)
    end
    # save the plot
    savefig(save_filename)
    end

    # maybe show the plot
    if showPlot
    show()
    end

    # return the figure object
    return fig

end
