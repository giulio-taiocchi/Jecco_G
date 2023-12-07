
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  64,
    u_outer_min      =  0.1,
    u_outer_max      =  1.5,
    u_outer_domains  =  2,
    u_outer_nodes    =  24,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS4_3_1.QNM_1D(
    energy_dens  = 1.0,
    AH_pos = 1.32,
)

evoleq = AffineNull(
    #gaugecondition = ConstantAH(u_AH = 1.0),
)

io = InOut(
    out_boundary_every  = 10,
    out_gauge_every     = 10,
)

integration = Integration(
    dt              = 0.002,
    tmax            = 30.0,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
