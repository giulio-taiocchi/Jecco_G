
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -5.0,
    x_max            =  5.0,
    x_nodes          =  32,
    y_min            = -5.0,
    y_max            =  5.0,
    y_nodes          =  32,
    u_outer_min      =  0.2,
    u_outer_max      =  1.03,
    u_outer_domains  =  1,
    u_outer_nodes    =  36,
    u_inner_nodes    =  12,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS4_3_1.BoostedBBnumerical(
    AH_pos = 1.001,
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.00),
    
)

io = InOut(
    out_boundary_every  = 1,
    out_bulk_every      = 1,
    #out_gauge_every     = 10,
    remove_existing     = true,
)

integration = Integration(
#0.002
    dt              = 0.01,
    tmax            = 50.0,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
