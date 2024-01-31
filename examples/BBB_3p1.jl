
using Jecco.AdS4_3_1

grid = SpecCartGrid3D(
    x_min            = -1.0,
    x_max            =  1.0,
    x_nodes          =  24,
    y_min            = -1.0,
    y_max            =  1.0,
    y_nodes          =  24,
    u_outer_min      =  0.1,
    u_outer_max      =  1.0,
    u_outer_domains  =  1,
    u_outer_nodes    =  128,
    u_inner_nodes    =  36,
    fd_order         =  4,
    sigma_diss       =  0.2,
)


id = AdS4_3_1.BoostedBBnumerical(
    AH_pos = 0.85,
    a3_ampx	= -3/4,
    a3_translx	= -7/4,
    A = -1.0,
    B = 1.0,
    phase_a = 1.0,
    phase_fx = 0.5,
)

evoleq = AffineNull(
    gaugecondition = ConstantAH(u_AH = 1.00),# this u_AH is a not used parameter
    
)

io = InOut(
    out_boundary_every  = 1,
    out_bulk_every      = 1,
    out_bulkconstrained_every = 1,
    #out_gauge_every     = 1,
    remove_existing     = true,
)

integration = Integration(

    dt              = 0.002,
    tmax            = 50,
    ODE_method      = AdS4_3_1.AB3(),
    filter_poststep = true,
)

run_model(grid, id, evoleq, integration, io)
