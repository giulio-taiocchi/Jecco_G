
#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_outer("B1")

should expand to

  B1t = B1_x -  (Fx + xi_x) * B1p

etc.

=#
macro tilde_outer(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx + xi_x) * $fp) )
end
macro hat_outer(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy + xi_y)  * $fp) )
end
macro bar_outer(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx + xi_x)  * ( -2*($fp_x) + (Fx + xi_x) * ($fpp) )) )
end
macro star_outer(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy + xi_y)  * ( -2*($fp_y) + (Fy + xi_y) * ($fpp) )) )
end
macro cross_outer(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx + xi_x) * ($fp_y) -
                 (Fy + xi_y)  * ( $fp_x - (Fx + xi_x)  * ($fpp) ) ) )
end


# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (u, xi, B, Bp, G, Gp) = vars
	
	ABCS[1] = 4*u^4

	ABCS[2] = 8*u^3

	ABCS[3] = Gp^2 + Bp^2*cosh(G)^2

	ABCS[4] =0


    nothing
end

# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B11 d_u Fx + B12 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
        u, xi, xi_x, xi_y,
        B     ,        G      ,        S      ,
        Bp    ,        Gp     ,        Sp     ,
        Bpp   ,        Gpp    ,        Spp    ,
        B_x   ,        G_x    ,        S_x    ,
        B_y   ,        G_y    ,        S_y    ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,
        Bp_y  ,        Gp_y   ,        Sp_y
    ) = vars

 



x0 = exp(B)
x1 = S ^ 2
x2 = 2 * x1
x3 = u ^ 4 * x2
x4 = cosh(G)
x5 = x4 ^ 2
x6 = Bp * x5
x7 = u ^ 2
x8 = x2 * x7
x9 = x0 * x8
x10 = Bp * x4 * sinh(G)
x11 = Bp ^ 2
x12 = 2 * Bp
x13 = Bpp * S + Sp * x12
x14 = S * x11 + x13
x15 = 2 * S
x16 = x15 * x5
x17 = 2 * G
x18 = sinh(x17)
x19 = x1 * x18
x20 = Gp * x12
x21 = x19 * x20
x22 = Sp ^ 2
x23 = Gp ^ 2 * x1
x24 = S * Spp
x25 = -4 * x22 + 2 * x23 + 4 * x24
x26 = S * x11 - x13
x27 = S * x20
x28 = cosh(x17)
x29 = Bp * x28
x30 = Gp * x15 * x29
x31 = 4 * Sp
x32 = Gp * x31 + Gpp * x15
x33 = x18 * x26 + x27 - x30 - x32
x34 = S * x0
x35 = 4 * x34
x36 = Gp * x2
x37 = x0 * x36
x38 = Bp * x19
x39 = S_x * x35
x40 = x2 * x5
x41 = Bp_x * x0
x42 = x0 * x2
x43 = x12 * x19
x44 = S * x5
x45 = Gp * x38
x46 = -2 * x22 + x23 + 2 * x24
x47 = x34 * (x14 * x18 - x27 + x30 - x32)
x48 = 4 * S
AA[1,1] = x0 * x3
AA[1,2] = 0
BB[1,1] = x9 * (2 * u - x6)
BB[1,2] = x8 * (Gp + x10)
CC[1,1] = x0 * (x14 * x16 + x21 + x25)
CC[1,2] = S * x33
SS[1] = -B_x * x42 * x6 - B_y * x36 - B_y * x38 + 2 * Bp * G_y * x1 * x28 + 2 * Bp * S * S_y * x18 + Bp_y * x1 * x18 - G_x * x0 * x43 - G_x * x37 + 4 * Gp * S * S_y + 2 * Gp_y * x1 + S * x33 * xi_y + 4 * S_x * Sp * x0 - Sp_x * x35 + 2 * x0 * xi_x * (x14 * x44 + x45 + x46) - x39 * x6 - x40 * x41
AA[2,1] = 0
AA[2,2] = x3
BB[2,1] = x9 * (Gp - x10)
BB[2,2] = x1 * x7 * (Bp + 4 * u + x29)
CC[2,1] = x47
CC[2,2] = x16 * x26 - x21 + x25
SS[2] = -B_x * x0 * x38 + B_x * x37 - B_y * x2 * x6 + Bp_y * x40 - G_x * x29 * x42 - G_y * x36 + G_y * x43 + Gp * x39 + Gp_x * x42 + S_y * x31 + S_y * x48 * x6 - Sp_y * x48 - x10 * x39 - x19 * x41 + x47 * xi_x + 2 * xi_y * (x26 * x44 - x45 + x46)
    
    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
          u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,        G      ,        S      ,    Fx     ,    Fy     ,
        Bp    ,        Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,        Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,        G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,        G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,        Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,        G_xx   ,        S_xx   ,
        B_yy  ,        G_yy   ,        S_yy   ,
                        G_xy   ,        S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Bp")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Bp")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")

    
   
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


    x0 = exp(B)
x1 = 8 * x0
x2 = S ^ 2
x3 = cosh(G)
x4 = 4 * x3
x5 = sinh(G)
x6 = 2 * Gh
x7 = 2 * x0 * x2
x8 = 4 * Sh
x9 = Fyh + xi_yy
x10 = 4 * Sp
x11 = Fy + xi_y
x12 = 4 * x11
x13 = x11 ^ 2
x14 = 4 * Spp
x15 = Spp * x11
x16 = 4 * St
x17 = Fxh + Fyt + 2 * xi_xy
x18 = Fx + xi_x
x19 = Spp * x18 + Spt
x20 = 2 * x11
x21 = -Fxp
x22 = Gpp * x11
x23 = 2 * x18
x24 = 2 * Fy + 2 * xi_y
x25 = Gpp * x18 + Gpt
x26 = 2 * x9
x27 = 2 * Bh
x28 = 2 * Bph
x29 = Fxt + xi_xx
x30 = 4 * x18
x31 = x18 ^ 2
x32 = 2 * x29
x33 = 2 * x31
x34 = 2 * Bt
x35 = 2 * Fx + 2 * xi_x
ABCS[1] = 0
ABCS[2] = -S ^ 3 * u ^ 2 * x1
ABCS[3] = Sp * x1 * x2
ABCS[4] = ( -12 * S ^ 4 * x0 + S * x0 * (x3 * (-Gh * x16 - Gt * x8) + x5 * (-8 * Sc + 4 * Sp * x17 + 8 * Spt * x11 - 8 * x15 * x18 - x19 * (8 * Fy + 8 * xi_y))) + S * (Gh * x5 * x8 - x3 * (Bh * x8 + Sph * x12 - 4 * Ss + x10 * x9 - x12 * (Sph + x15) - x13 * x14)) - Sh ^ 2 * x4 + Sh * St * x1 * x5 + x2 * (x3 * (2 * Bh ^ 2 + Bp * x26 - 2 * Bpp * x13 - 2 * Bs + Fyp ^ 2 + Fyp * x27 - 2 * Fyph + 2 * Gh ^ 2 + x11 * x28 - x11 * (Bpp * x20 + x28)) + x5 * (-Gp * x26 - Gph * x20 + 2 * Gpp * x13 + 2 * Gs + x24 * (Gph + x22) - x6 * (Fyp + x27))) + x3 * x7 * (-2 * Gc + Gh * (-Bt - x21) + Gp * x17 + Gpt * x20 + Gt * (Bh + Fyp) - x22 * x23 - x24 * x25) + x5 * x7 * (-Fxp * Fyp + Fxph + Fypt - Gt * x6) + (S * (Gt * x16 * x5 + x3 * (Bt * x16 + 4 * Sb - Spt * x30 - x10 * x29 - x14 * x31 + x19 * x30)) - St ^ 2 * x4 + x2 * (x3 * (2 * Bb - Bp * x32 - Bpp * x33 - Bpt * x23 + 2 * Bt ^ 2 + Fxp ^ 2 - Fxp * x34 - 2 * Fxpt + 2 * Gt ^ 2 + x35 * (Bpp * x18 + Bpt)) + x5 * (2 * Gb - Gp * x32 - Gpp * x33 - Gpt * x23 + 2 * Gt * (x21 + x34) + x25 * x35))) * exp(2 * B))

    nothing
end


# this is another coupled equation, for Bd and Gd. the notation used is
#
# ( A11 d_uu Bd + A12 d_uu Gd + B11 d_u Bd + B12 d_u Gd + C11 Bd + C12 Gd ) = -S1
# ( A21 d_uu Bd + A22 d_uu Gd + B21 d_u Bd + B22 d_u Gd + C21 Bd + C22 Gd ) = -S2

function BdGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Outer)
    (
       u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,      G      ,        S      ,    Fx     ,    Fy     ,  Sd,
        Bp    ,      Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,      Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,      G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,      G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,      Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,      Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,      G_xx   ,        S_xx   ,
        B_yy  ,      G_yy   ,        S_yy   ,
                      G_xy   ,        S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")


    
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    tanhG   = tanh(G)
    sinhG   = sinh(G)
    sechG   = sech(G)


    x0 = exp(B)
x1 = 8 * x0
x2 = S ^ 3
x3 = u ^ 2 * x2
x4 = tanh(G)
x5 = S ^ 2
x6 = x1 * x5
x7 = Bp * x2
x8 = sech(G)
x9 = Fyp ^ 2
x10 = S * x8
x11 = 4 * x0
x12 = exp(2 * B)
x13 = Fxp * St
x14 = Fxpt * S
x15 = Fxp ^ 2 * S
x16 = sinh(G)
x17 = 4 * Sh
x18 = cosh(G)


AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x1 * x3
BB[1,2] = 0
CC[1,1] = x6 * (Gp * S * x4 + Sp)
CC[1,2] = x1 * x4 * x7
SS[1] = ( Bp * Sd * x6 + 8 * Fyp * Sh * x8 + x10 * x11 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp) + x10 * (-4 * Fyph + 2 * x9) - x12 * x8 * (8 * x13 - 4 * x14 + 2 * x15))
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x11 * x3
CC[2,1] = -2 * x0 * x7 * sinh(2 * G)
CC[2,2] = Sp * x11 * x5
SS[2] = ( -Fyp * x16 * x17 + 4 * Gp * Sd * x0 * x5 - 2 * S * x0 * x18 * (Bh * Fxp - Bp * Fxh + Bp * Fyt - Bt * Fyp - Fxp * Fyp + Fxph + Fypt) - S * x16 * (-2 * Fyph + x9) + x0 * x18 * (Fxp * x17 + 4 * Fyp * St) - x12 * x16 * (4 * x13 - 2 * x14 + x15) )

    nothing
end


function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Outer)
    (
         u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd, Bd, Gd, 
        Bp  ,  Gp  ,  Sp   , Fxp   , Fyp   ,
        Bpp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,
        B_x ,  G_x ,  S_x  , Fx_x  , Fy_x  ,
        B_y ,  G_y ,  S_y  , Fx_y  , Fy_y  ,
        Bp_x,  Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        Bp_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,
        B_xx,  G_xx,  S_xx ,
        B_yy,  G_yy,  S_yy ,
                G_xy,  S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")

    @star_outer("B")
    @star_outer("G")
    @star_outer("S")

    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("Fxp")
    @hat_outer("Fyp")

    @cross_outer("G")
    @cross_outer("S")


    #expB1   = exp(B1)
    sinh2G  = sinh(*(2, G))
    cosh2G  = cosh(*(2, G))
    coshGsq = cosh(G)^2
    coshG   = cosh(G)
    sinhG   = sinh(G)


   x0 = exp(B)
x1 = S ^ 4
x2 = 2 * x1
x3 = cosh(G)
x4 = 4 * x3
x5 = sinh(G)
x6 = Sh * x5
x7 = Fy * xi_y
x8 = Fy ^ 2
x9 = 2 * Spp * x3
x10 = xi_y ^ 2
x11 = 4 * S
x12 = S ^ 2
x13 = 8 * Gpp
x14 = 4 * Gpp
x15 = 2 * Bp
x16 = 4 * Bpp
x17 = Gp * S
x18 = 2 * Sp
x19 = 2 * S
x20 = Bt * S
x21 = 2 * Gt
x22 = 2 * x12
x23 = 4 * Spp
x24 = Fx * x23
x25 = xi_x * (Fy + xi_y)
x26 = 2 * Gp
x27 = Fx * x13


ABCS[1] = u ^ 4 * x0 * x2
ABCS[2] = 4 * u ^ 3 * x0 * x1
ABCS[3] = 0
ABCS[4] = (Sh ^ 2 * x4 + x0 * (-8 * St * x6 + x11 * (x3 * (Gh * St + Gt * Sh) + x5 * (-Fxh * Sp + Fy * x24 - Fyt * Sp + 2 * Sc - x18 * xi_xy + x23 * x25 + x24 * xi_y)) + x12 * (-8 * Sd * Sp + x3 * (-Bh * x21 + 2 * Bt * Gh - Fxh * x26 + Fy * x27 - Fyt * x26 + 4 * Gc - 4 * Gp * xi_xy + x13 * x25 + x27 * xi_y) - x5 * (2 * Fxp * Fyp - 4 * Gh * Gt)) + x2 * (Bd * Bp * x3 ^ 2 + Gd * Gp)) + x11 * (Bh * Sh * x3 + Fyh * Sp * x3 - Gh * x6 + Sp * x3 * xi_yy - Spp * x4 * x7 - Ss * x3 - x10 * x9 - x8 * x9) + x12 * (x3 * (-2 * Bh ^ 2 + 8 * Bpp * x7 + 2 * Bs - Fyh * x15 + Fyp ^ 2 - 2 * Gh ^ 2 + x10 * x16 - x15 * xi_yy + x16 * x8) + x5 * (4 * Bh * Gh + 2 * Fyh * Gp + 2 * Gp * xi_yy - 2 * Gs - x10 * x14 - x13 * x7 - x14 * x8)) + (-x19 * x5 * (-Fxt * x17 + Gb * S + x21 * (St + x20)) + x19 * xi_xx * (x17 * x5 + x3 * (Bp * S + x18)) + x3 * (-Bb * x22 + 2 * Bp * Fxt * x12 - Bt ^ 2 * x22 + Fxp ^ 2 * x12 + 4 * Fxt * S * Sp - Gt ^ 2 * x22 - Sb * x11 + 4 * St ^ 2 - 4 * St * x20)) * exp(2 * B))

    nothing
end

function xi_t_eq_coeff(vars::Tuple, ::Outer)
    (
        kappa, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B   ,  G   ,  S    , Fx    , Fy    , Sd ,  Bd  , Gd,  A   ,
        Bp  ,  Gp  ,  Sp   , Fxp   , Fyp   , Sdp,  Bdp , Gdp, Ap  ,
        Bpp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,                   App ,
        B_x ,   G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x, Bd_x, Gd_x,A_x ,
	B_y ,  G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y, Bd_y,Gd_y, A_y ,
        Bp_x,  Gp_x,  Sp_x , Fxp_x , Fyp_x ,                   Ap_x,
        Bp_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,                   Ap_y,
                                 Fy_xx ,       A_xx,
                                          Fx_yy ,               A_yy,
                                          Fx_xy , Fy_xy ,       A_xy,Fx_xx,Fy_yy,
                                          B_xx,B_yy,B_xy,G_xx,G_yy,G_xy,S_xx,S_yy,S_xy
    ) = vars

    @tilde_outer("B")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")
    @tilde_outer("Bd")
    @tilde_outer("Gd")
    @tilde_outer("A")

    @hat_outer("B")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")
    @hat_outer("Bd")
    @hat_outer("Gd")
    @hat_outer("A")

    @bar_outer("B")
    @bar_outer("G")
    @bar_outer("S")
    @bar_outer("Fx")
    @bar_outer("Fy")
    #@bar_outer("Sd")
    #@bar_outer("Bd")
    #@bar_outer("Gd")
    @bar_outer("A")
    
    @tilde_outer("Sd")
    @tilde_outer("Bd")
    @tilde_outer("Gd")
    

    @star_outer("A")
    @star_outer("B")
    @star_outer("G")
    @star_outer("S")
    @star_outer("Fx")
    @star_outer("Fy")
    #@star_outer("Sd")
    #@star_outer("Bd")
    #@star_outer("Gd")
   

    @tilde_outer("Bp")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")
    @tilde_outer("Ap")

    @hat_outer("Bp")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")
    @hat_outer("Ap")

    @cross_outer("A")
    @cross_outer("Fx")
    @cross_outer("Fy")
    @cross_outer("S")
    @cross_outer("G")
    

   

    x0 = cosh(G)
x1 = 8 * x0
x2 = S ^ 4
x3 = exp(3 * B)
x4 = x2 * x3
x5 = x1 * x4
x6 = exp(B)
x7 = x2 * x6
x8 = x1 * x7
x9 = exp(2 * B)
x10 = x2 * x9
x11 = sinh(G)
x12 = 16 * x11
x13 = Gh * S
x14 = S * x6
x15 = Fxp * x14
x16 = Gp * S
x17 = 2 * Fy
x18 = Fx * Sp
x19 = 2 * x6
x20 = tanh(G)
x21 = Sp * x20
x22 = Fyp * S
x23 = Bp * S
x24 = Gt * x14
x25 = 2 * x16
x26 = Fx * x6
x27 = x25 * x26
x28 = x16 * x20
x29 = x19 * xi_x
x30 = S ^ 3
x31 = x1 * x30
x32 = x31 * x9
x33 = Bh * S
x34 = Sp * x17
x35 = 2 * xi_y
x36 = x31 * x6
x37 = 16 * x0
x38 = S * x37
x39 = Sph * x38
x40 = Fy * Sp
x41 = Sh * x37
x42 = x12 * x13
x43 = Sp * xi_y
x44 = S ^ 2
x45 = x1 * x44
x46 = x11 * x44
x47 = 8 * x46
x48 = 32 * x11
x49 = St * x6
x50 = x23 * x41
x51 = x22 * x37
x52 = Fy * x16
x53 = Sh * x12
x54 = S * x0
x55 = Spp * x54
x56 = 64 * xi_y
x57 = Fy * x56
x58 = x16 * xi_y
x59 = Bp * x45
x60 = 32 * Sd
x61 = x30 * x6
x62 = Bp * x61
x63 = 24 * Fy
x64 = x0 * x44
x65 = Bph * x64
x66 = 24 * xi_y
x67 = Fypp * x45
x68 = Gph * x46
x69 = Fy ^ 2
x70 = 32 * x69
x71 = Gp * x47
x72 = Gh * x47
x73 = xi_y ^ 2
x74 = 32 * x73
x75 = x44 * x6
x76 = Fx * Sph
x77 = x12 * x14
x78 = Sh * Sp
x79 = 16 * Fx
x80 = x11 * x6
x81 = x79 * x80
x82 = Fy * Spt
x83 = x12 * x49
x84 = x37 * x49
x85 = 16 * xi_x
x86 = x11 * x85
x87 = x6 * x86
x88 = x0 * x23
x89 = 32 * x88
x90 = 4 * x44
x91 = Fyp ^ 2
x92 = x47 * x6
x93 = x37 * x75
x94 = x37 * x9
x95 = S * x94
x96 = Bp * Fy
x97 = Fyp * x45
x98 = Sp * x23
x99 = x37 * x98
x100 = Bp * xi_y
x101 = Bpp * Fy
x102 = Fy * Gp
x103 = 24 * x46
x104 = Fyp * x103
x105 = 24 * Gh
x106 = x105 * x64
x107 = Gpp * x46
x108 = Gp * xi_y
x109 = 64 * Fy
x110 = Fx * Spp
x111 = x11 * x14
x112 = x110 * x111
x113 = Sp * x22
x114 = x6 * x79
x115 = Sp * x0
x116 = x115 * x13
x117 = Sh * x0 * x16
x118 = x12 * x15
x119 = x24 * x37
x120 = x109 * xi_x
x121 = Spp * x111
x122 = x6 * x85
x123 = x56 * xi_x
x124 = Bpp * x64
x125 = x12 * x75
x126 = Bp * x125
x127 = x46 * x6
x128 = Bp * x127
x129 = x45 * x6
x130 = Gh * x129
x131 = 8 * Fx
x132 = Fypp * x127
x133 = Fx * Gph
x134 = x0 * x75
x135 = 24 * x134
x136 = x79 * x9
x137 = Spt * x54
x138 = St * x115
x139 = Gp * x129
x140 = Fxp * Fyp
x141 = Fxpp * x92
x142 = Fy * Gpt
x143 = 8 * xi_x
x144 = x12 * x9
x145 = Gt * S
x146 = x85 * x9
x147 = 32 * x46
x148 = x23 * x6
x149 = x148 * x48
x150 = x149 * x18
x151 = Sp * xi_x
x152 = x149 * x151
x153 = x45 * x9
x154 = x47 * x9
x155 = Bp * Gp
x156 = x12 * x44
x157 = x155 * x156
x158 = Bp ^ 2
x159 = x158 * x44
x160 = x0 * x159
x161 = Fy * xi_y
x162 = Gp ^ 2
x163 = x162 * x64
x164 = 40 * x163
x165 = Fyp * x128
x166 = Bp * Fx
x167 = Gh * x135
x168 = St * x88
x169 = Fxp * x125
x170 = Gt * x129
x171 = Bp * xi_x
x172 = S * x115
x173 = Bt * x172
x174 = Bt * x129
x175 = Fxp * x172
x176 = Fx * Gpp
x177 = x134 * x176
x178 = Fx * Gp
x179 = Fyp * x135
x180 = x105 * x127
x181 = x11 * x136
x182 = St * x16
x183 = Sp * x145
x184 = Fxp * x135
x185 = 24 * Gt * x127
x186 = Gpp * x134
x187 = Gp * xi_x
x188 = x86 * x9
x189 = 3 * x160
x190 = 20 * x163
x191 = x64 * x9
x192 = 24 * x191
x193 = Bp * x192
x194 = Bpt * x192
x195 = x156 * x9
x196 = Fxpp * x191
x197 = Fx * x9
x198 = Gpt * x103
x199 = x71 * x9
x200 = x9 * xi_x
x201 = 32 * x134
x202 = x166 * x201
x203 = x171 * x201
x204 = 3 * G
x205 = x159 * cosh(x204)
x206 = Fxp ^ 2
x207 = x9 * x90
x208 = 40 * x166
x209 = Bt * x191
x210 = 40 * x171
x211 = Fxp * x191
x212 = x46 * x9
x213 = Gt * x212
x214 = Fx ^ 2
x215 = x214 * x94
x216 = xi_x ^ 2
x217 = x26 * x46
x218 = x158 * x217
x219 = x6 * xi_x
x220 = x219 * x46
x221 = x158 * x220
x222 = Bpp * Fx
x223 = x103 * x9
x224 = Bt * x223
x225 = Fxp * x223
x226 = 40 * x162
x227 = x217 * x226
x228 = Gt * x192
x229 = x220 * x226
x230 = Bpp * x44
x231 = Bpp * x216
x232 = Gpp * x216
x233 = x214 * x9
x234 = 48 * x155 * x46
x235 = x216 * x9
x236 = x159 * sinh(x204)
x237 = x236 * x26
x238 = x197 * xi_x
x239 = x219 * x236
x240 = 35 * x160
x241 = x190 * x9
x242 = x205 * x9
x243 = 2 * x0
x244 = 2 * x11
x245 = 3 * x11
x246 = Fy + xi_y
x247 = S ^ 6 * x9
x248 = S ^ 5
x249 = 4 * x247
x250 = 8 * x10
x251 = x11 * x250
x252 = 8 * x9
x253 = A * x252
x254 = Fxp * x10
x255 = 4 * x11
x256 = x10 * x255
x257 = 12 * x0
x258 = x246 ^ 2
x259 = x258 * x7
x260 = x0 ^ 2
x261 = x246 ^ 4
x262 = x0 ^ 4
x263 = x23 * x262
x264 = 4 * Sp
x265 = x246 ^ 3
x266 = A * x244
x267 = x10 * x266
x268 = Sp * x246
x269 = 2 * x261
x270 = Bp * x246
x271 = Fyp * x270
x272 = 4 * x0
x273 = x272 * x7
x274 = A * x273
x275 = A * x272
x276 = x275 * x61
x277 = Fx + xi_x
x278 = x277 ^ 2
x279 = x278 * x4
x280 = 2 * x262
x281 = x265 * x44
x282 = Bp * x281
x283 = x280 * x282
x284 = 2 * G
x285 = sinh(x284)
x286 = sinh(4 * G)
x287 = x285 ^ 2
x288 = Bpp * x258
x289 = A * x7
x290 = x243 * x289
x291 = x288 * x290
x292 = x0 ^ 3
x293 = 4 * x270 * x292
x294 = Gp * x277
x295 = x10 * x272
x296 = Ah * x295
x297 = x259 * x272
x298 = Fxp * x246
x299 = Ap * x256
x300 = Fyp * x277
x301 = 4 * x115
x302 = Ap * x301
x303 = x258 * x61
x304 = Gp * x246
x305 = At * x295
x306 = Sp * x277
x307 = x265 * x285
x308 = x144 * x30
x309 = Sd * x308
x310 = Gd * Gp
x311 = Fypp * x246
x312 = Fyph + x311
x313 = Bh + x270
x314 = Gh + x304
x315 = x11 ^ 2
x316 = x246 * x260
x317 = x252 * x44
x318 = x315 * x317
x319 = x246 * x318
x320 = x260 * x317
x321 = x298 * x300
x322 = x285 * x75
x323 = x155 * x266
x324 = A * x301
x325 = A * x243
x326 = x10 * x325
x327 = Gp * x298
x328 = x3 * x30
x329 = x275 * x328
x330 = Gp * x300
x331 = A * x264
x332 = Fyp * x246
x333 = Fyh + x332 + xi_yy
x334 = 8 * Gd
x335 = x11 * x334
x336 = x278 * x328
x337 = App * x246
x338 = x251 * x277
x339 = x259 * x292
x340 = 4 * Bd * Bp
x341 = x11 ^ 3
x342 = x11 * x90
x343 = x292 * x342
x344 = exp(4 * B)
x345 = x277 ^ 4
x346 = x344 * x345
x347 = 8 * x292
x348 = x19 * x260
x349 = x281 * x348
x350 = x277 ^ 3
x351 = x3 * x350
x352 = 2 * x351
x353 = Fxh + x298 + xi_xy
x354 = Fyt + x300 + xi_xy
x355 = Sph + Spp * x246
x356 = Sh + x268
x357 = x260 * x265
x358 = x356 * x357
x359 = x356 ^ 2
x360 = A * x0
x361 = x360 * x90
x362 = x313 ^ 2
x363 = Fxph + Fxpp * x246
x364 = x314 ^ 2
x365 = Fxpp * x277 + Fxpt
x366 = Fypp * x277 + Fypt
x367 = Gt + x294
x368 = Bp * x277
x369 = Bt + x368
x370 = 2 * x260
x371 = Sdh + Sdp * x246
x372 = Aph + x337
x373 = Bdh + Bdp * x246
x374 = A * x158
x375 = x285 * x44
x376 = x155 * x346
x377 = x287 * x346
x378 = x3 * x375
x379 = x246 * x277
x380 = 2 * A
x381 = Fxp * x277
x382 = Bp * x336
x383 = x278 * x3
x384 = x11 * x30
x385 = A * x384
x386 = x1 * x10 * x277
x387 = x277 * x308
x388 = 4 * x23
x389 = 6 * x285
x390 = x258 * x6
x391 = x30 * x94
x392 = x277 * x304
x393 = Fxt + x381 + xi_xx
x394 = x272 * x4
x395 = x393 * x4
x396 = Bd * x37
x397 = x22 * x356
x398 = x258 * x260
x399 = 8 * x398
x400 = Sd * x3
x401 = x31 * x400
x402 = kappa * x5
x403 = x279 * x292
x404 = x314 * x325
x405 = Fyp * x367
x406 = x0 * x246
x407 = 4 * x355
x408 = Gph + Gpp * x246
x409 = x279 * x341
x410 = 8 * Bd * Gp
x411 = x286 / 2
x412 = Bpp * x277
x413 = x260 * x44
x414 = x352 * x413
x415 = x246 * x61
x416 = 8 * x11
x417 = kappa * x416
x418 = Bd * Gp
x419 = x11 * x399 * x7
x420 = cosh(x284)
x421 = x259 * x420
x422 = x258 * x315
x423 = Bpp * x246
x424 = Bph + x423
x425 = St + x306
x426 = x425 ^ 2
x427 = x3 * x426
x428 = x369 ^ 2
x429 = x367 ^ 2
x430 = x23 * x287
x431 = x353 * x44
x432 = x252 * x431
x433 = x298 * x432
x434 = x398 * x44
x435 = 3 * x9
x436 = x252 * x354
x437 = x300 * x436 * x44
x438 = Sdp * x277 + Sdt
x439 = x246 * x308
x440 = App * x277 + Apt
x441 = Bdp * x277 + Bdt
x442 = x277 * x441
x443 = x246 * x251
x444 = A * x256
x445 = Bp * x354
x446 = Spp * x277 + Spt
x447 = x270 * x313
x448 = x315 * x90
x449 = x448 * x6
x450 = Bp * x265 * x294 * x449
x451 = x298 * x383
x452 = x258 * x449
x453 = x19 * x434
x454 = x260 * x267
x455 = x1 * x341
x456 = x159 * x287
x457 = x246 * x424
x458 = x277 * x408
x459 = Gpp * x277 + Gpt
x460 = x246 * x459
x461 = x277 * x369
x462 = x4 * x461
x463 = x277 * x367 * x4
x464 = Bd * x463
x465 = x344 * x350
x466 = x287 * x44
x467 = Bp * x466
x468 = x467 / 2
x469 = x0 * x265
x470 = Bp * x314
x471 = x243 * x46
x472 = x258 * x471
x473 = Fyp * x314
x474 = x260 * x278
x475 = x44 * x474
x476 = x475 * x9
x477 = x246 * x425
x478 = Gp * x265
x479 = x313 * x471
x480 = x277 * x384
x481 = 8 * x367
x482 = x6 * x90
x483 = 2 * x434
x484 = Bp * x333
x485 = x11 * x54
x486 = x127 * x243
x487 = x258 * x486
x488 = x10 * x404
x489 = x368 * x488
x490 = x10 * x360
x491 = x270 * x367
x492 = x304 * x367
x493 = x11 * x207
x494 = x280 * x44
x495 = x278 * x9
x496 = x494 * x495
x497 = x383 * x448
x498 = 2 * x3 * x475
x499 = x115 * x14 * x255 * x258
x500 = x0 * x341
x501 = x11 * x23
x502 = x313 * x356
x503 = x255 * x314
x504 = x289 * x313
x505 = Bpt + x412
x506 = x344 * x474
x507 = S * x425
x508 = Fxp * x507
x509 = Fxp * x314
x510 = x314 * x351
x511 = x258 * x277
x512 = x255 * x285
x513 = Bd * x512
x514 = Bd * x314
x515 = x265 * x292 * x80
x516 = x158 * x277
x517 = x343 * x351
x518 = x80 * x90
x519 = x0 * x258 * x518
x520 = x277 * x406
x521 = x518 * x520
x522 = x11 * x285
x523 = x1 * x268
x524 = 6 * x9
x525 = x465 * x505
x526 = Gh + x102 + x108
x527 = 6 * Gp
x528 = x285 * x46
x529 = x383 * x471
x530 = x277 * x285
x531 = x258 * x9
x532 = St + x151 + x18
x533 = 16 * S
x534 = x356 * x9
x535 = x277 * x533 * x534
x536 = x246 * x315
x537 = x425 * x536
x538 = x258 * x314
x539 = 4 * x398
x540 = x172 * x255 * x383
x541 = x267 * x316
x542 = Bp * x343
x543 = x300 * x406
x544 = x275 * x30
x545 = x367 * x534
x546 = x314 * x367
x547 = x369 * x465
x548 = 2 * x44 * x506
x549 = Fxp * x369
x550 = x483 * x9
x551 = x246 * x383
x552 = x277 * x536
x553 = x207 * x552
x554 = Bpp * xi_x + Bpt + x222
x555 = x292 * x554
x556 = x19 * x265 * x46
x557 = x258 * x425
x558 = x252 * x315
x559 = x539 * x9
x560 = x278 * x356
x561 = 4 * x9
x562 = x474 * x561
x563 = x356 * x562
x564 = x406 * (x131 + x143)
x565 = x3 * x425
x566 = x260 * x425
x567 = x23 * x292
x568 = x255 * x567
x569 = x0 * x342
x570 = x278 * x344
x571 = x471 * x570
x572 = Fxp * x367
x573 = x472 * x9
x574 = x471 * x495
x575 = x3 * x393
x576 = 2 * x476
x577 = x3 * x474
x578 = x383 * x90
x579 = x315 * x354
x580 = x292 * x424
x581 = x244 * x44
x582 = x351 * x581
x583 = x557 * x561
x584 = Sh + x40 + x43
x585 = Bt + x166 + x171
x586 = S * x383
x587 = x285 * x368
x588 = x313 * x368
x589 = 6 * x127
x590 = x494 * x531
x591 = Bph + Bpp * xi_y + x101
x592 = x260 * x277
x593 = x353 * x533
x594 = x533 * x9
x595 = x11 * x560
x596 = x270 * x383
x597 = Bh + x100 + x96
x598 = Gt + x178 + x187
x599 = x11 * x369
x600 = x11 * x379
x601 = x1 * x600
x602 = 4 * x506
x603 = x369 * x507
x604 = x313 * x367
x605 = x367 * x570
x606 = x532 * x551
x607 = x368 * x597
x608 = x207 * x258
x609 = x368 * x598
x610 = x17 + x35
x611 = Spp * x258 + Ss + x355 * x610
x612 = Bs + x288 + x424 * x610
x613 = Gpp * x258 + Gs + x408 * x610
x614 = S * x611
x615 = 2 * Spt
x616 = Fx * (x110 + x615) + Sb + Spp * x216 + xi_x * (2 * x110 + x615)
x617 = 2 * Bpt
x618 = Bb + Fx * (x222 + x617) + x231 + xi_x * (2 * x222 + x617)
x619 = 2 * Gpt
x620 = Fx * (x176 + x619) + Gb + x232 + xi_x * (2 * x176 + x619)
x621 = S * x616
x622 = 3 * x110
x623 = 3 * Spp
x624 = Fy * x622 + Sc + x76 + x82 + xi_x * (Fy * x623 + Sph + x623 * xi_y) + xi_y * (Spt + x622)
x625 = 3 * x176
x626 = 3 * Gpp
x627 = Fy * x625 + Gc + x133 + x142 + xi_x * (Fy * x626 + Gph + x626 * xi_y) + xi_y * (Gpt + x625)
axx = -x5
ayy = -x8
axy = x10 * x12
bx = -x32 * (Bt * x14 + Fx * x19 * x23 - x13 + x15 - x16 * x17 + x17 * x21 - x18 * x19 - x20 * x22 + x20 * x24 + x20 * x27 + x29 * (-Sp + x23 + x28) + xi_y * (2 * Sp * x20 - x25))
by = x36 * (-x13 * x20 + x15 * x20 + x17 * x23 - x17 * x28 - 2 * x21 * x26 - x22 + x24 + x27 + x29 * (x16 - x21) + x33 + x34 + x35 * (Sp + x23 - x28))
cc = -x75 * (Bb * x153 + Bh ^ 2 * x45 - Bs * x45 + Bt ^ 2 * x153 + Bt * Fxp * x153 + Bt * Gt * x195 + Bt * St * x95 - Bt * x130 + Fx * x194 + 2 * Fx * x242 * xi_x - Fxh * x126 - Fxh * x139 + Fxp * Gt * x154 - Fxp * x130 - Fxph * x92 + Fxpt * x153 + Fxt * x193 + Fxt * x199 - Fy * x141 + Fy * x150 + Fy * x152 - Fy * x227 - Fy * x229 + Fy * x39 - Fy * x50 + Fy * x67 + Fyh * x59 + Fyh * x71 + Fyp * x72 + Fyph * x45 - Fypt * x92 - Fyt * x126 - Fyt * x139 + Gb * x154 - Gc * x93 + Gh ^ 2 * x45 - Gh * Gt * x125 - Gp * x93 * xi_xy - Gph * x135 * xi_x + Gpp * x195 * x214 - Gpt * x134 * x66 + Gs * x47 + Gt ^ 2 * x153 - Gt * x6 * x97 + Sb * x95 - Sc * x14 * x48 - Sh ^ 2 * x37 + Sh * x42 + Sh * x48 * x49 - Sp * x60 * x75 - Sph * x14 * x86 - Spt * x77 * xi_y + Ss * x38 - St ^ 2 * x94 + St * x144 * x145 + x0 * x206 * x207 + x0 * x90 * x91 - x100 * x102 * x147 - x100 * x169 - x100 * x170 - x100 * x72 - x100 * x97 - x101 * x56 * x64 + x102 * x104 + x102 * x106 - x102 * x174 - x102 * x184 - x102 * x185 - x102 * x202 - x102 * x203 + x104 * x108 + x106 * x108 + x107 * x57 + x107 * x70 + x107 * x74 - x108 * x174 - x108 * x184 - x108 * x185 - x108 * x202 - x108 * x203 - x109 * x112 - x109 * x177 - x112 * x56 + x113 * x81 + x113 * x87 - x114 * x116 - x114 * x117 - x116 * x122 - x117 * x122 + x118 * x40 + x118 * x43 - x119 * x40 - x119 * x43 - x120 * x121 - x120 * x186 - x121 * x123 - x123 * x186 - x124 * x70 - x124 * x74 - 32 * x128 * xi_xy - x13 * x84 - x131 * x132 + x131 * x196 - x132 * x143 - x133 * x135 - x135 * x142 + x136 * x137 - x136 * x138 + x136 * x168 + x136 * x173 - x136 * x175 + x137 * x146 - x138 * x146 - x140 * x92 - x141 * xi_y + x143 * x196 + x146 * x168 + x146 * x173 - x146 * x175 + x147 * x176 * x200 + x150 * xi_y + x152 * xi_y - x157 * x69 - x157 * x73 + 6 * x160 * x161 + 70 * x160 * x238 + x161 * x164 + x164 * x238 - x165 * x79 - x165 * x85 - x166 * x167 - x167 * x171 - x169 * x96 + x17 * x205 * xi_y - x17 * x218 - x17 * x221 - x17 * x237 - x17 * x239 - x170 * x96 + 96 * x171 * x178 * x212 - x177 * x56 - x178 * x179 - x178 * x180 + x178 * x224 + x178 * x225 + x178 * x228 - x179 * x187 - x18 * x200 * x89 - x180 * x187 + x181 * x182 + x181 * x183 + x182 * x188 + x183 * x188 + x187 * x224 + x187 * x225 + x187 * x228 + x189 * x69 + x189 * x73 + x190 * x69 + x190 * x73 + 32 * x191 * x222 * xi_x + x193 * xi_xx + x194 * xi_x + x195 * x232 + x197 * x198 + x198 * x200 + x199 * xi_xx + x205 * x69 + x205 * x73 + x208 * x209 + x208 * x211 + x208 * x213 + x209 * x210 + x210 * x211 + x210 * x213 + x214 * x241 + x214 * x242 + x215 * x230 - x215 * x98 + x216 * x241 + x216 * x242 - x216 * x94 * x98 - x218 * x35 - x221 * x35 - x227 * xi_y - x229 * xi_y + x231 * x44 * x94 + x233 * x234 + x233 * x240 + x234 * x235 + x235 * x240 - x237 * x35 - x239 * x35 - x24 * x41 + 8 * x33 * (Bp * S * x0 * x246 + Fx * Gp * S * x0 * x6 + Gp * S * x0 * x6 * xi_x + Gt * S * x0 * x6 - Sh * x243 - x0 * x22 - x0 * x34 - x115 * x35 - x13 * x244 - x245 * x52 - x245 * x58) + x39 * xi_y - x40 * x41 + x40 * x42 - x40 * x51 + x40 * x83 - x40 * x89 * xi_y - x41 * x43 + x42 * x43 - x43 * x51 + x43 * x83 - x50 * xi_y + x52 * x53 - x52 * x84 + x53 * x58 + x55 * x57 + x55 * x70 + x55 * x74 - x58 * x84 + x59 * xi_yy - x60 * x62 - x63 * x65 + x63 * x68 - x65 * x66 + x66 * x68 + x67 * xi_y - x69 * x99 - 48 * x7 + x71 * xi_yy - x72 * x96 - x73 * x99 - x76 * x77 - x77 * x82 + x78 * x81 + x78 * x87 - x96 * x97) / 2
SS =  (4 * A * Bp * Fxp * x0 * x2 * x277 * x3 + A * Bp * Fxp * x11 * x2 * x277 * x285 * x3 + A * Bp * Fyp * x0 * x2 * x277 * x285 * x9 + 2 * A * Bp * Fyp * x2 * x246 * x292 * x6 + 4 * A * Bp * Gp * x11 * x2 * x258 * x260 * x6 + 2 * A * Bp * Gp * x11 * x2 * x278 * x3 * x420 + 2 * A * Bp * Gp * x11 * x2 * x278 * x3 + 4 * A * Bp * Sp * x0 * x278 * x3 * x30 * x315 + 8 * A * Bp * Sp * x11 * x246 * x277 * x30 * x9 + 4 * A * Bp * Sp * x258 * x292 * x30 * x6 + A * Bp * x0 * x2 * x246 * x285 * x369 * x9 + 2 * A * Bp * x0 * x2 * x246 * x367 * x420 * x9 + 2 * A * Bp * x0 * x2 * x277 * x285 * x3 * x367 + A * Bp * x0 * x2 * x277 * x285 * x313 * x9 + 4 * A * Bp * x0 * x2 * x277 * x3 * x369 + 6 * A * Bp * x0 * x2 * x3 * x393 + 2 * A * Bp * x0 * x2 * x333 * x6 + 2 * A * Bp * x0 * x246 * x285 * x30 * x425 * x9 + 4 * A * Bp * x0 * x246 * x30 * x356 * x6 + 2 * A * Bp * x11 * x2 * x246 * x314 * x420 * x6 + 4 * A * Bp * x11 * x2 * x246 * x314 * x6 + 2 * A * Bp * x11 * x2 * x277 * x285 * x314 * x9 + 4 * A * Bp * x11 * x2 * x277 * x3 * x367 + 2 * A * Bp * x11 * x246 * x285 * x30 * x356 * x6 + 4 * A * Bp * x11 * x260 * x277 * x30 * x356 * x9 + 2 * A * Bp * x2 * x246 * x292 * x313 * x6 + 2 * A * Bp * x2 * x277 * x292 * x3 * x369 + 4 * A * Bp * x277 * x292 * x3 * x30 * x425 + 2 * A * Bpp * x0 * x2 * x278 * x3 + A * Bpp * x11 * x2 * x278 * x285 * x3 + 2 * A * Bpp * x2 * x258 * x292 * x6 + 2 * A * Fxp * Gp * x11 * x2 * x277 * x3 + 4 * A * Fxp * Sp * x11 * x246 * x30 * x9 + 2 * A * Fxp * x0 * x2 * x3 * x369 + 2 * A * Fxp * x11 * x2 * x3 * x367 + 2 * A * Fyp * Gp * x11 * x2 * x246 * x6 + 4 * A * Fyp * Sp * x11 * x277 * x30 * x9 + 2 * A * Fyp * x11 * x2 * x314 * x6 + 8 * A * Gp * Sp * x0 * x246 * x277 * x30 * x9 + 2 * A * Gp * x0 * x2 * x246 * x314 * x6 + 2 * A * Gp * x0 * x2 * x277 * x3 * x367 + 2 * A * Gp * x11 * x2 * x277 * x3 * x369 + 2 * A * Gp * x11 * x2 * x3 * x393 + 2 * A * Gp * x11 * x2 * x333 * x6 - A * Sd * Sp * x250 + 4 * A * Sp * x0 * x246 * x356 * x44 * x6 + 4 * A * Sp * x0 * x277 * x3 * x425 * x44 + 2 * A * x0 * x158 * x2 * x258 * x6 + 2 * A * x0 * x158 * x2 * x278 * x3 + A * x0 * x162 * x2 * x258 * x6 + A * x0 * x162 * x2 * x278 * x3 + A * x0 * x2 * x206 * x3 + A * x0 * x2 * x246 * x285 * x505 * x9 + 2 * A * x0 * x2 * x3 * x365 + 2 * A * x0 * x2 * x3 * x428 + 2 * A * x0 * x2 * x3 * x429 + 2 * A * x0 * x2 * x3 * x618 + 2 * A * x0 * x2 * x312 * x6 + 2 * A * x0 * x2 * x313 * x367 * x9 + 2 * A * x0 * x2 * x362 * x6 + 2 * A * x0 * x2 * x364 * x6 + A * x0 * x2 * x6 * x91 + 4 * A * x0 * x3 * x30 * x369 * x425 + 4 * A * x0 * x3 * x30 * x616 + 4 * A * x0 * x30 * x6 * x611 + A * x11 * x158 * x2 * x258 * x285 * x6 + A * x11 * x158 * x2 * x278 * x285 * x3 + A * x11 * x2 * x246 * x285 * x424 * x6 + 2 * A * x11 * x2 * x246 * x408 * x6 + 2 * A * x11 * x2 * x260 * x277 * x424 * x9 + 2 * A * x11 * x2 * x277 * x3 * x459 + 4 * A * x11 * x2 * x3 * x367 * x369 + 2 * A * x11 * x2 * x3 * x620 + 2 * A * x11 * x2 * x6 * x613 + 4 * A * x11 * x246 * x30 * x446 * x9 + 4 * A * x11 * x277 * x30 * x355 * x9 + 4 * A * x11 * x3 * x30 * x367 * x425 + 4 * A * x11 * x30 * x314 * x356 * x6 + 8 * A * x11 * x356 * x425 * x44 * x9 - A * x11 * x4 * x505 * x530 - A * x155 * x255 * x260 * x279 + 2 * A * x2 * x277 * x292 * x3 * x505 - 12 * A * x247 - A * x268 * x425 * x493 - A * x293 * x356 * x61 - A * x295 * x627 - A * x306 * x356 * x493 - A * x368 * x369 * x4 * x522 - A * x406 * x407 * x61 + 4 * Ab * x0 * x2 * x3 - Ac * x251 + 4 * Ah * Bp * x0 * x2 * x246 * x315 * x6 + 4 * Ah * Fyp * x0 * x2 * x6 + 4 * Ah * Gp * x11 * x2 * x246 * x6 + 8 * Ah * Sp * x11 * x277 * x30 * x9 + 4 * Ah * x11 * x2 * x314 * x6 - Ah * x254 * x255 - Ah * x268 * x36 - Ah * x273 * x313 - Ah * x293 * x7 + 4 * Ap * Bp * x0 * x2 * x278 * x3 - Ap * Bp * x297 + 4 * Ap * Fxp * x0 * x2 * x277 * x3 + 4 * Ap * Fyp * x0 * x2 * x246 * x6 + 4 * Ap * Gp * x11 * x2 * x258 * x6 + 4 * Ap * Gp * x11 * x2 * x278 * x3 + 8 * Ap * Sd * x248 * x9 + 8 * Ap * Sp * x11 * x246 * x277 * x30 * x9 + 4 * Ap * x11 * x2 * x353 * x9 + 4 * Ap * x11 * x2 * x354 * x9 - Ap * x273 * x333 - Ap * x304 * x386 - Ap * x393 * x394 + 4 * App * x0 * x2 * x258 * x6 - App * x272 * x279 + 4 * As * x0 * x2 * x6 + 4 * At * Bp * x2 * x277 * x292 * x3 + 4 * At * Fxp * x0 * x2 * x3 - At * Fyp * x256 + 4 * At * Gp * x11 * x2 * x277 * x3 + 8 * At * Sp * x11 * x246 * x30 * x9 + 4 * At * x0 * x2 * x3 * x369 + 4 * At * x11 * x2 * x3 * x367 - At * x3 * x306 * x31 - At * x315 * x368 * x394 - Bd ^ 2 * x249 * x260 + 8 * Bd * Bp * x11 * x2 * x246 * x260 * x277 * x9 + 8 * Bd * Gp * x0 * x2 * x258 * x285 * x6 + 8 * Bd * Gp * x11 * x2 * x260 * x278 * x3 + 8 * Bd * Gp * x11 * x2 * x278 * x3 * x420 + 8 * Bd * Gp * x2 * x258 * x341 * x6 + 16 * Bd * Sd * x248 * x9 + 8 * Bd * Sp * x0 * x278 * x3 * x30 * x315 + 16 * Bd * Sp * x0 * x278 * x3 * x30 + 8 * Bd * Sp * x258 * x292 * x30 * x6 - Bd * Sp * x336 * x347 - Bd * Sp * x36 * x422 + 4 * Bd * x0 * x2 * x246 * x285 * x369 * x9 + 8 * Bd * x0 * x2 * x246 * x367 * x420 * x9 + 8 * Bd * x0 * x2 * x246 * x367 * x9 + 8 * Bd * x0 * x2 * x277 * x285 * x3 * x367 + 4 * Bd * x0 * x2 * x277 * x285 * x313 * x9 + 8 * Bd * x0 * x2 * x277 * x314 * x9 + 8 * Bd * x0 * x246 * x285 * x30 * x425 * x9 + 8 * Bd * x11 * x2 * x246 * x314 * x420 * x6 + 8 * Bd * x11 * x2 * x277 * x285 * x314 * x9 + 8 * Bd * x11 * x2 * x353 * x9 + 8 * Bd * x11 * x2 * x354 * x9 + 8 * Bd * x11 * x246 * x285 * x30 * x356 * x6 + 16 * Bd * x11 * x260 * x277 * x30 * x356 * x9 + 8 * Bd * x2 * x246 * x292 * x313 * x6 + 8 * Bd * x2 * x277 * x292 * x3 * x369 - Bd * x251 * x316 * x369 - Bd * x260 * x313 * x338 - Bd * x268 * x387 + 16 * Bd * x277 * x292 * x3 * x30 * x425 - Bd * x285 * x367 * x443 - 8 * Bd * x285 * x480 * x565 - 16 * Bd * x292 * x356 * x415 - Bd * x32 * x356 * x530 - Bd * x439 * x566 + 4 * Bp * Fxp * x0 * x246 * x278 * x3 * x341 * x44 + 2 * Bp * Fxp * x258 * x262 * x277 * x44 * x9 + 2 * Bp * Fxp * x262 * x344 * x350 * x44 + Bp * Fyp * x0 * x11 * x246 * x278 * x285 * x44 * x9 + 4 * Bp * Fyp * x11 * x258 * x277 * x292 * x44 * x6 + Bp * Fyp * x265 * x287 * x44 / 2 + 8 * Bp * Gd * x11 * x2 * x258 * x6 + 8 * Bp * Gd * x11 * x2 * x260 * x278 * x3 + 8 * Bp * Gd * x2 * x258 * x341 * x6 - Bp * Gd * x419 + 4 * Bp * Gp * x11 * x292 * x344 * x345 * x44 + 4 * Bp * Gp * x246 * x3 * x315 * x350 * x420 * x44 + 4 * Bp * Gp * x246 * x3 * x315 * x350 * x44 + Bp * Gp * x261 * x285 * x44 + Bp * Gp * x261 * x286 * x44 / 2 + 2 * Bp * Gp * x265 * x277 * x287 * x44 * x6 + 8 * Bp * S * Sp * x0 * x246 * x3 * x341 * x350 + 8 * Bp * S * Sp * x11 * x265 * x277 * x292 * x6 + 4 * Bp * S * Sp * x246 * x285 * x3 * x350 + 4 * Bp * S * Sp * x260 * x261 + Bp * S * Sp * x261 * x287 + 4 * Bp * S * Sp * x262 * x344 * x345 + 8 * Bp * S * x0 * x11 * x258 * x277 * x356 * x6 + 8 * Bp * S * x0 * x258 * x277 * x341 * x584 * x6 + 12 * Bp * S * x11 * x246 * x278 * x292 * x3 * x425 + 4 * Bp * S * x11 * x258 * x277 * x292 * x584 * x6 + 4 * Bp * S * x11 * x265 * x292 * x425 * x6 + 4 * Bp * S * x11 * x292 * x3 * x350 * x584 + 4 * Bp * S * x246 * x262 * x278 * x356 * x9 + 4 * Bp * S * x258 * x260 * x277 * x425 * x9 + Bp * S * x258 * x277 * x287 * x425 * x9 + 4 * Bp * S * x260 * x344 * x350 * x425 + 4 * Bp * S * x262 * x265 * x356 + Bp * S * x287 * x344 * x350 * x425 + 8 * Bp * Sd * x0 * x258 * x30 * x315 * x6 - Bp * Sd * x248 * x253 + 8 * Bp * Sd * x278 * x292 * x3 * x30 + 4 * Bp * x0 * x11 * x246 * x277 * x333 * x44 * x6 + Bp * x0 * x11 * x246 * x278 * x285 * x313 * x44 * x9 + Bp * x0 * x11 * x258 * x277 * x285 * x369 * x44 * x9 + 4 * Bp * x0 * x11 * x258 * x353 * x44 * x6 + 3 * Bp * x0 * x11 * x265 * x285 * x313 * x44 + 4 * Bp * x0 * x11 * x278 * x3 * x353 * x44 + 3 * Bp * x0 * x11 * x285 * x344 * x350 * x369 * x44 + 4 * Bp * x0 * x11 * x344 * x350 * x367 * x44 + 8 * Bp * x0 * x246 * x278 * x341 * x44 * x526 * x9 + 6 * Bp * x11 * x246 * x278 * x292 * x3 * x369 * x44 + 4 * Bp * x11 * x246 * x278 * x292 * x44 * x526 * x9 + 6 * Bp * x11 * x258 * x277 * x292 * x313 * x44 * x6 + 2 * Bp * x11 * x265 * x292 * x369 * x44 * x6 + 4 * Bp * x11 * x265 * x292 * x44 * x526 + 2 * Bp * x11 * x292 * x3 * x313 * x350 * x44 + 8 * Bp * x246 * x260 * x277 * x354 * x44 * x9 + 2 * Bp * x246 * x260 * x278 * x3 * x367 * x44 + 4 * Bp * x246 * x260 * x278 * x313 * x44 * x9 + 3 * Bp * x246 * x278 * x287 * x3 * x367 * x44 + 2 * Bp * x258 * x260 * x277 * x314 * x420 * x44 * x6 + 4 * Bp * x258 * x260 * x277 * x369 * x44 * x9 + 2 * Bp * x258 * x260 * x393 * x44 * x9 + (3 // 2) * Bp * x258 * x277 * x286 * x367 * x44 * x9 + 8 * Bp * x258 * x277 * x314 * x315 * x44 * x6 + 4 * Bp * x258 * x277 * x315 * x420 * x44 * x526 * x6 + 4 * Bp * x260 * x265 * x313 * x44 + 2 * Bp * x260 * x265 * x367 * x44 * x6 + 2 * Bp * x260 * x278 * x344 * x393 * x44 + 2 * Bp * x260 * x3 * x314 * x350 * x420 * x44 + 4 * Bp * x260 * x344 * x350 * x369 * x44 + Bp * x265 * x287 * x367 * x44 * x6 + 4 * Bp * x265 * x315 * x420 * x44 * x598 * x6 - Bp * x265 * x367 * x420 * x449 - Bp * x278 * x315 * x401 - Bp * x279 * x335 + Bp * x286 * x344 * x350 * x367 * x44 / 2 - Bp * x292 * x380 * x381 * x4 - Bp * x292 * x556 * x585 - Bp * x292 * x582 * x597 + 4 * Bp * x3 * x314 * x315 * x350 * x420 * x44 - Bp * x300 * x390 * x500 * x90 - Bp * x300 * x454 - Bp * x334 * x409 - Bp * x351 * x420 * x448 * x526 - Bp * x353 * x444 - Bp * x494 * x547 + 4 * Bpp * x0 * x246 * x3 * x341 * x350 * x44 + 4 * Bpp * x11 * x265 * x277 * x292 * x44 * x6 + 2 * Bpp * x246 * x285 * x3 * x350 * x44 + 2 * Bpp * x260 * x261 * x44 + Bpp * x261 * x287 * x44 / 2 + 2 * Bpp * x262 * x344 * x345 * x44 - Bpp * x380 * x403 + 8 * Fxc * x246 * x315 * x44 * x9 - Fxc * x316 * x317 + 8 * Fxp * Fyp * x246 * x260 * x277 * x44 * x9 + 3 * Fxp * Gp * x258 * x277 * x285 * x44 * x9 + Fxp * Gp * x285 * x344 * x350 * x44 - Fxp * Gp * x349 + 4 * Fxp * S * Sp * x258 * x260 * x277 * x9 + 8 * Fxp * S * Sp * x258 * x277 * x315 * x9 + 4 * Fxp * S * Sp * x260 * x344 * x350 + 8 * Fxp * S * x0 * x11 * x258 * x356 * x6 + 8 * Fxp * S * x0 * x11 * x278 * x3 * x356 + 8 * Fxp * S * x246 * x277 * x285 * x3 * x425 + 16 * Fxp * Sd * x0 * x277 * x3 * x30 + 4 * Fxp * x0 * x11 * x246 * x277 * x3 * x369 * x44 + 8 * Fxp * x246 * x260 * x277 * x313 * x44 * x9 + 8 * Fxp * x246 * x260 * x353 * x44 * x9 + 4 * Fxp * x246 * x277 * x3 * x315 * x367 * x44 + 8 * Fxp * x246 * x354 * x44 * x9 + 4 * Fxp * x258 * x314 * x315 * x44 * x6 + 6 * Fxp * x260 * x278 * x3 * x314 * x44 + 8 * Fxp * x277 * x315 * x333 * x44 * x9 - Fxp * x306 * x329 - Fxp * x313 * x383 * x569 - Fxp * x313 * x519 - Fxp * x465 * x468 + 8 * Fxs * x260 * x277 * x44 * x9 - Fxs * x277 * x318 + 8 * Fyb * x246 * x260 * x44 * x9 - Fyb * x319 + 8 * Fyc * x277 * x315 * x44 * x9 - Fyc * x277 * x320 + 3 * Fyp * Gp * x246 * x278 * x285 * x44 * x9 + Fyp * Gp * x265 * x285 * x44 - Fyp * Gp * x414 + 4 * Fyp * S * Sp * x246 * x260 * x278 * x9 + 8 * Fyp * S * Sp * x246 * x278 * x315 * x9 + 4 * Fyp * S * Sp * x260 * x265 + 8 * Fyp * S * x0 * x11 * x258 * x425 * x6 + 8 * Fyp * S * x0 * x11 * x278 * x3 * x425 + 8 * Fyp * S * x246 * x277 * x285 * x356 * x6 + 16 * Fyp * Sd * x0 * x246 * x30 * x6 + 4 * Fyp * x0 * x11 * x258 * x369 * x44 * x6 + 4 * Fyp * x0 * x11 * x278 * x3 * x369 * x44 + 4 * Fyp * x246 * x277 * x314 * x315 * x44 * x6 + 8 * Fyp * x246 * x315 * x393 * x44 * x9 + 2 * Fyp * x258 * x260 * x313 * x44 + 6 * Fyp * x258 * x260 * x367 * x44 * x6 + 8 * Fyp * x260 * x277 * x354 * x44 * x9 + 2 * Fyp * x260 * x278 * x313 * x44 * x9 - Fyp * x268 * x276 + 8 * Fyp * x277 * x353 * x44 * x9 + 4 * Fyp * x278 * x3 * x315 * x367 * x44 - Fyp * x283 - Fyp * x290 * x313 + 16 * Fypp * x246 * x278 * x315 * x44 * x9 - Gd ^ 2 * x249 + 8 * Gd * Gp * x11 * x2 * x246 * x277 * x9 + 8 * Gd * x0 * x2 * x277 * x313 * x9 + 8 * Gd * x0 * x2 * x353 * x9 + 8 * Gd * x0 * x2 * x354 * x9 - Gd * x1 * x10 * x246 * x369 + 16 * Gd * x11 * x246 * x30 * x356 * x6 + 16 * Gd * x11 * x277 * x3 * x30 * x425 - Gd * x277 * x356 * x391 - Gd * x391 * x477 + 8 * Gp * Sd * x11 * x258 * x30 * x6 + 8 * Gp * Sd * x11 * x278 * x3 * x30 + 4 * Gp * x0 * x11 * x246 * x277 * x353 * x44 * x9 + 4 * Gp * x0 * x11 * x246 * x277 * x354 * x44 * x9 - Gp * x10 * x325 * x353 - Gp * x11 * x303 * x331 + 2 * Gp * x246 * x260 * x278 * x314 * x44 * x9 + 4 * Gp * x246 * x278 * x314 * x315 * x44 * x9 + 2 * Gp * x258 * x260 * x277 * x367 * x44 * x9 + 2 * Gp * x258 * x260 * x353 * x44 * x6 + Gp * x258 * x277 * x285 * x369 * x44 * x9 + 4 * Gp * x258 * x277 * x313 * x315 * x44 * x6 + 4 * Gp * x258 * x277 * x315 * x367 * x44 * x9 + Gp * x258 * x285 * x333 * x44 + Gp * x258 * x285 * x393 * x44 * x9 + 4 * Gp * x258 * x315 * x354 * x44 * x6 + 2 * Gp * x260 * x265 * x314 * x44 + 2 * Gp * x260 * x278 * x3 * x354 * x44 + 2 * Gp * x260 * x344 * x350 * x367 * x44 - Gp * x264 * x383 * x385 - Gp * x270 * x352 * x466 + Gp * x278 * x285 * x333 * x44 * x9 + Gp * x278 * x285 * x344 * x393 * x44 + 4 * Gp * x278 * x3 * x315 * x353 * x44 + Gp * x285 * x344 * x350 * x369 * x44 - Gp * x314 * x351 * x471 - Gp * x326 * x354 - Gp * x353 * x452 - Gp * x578 * x579 + 8 * S * Sp * x0 * x11 * x258 * x277 * x313 * x6 + 4 * S * Sp * x0 * x11 * x265 * x314 + 4 * S * Sp * x0 * x11 * x344 * x350 * x367 + 8 * S * Sp * x246 * x277 * x315 * x353 * x9 + 8 * S * Sp * x246 * x277 * x315 * x354 * x9 + 6 * S * Sp * x246 * x278 * x285 * x314 * x9 + 4 * S * Sp * x258 * x260 * x277 * x369 * x9 + 4 * S * Sp * x258 * x260 * x333 + 4 * S * Sp * x258 * x260 * x393 * x9 + 6 * S * Sp * x258 * x277 * x285 * x367 * x9 + 4 * S * Sp * x260 * x278 * x333 * x9 + 4 * S * Sp * x260 * x278 * x344 * x393 + 4 * S * Sp * x260 * x344 * x350 * x369 - S * Sp * x389 * x451 + 8 * S * x0 * x11 * x246 * x277 * x3 * x369 * x425 + 8 * S * x0 * x11 * x246 * x277 * x3 * x616 + 8 * S * x0 * x11 * x246 * x277 * x6 * x611 + 8 * S * x0 * x11 * x258 * x6 * x624 + 8 * S * x0 * x11 * x278 * x3 * x624 - S * x11 * x277 * x523 * x575 + 4 * S * x246 * x260 * x278 * x355 * x9 + 16 * S * x246 * x260 * x353 * x425 * x9 + 8 * S * x246 * x277 * x3 * x315 * x367 * x425 + 8 * S * x246 * x277 * x314 * x315 * x356 * x6 + 8 * S * x246 * x278 * x315 * x355 * x9 + 16 * S * x246 * x315 * x354 * x425 * x9 + 4 * S * x258 * x260 * x277 * x446 * x9 + 4 * S * x258 * x260 * x313 * x356 + 4 * S * x258 * x260 * x314 * x425 * x6 + 4 * S * x258 * x260 * x356 * x367 * x6 + 8 * S * x258 * x277 * x315 * x446 * x9 - S * x260 * x264 * x510 + 4 * S * x260 * x265 * x355 + 16 * S * x260 * x277 * x354 * x356 * x9 + 4 * S * x260 * x278 * x3 * x314 * x425 + 4 * S * x260 * x278 * x3 * x356 * x367 + 4 * S * x260 * x278 * x313 * x356 * x9 + 4 * S * x260 * x344 * x350 * x446 - S * x264 * x313 * x357 - S * x268 * x313 * x562 - 4 * S * x268 * x367 * x577 + 16 * S * x277 * x315 * x353 * x356 * x9 - S * x389 * x446 * x551 - S * x545 * x601 + 32 * Sd ^ 2 * x2 * x9 + 8 * Sd * Sp * x0 * x258 * x44 * x6 + 8 * Sd * Sp * x0 * x278 * x3 * x44 + 8 * Sd * x0 * x246 * x30 * x313 * x6 + 8 * Sd * x0 * x246 * x30 * x367 * x9 + 8 * Sd * x0 * x277 * x30 * x314 * x9 + 16 * Sd * x11 * x246 * x425 * x44 * x9 + 16 * Sd * x11 * x277 * x356 * x44 * x9 + 8 * Sd * x11 * x30 * x353 * x9 + 8 * Sd * x11 * x30 * x354 * x9 - Sd * x195 * x268 * x277 - Sd * x246 * x356 * x93 - Sd * x258 * x347 * x62 - Sd * x314 * x415 * x416 - Sd * x333 * x36 - Sd * x391 * x392 + 4 * Sp * x0 * x11 * x265 * x425 * x6 + 4 * Sp * x0 * x11 * x3 * x350 * x356 - 2 * Sp * x15 * x307 + 6 * Sp * x246 * x278 * x285 * x3 * x425 + 6 * Sp * x258 * x277 * x285 * x356 * x6 + 8 * kappa * x0 * x2 * x246 * x313 * x6 + 8 * kappa * x0 * x2 * x246 * x367 * x9 + 8 * kappa * x0 * x2 * x277 * x314 * x9 + 8 * kappa * x11 * x2 * x353 * x9 + 8 * kappa * x11 * x2 * x354 * x9 - kappa * x333 * x8 + 4 * x0 * x11 * x246 * x277 * x3 * x428 * x44 + 4 * x0 * x11 * x246 * x277 * x3 * x429 * x44 + 4 * x0 * x11 * x246 * x277 * x3 * x44 * x618 + 4 * x0 * x11 * x246 * x277 * x313 * x367 * x44 * x9 + 4 * x0 * x11 * x246 * x277 * x362 * x44 * x6 + 4 * x0 * x11 * x246 * x277 * x364 * x44 * x6 + 6 * x0 * x11 * x246 * x278 * x408 * x44 * x9 + x0 * x11 * x258 * x277 * x285 * x44 * x505 * x9 + 8 * x0 * x11 * x258 * x277 * x424 * x44 * x6 + 6 * x0 * x11 * x258 * x277 * x44 * x459 * x9 + 4 * x0 * x11 * x258 * x313 * x314 * x44 + 4 * x0 * x11 * x258 * x314 * x367 * x44 * x6 + 4 * x0 * x11 * x278 * x3 * x314 * x367 * x44 + 4 * x0 * x11 * x278 * x313 * x314 * x44 * x9 + 3 * x0 * x11 * x285 * x344 * x350 * x44 * x505 + 4 * x0 * x158 * x246 * x3 * x341 * x350 * x44 + 4 * x0 * x158 * x265 * x277 * x341 * x44 * x6 + 4 * x0 * x2 * x246 * x285 * x441 * x9 + 8 * x0 * x2 * x246 * x372 * x6 + 8 * x0 * x2 * x246 * x373 * x6 + 8 * x0 * x2 * x277 * x3 * x440 + 16 * x0 * x246 * x30 * x371 * x6 + 4 * x0 * x258 * x277 * x341 * x44 * x591 * x6 - x0 * x277 * x298 * x314 * x493 + 16 * x0 * x277 * x3 * x30 * x438 - x0 * x294 * x538 * x589 - x0 * x367 * x599 * x608 - 6 * x0 * x383 * x46 * x492 - x1 * x356 * x557 * x80 - x1 * x383 * x477 * x501 - x1 * x565 * x595 - x11 * x129 * x258 * x588 + 24 * x11 * x2 * x246 * x277 * x9 + 4 * x11 * x2 * x246 * x285 * x373 * x6 + 8 * x11 * x2 * x260 * x277 * x373 * x9 + 6 * x11 * x246 * x278 * x292 * x3 * x44 * x505 + 2 * x11 * x258 * x277 * x292 * x44 * x591 * x6 + 2 * x11 * x265 * x292 * x44 * x505 * x6 + 2 * x11 * x292 * x3 * x350 * x44 * x591 - x11 * x292 * x608 * x609 - x11 * x369 * x523 * x586 - x11 * x410 * x421 - x11 * x427 * x564 - x11 * x45 * x505 * x551 - x111 * x265 * x272 * x446 - x111 * x277 * x333 * x523 - x113 * x285 * x352 - x12 * x464 - x14 * x264 * x357 * x367 - x14 * x306 * x314 * x539 - 8 * x14 * x306 * x315 * x538 - x14 * x355 * x389 * x511 - x14 * x502 * x601 - x140 * x258 * x322 - x140 * x267 - x140 * x278 * x378 - x148 * x265 * x306 * x455 - x153 * x258 * x341 * x609 - x153 * x600 * x627 - x155 * x261 * x343 - x158 * x207 * x278 * x398 + 2 * x158 * x246 * x285 * x3 * x350 * x44 - x158 * x246 * x517 + 4 * x158 * x258 * x262 * x278 * x44 * x9 + 2 * x158 * x260 * x261 * x315 * x44 + 2 * x158 * x260 * x315 * x344 * x345 * x44 + 2 * x158 * x261 * x262 * x44 + 2 * x158 * x262 * x344 * x345 * x44 + 2 * x158 * x265 * x277 * x285 * x44 * x6 - x159 * x260 * x269 - x159 * x261 * x287 - x159 * x346 * x370 - x162 * x267 * x379 - x19 * x258 * x292 * x46 * x607 - x19 * x281 * x285 * x412 + 8 * x2 * x277 * x292 * x3 * x441 + 4 * x206 * x258 * x315 * x44 * x9 + 4 * x206 * x258 * x44 * x9 + x206 * x260 * x278 * x344 * x44 - x206 * x378 * x379 - x206 * x434 * x435 - x207 * x315 * x321 - x207 * x457 * x474 - 16 * x22 * x277 * x537 * x9 - x22 * x306 * x389 * x390 - x23 * x246 * x563 - x23 * x260 * x264 * x346 - x23 * x455 * x606 - x230 * x262 * x269 - x230 * x346 * x370 - x230 * x377 / 2 - x244 * x289 * x304 * x313 + 8 * x246 * x260 * x277 * x363 * x44 * x9 + 8 * x246 * x260 * x277 * x366 * x44 * x9 + 2 * x246 * x262 * x278 * x424 * x44 * x9 + 8 * x246 * x277 * x3 * x315 * x367 * x369 * x44 + 4 * x246 * x277 * x3 * x315 * x44 * x620 + 4 * x246 * x277 * x315 * x44 * x6 * x613 - x246 * x285 * x514 * x8 - x246 * x313 * x513 * x7 - x246 * x314 * x417 * x7 + x246 * x315 * x356 * x425 * x9 * (x79 + x85) - x246 * x347 * x373 * x7 - x246 * x356 * x430 * x495 - x251 * x316 * x441 - x252 * x397 * x474 - x252 * x398 * x508 - x253 * x384 * x624 - x254 * x404 - x255 * x425 * x54 * x605 - x257 * x259 - x257 * x279 + 4 * x258 * x260 * x277 * x44 * x505 * x9 + 2 * x258 * x260 * x312 * x44 + 2 * x258 * x260 * x314 * x369 * x44 * x6 + 4 * x258 * x260 * x359 + 4 * x258 * x260 * x426 * x9 + 4 * x258 * x260 * x44 * x6 * x627 + 2 * x258 * x260 * x44 * x612 + x258 * x260 * x44 * x91 + 8 * x258 * x315 * x365 * x44 * x9 - x258 * x324 * x62 - x258 * x356 * x503 * x54 - x258 * x456 * x495 - x258 * x482 * x500 * x607 - x259 * x323 - x260 * x264 * x425 * x465 + 2 * x260 * x278 * x3 * x314 * x369 * x44 + 4 * x260 * x278 * x3 * x44 * x627 + 2 * x260 * x278 * x344 * x365 * x44 + 4 * x260 * x278 * x344 * x426 + 4 * x260 * x278 * x359 * x9 + 2 * x260 * x278 * x44 * x612 * x9 + 4 * x260 * x344 * x350 * x44 * x505 - x261 * x263 * x264 + 2 * x262 * x265 * x424 * x44 - x263 * x277 * x583 - 4 * x263 * x425 * x465 - x264 * x358 + x265 * x285 * x408 * x44 + x265 * x287 * x424 * x44 - x265 * x356 * x430 - x266 * x328 * x425 * x587 - x266 * x367 * x368 * x4 * x420 - x267 * x285 * x491 - x267 * x294 * x314 - x267 * x363 - x267 * x366 - x267 * x492 - x268 * x315 * x481 * x586 - x268 * x347 * x351 * x501 - x268 * x558 * x560 - x268 * x563 - x270 * x277 * x569 * x575 - x270 * x285 * x290 * x314 - 3 // 2 * x270 * x286 * x314 * x44 * x495 - x270 * x369 * x454 - x270 * x383 * x45 * x599 - x270 * x385 * x561 * x566 - x270 * x420 * x497 * x598 - x270 * x432 * x592 - x270 * x500 * x578 * x585 - x270 * x504 * x522 - x271 * x274 - x271 * x290 * x315 - x271 * x496 - x272 * x279 * x310 - x274 * x447 - x276 * x502 - x277 * x319 * x546 - x277 * x329 * x446 - x277 * x37 * x400 * x425 * x44 - x277 * x505 * x590 - x278 * x285 * x418 * x5 + 8 * x278 * x312 * x315 * x44 * x9 + 4 * x278 * x315 * x44 * x9 * x91 - x278 * x413 * x435 * x91 + 4 * x278 * x44 * x9 * x91 - x282 * x287 * x313 - x282 * x314 * x411 - x282 * x348 * x367 * x420 - x283 * x313 + x285 * x344 * x350 * x44 * x459 - 3 * x287 * x368 * x538 * x75 - 2 * x289 * x292 * x457 - x290 * x612 - x291 * x315 - x291 - x292 * x331 * x382 - x292 * x581 * x585 * x596 - x294 * x296 - x295 * x373 * x530 - x296 * x367 - x297 * x310 - x298 * x299 - x298 * x309 - x298 * x315 * x535 - x299 * x300 - x3 * x342 * x365 * x520 - x3 * x392 * x393 * x448 - x30 * x325 * x534 * x587 - x300 * x309 - x300 * x316 * x317 * x369 - x302 * x303 - x302 * x336 - x304 * x305 - x304 * x369 * x497 - x304 * x479 * x495 - x305 * x314 - x306 * x307 * x388 * x6 - x306 * x425 * x559 - x306 * x557 * x558 - 16 * x311 * x476 - x312 * x475 * x524 - x312 * x521 - 8 * x313 * x314 * x552 * x75 - x313 * x518 * x543 - x314 * x368 * x453 - x314 * x406 * x461 * x493 - x314 * x425 * x544 * x9 - x314 * x507 * x601 * x9 - x314 * x54 * x561 * x595 - x315 * x433 - x315 * x437 - x315 * x44 * x481 * x596 - x316 * x354 * x425 * x594 - x317 * x321 - x320 * x332 * x393 - x320 * x333 * x381 - x322 * x379 * x91 - x323 * x421 - x324 * x382 - x324 * x422 * x62 - x326 * x327 - x326 * x330 - x326 * x405 - x326 * x458 - x326 * x460 - x327 * x497 - x327 * x498 - x329 * x368 * x425 - x330 * x452 - x330 * x453 - x333 * x335 * x7 - x333 * x392 * x449 - x335 * x395 - x337 * x338 - x338 * x372 - x339 * x340 - x339 * x374 - x340 * x403 - x341 * x406 * x554 * x578 - x341 * x412 * x469 * x482 - x342 * x469 * x470 - x346 * x456 - x349 * x459 - x351 * x356 * x568 - x351 * x407 * x485 + 4 * x353 ^ 2 * x44 * x9 - x353 * x499 - x353 * x540 + 4 * x354 ^ 2 * x44 * x9 - x354 * x398 * x527 * x75 - x354 * x499 - x354 * x540 - 12 * x356 * x511 * x567 * x80 - x357 * x424 * x90 - x358 * x388 - x359 * x361 * x6 - x359 * x564 * x80 - x361 * x427 - x362 * x483 - x362 * x576 - x363 * x487 - x363 * x529 - x363 * x553 - x364 * x483 - x364 * x576 - x365 * x434 * x524 - x366 * x487 - x366 * x529 - x366 * x553 - x367 * x478 * x486 - x367 * x485 * x583 - x367 * x493 * x543 - x368 * x369 * x590 - x369 * x488 - x369 * x569 * x605 - x371 * x387 - x374 * x403 - x375 * x376 - x376 * x411 * x44 - x377 * x98 - x381 * x468 * x531 - x383 * x445 * x569 - x386 * x420 * x514 - x388 * x515 * x532 - x393 * x401 - x393 * x402 - x395 * x396 - x396 * x462 - x397 * x399 - x4 * x442 * x512 - x400 * x480 * x481 - x401 * x461 - x402 * x461 - x405 * x452 - x405 * x498 - x406 * x424 * x495 * x528 - x408 * x414 - x408 * x449 * x511 - x409 * x410 - x414 * x470 - x416 * x420 * x464 - x417 * x463 - x418 * x419 - x420 * x450 - x420 * x489 - x420 * x491 * x498 - x423 * x517 - 3 * x424 * x469 * x528 - x424 * x490 * x530 - x428 * x548 - x428 * x550 - x429 * x548 - x429 * x550 - x431 * x436 - x431 * x527 * x577 - x433 - x437 - x438 * x439 - x440 * x443 - x442 * x5 - x444 * x445 - x444 * x546 - x445 * x519 - x447 * x496 - x450 - x451 * x542 - x453 * x458 - x453 * x509 - x453 * x604 - x454 * x588 - x459 * x536 * x578 - x460 * x498 - x462 * x513 - x465 * x542 * x598 - x466 * x525 - x467 * x510 - x467 * x547 - x472 * x473 - x472 * x613 - x473 * x574 - x478 * x479 - x483 * x484 - x484 * x576 - x489 - 6 * x490 * x491 - x494 * x525 - x497 * x509 - x498 * x604 - x503 * x504 - x505 * x541 - 8 * x506 * x508 - x511 * x580 * x589 - x515 * x516 * x90 - x516 * x541 - x521 * x612 - x534 * x592 * x593 - x535 * x579 - x537 * x593 * x9 - x539 * x614 - x544 * x545 - x548 * x549 - x548 * x618 - x549 * x550 - x550 * x618 - x551 * x555 * x581 - x552 * x594 * x624 - x555 * x556 - x559 * x603 - x559 * x621 - x562 * x614 - x568 * x606 - x571 * x572 - x571 * x620 - x572 * x573 - x573 * x620 - x574 * x613 - x580 * x582 - x602 * x603 - x602 * x621)

    return axx, ayy, axy, bx, by, cc, SS
end
