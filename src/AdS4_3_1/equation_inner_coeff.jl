
#= tilde, hat, etc, definitions

We use these macros as shorthand notation. For instance

  @tilde_inner("B")

should expand to

  Bt = B_x -  (Fx * u + xi_x) * Bp

etc.

=#
macro tilde_inner(fname::String)
    ft    = Symbol(fname, "t")
    f_x   = Symbol(fname, "_x")
    fp    = Symbol(fname, "p")
    return esc( :($ft = $f_x - (Fx * u + xi_x) * $fp) )
end
macro hat_inner(fname::String)
    fh    = Symbol(fname, "h")
    f_y   = Symbol(fname, "_y")
    fp    = Symbol(fname, "p")
    return esc( :($fh = $f_y - (Fy * u + xi_y) * $fp) )
end

macro bar_inner(fname::String)
    fb    = Symbol(fname, "b")
    f_xx  = Symbol(fname, "_xx")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    return esc( :($fb = $f_xx + (Fx * u + xi_x) * ( -2*($fp_x) + (Fx * u + xi_x) * ($fpp) )) )
end

macro star_inner(fname::String)
    fs    = Symbol(fname, "s")
    f_yy  = Symbol(fname, "_yy")
    fpp   = Symbol(fname, "pp")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fs = $f_yy + (Fy * u + xi_y) * ( -2*($fp_y) + (Fy * u + xi_y)* ($fpp) )) )
end

macro cross_inner(fname::String)
    fc    = Symbol(fname, "c")
    f_xy  = Symbol(fname, "_xy")
    fpp   = Symbol(fname, "pp")
    fp_x  = Symbol(fname, "p_x")
    fp_y  = Symbol(fname, "p_y")
    return esc( :($fc = $f_xy  - (Fx * u + xi_x) * ($fp_y) -
                  (Fy * u + xi_y) * ( $fp_x -(Fx * u + xi_x) * ($fpp) ) ) )
end



# assuming
# (A d_uu + B d_u + C Id) f = -S

function S_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    ( u, xi, B, Bp, G, Gp) = vars

 
  
x0 = u ^ 6
x1 = u ^ 5
x2 = u ^ 4
x3 = 3 * u
x4 = (-B * x3 + Bp) ^ 2 * cosh(G * u ^ 3) ^ 2
ABCS[1] = 4 * x0
ABCS[2] = 24 * x1
ABCS[3] = x2 * (9 * G ^ 2 * x0 - 6 * G * Gp * x1 + Gp ^ 2 * x2 + x2 * x4 + 24)
ABCS[4] = - (x0 * (x4 + (-G * x3 + Gp) ^ 2) * (xi + 1 / u))

    nothing
end


# this is a coupled equation for Fx and Fy. the notation used is
#
# ( A11 d_uu Fx + A12 d_uu Fy + B d_u Fx + B2 d_u Fy + C11 Fx + C12 Fy ) = -S1
# ( A21 d_uu Fx + A22 d_uu Fy + B21 d_u Fx + B22 d_u Fy + C21 Fx + C22 Fy ) = -S2

function Fxy_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        u, xi, xi_x, xi_y,
        B     ,        G      ,       S      ,
        Bp    ,        Gp     ,       Sp     ,
        Bpp   ,        Gpp    ,       Spp    ,
        B_x   ,        G_x    ,       S_x    ,
        B_y   ,        G_y    ,       S_y    ,
        Bp_x  ,        Gp_x   ,       Sp_x   ,
        Bp_y  ,        Gp_y   ,       Sp_y   ,
       
    ) = vars
    
    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Sp")
   

    @hat_inner("Sp")
    

    @cross_inner("G")
    @cross_inner("S")

    u2 = u*u
    u3 = u*u2
    u4 = u2*u2
    u5 = u4*u
    u5 = u3*u2
    u6 = u3*u3
    u8 = u4*u4



    x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = u * xi
x4 = S * x0
x5 = x3 + x4 + 1
x6 = x5 ^ 2
x7 = 2 * x0 * x6
x8 = u ^ 2
x9 = Bp * x8
x10 = 3 * x1
x11 = -Bp
x12 = 3 * u
x13 = B * x12
x14 = x11 + x13
x15 = G * x0
x16 = 2 * x15
x17 = cosh(x16)
x18 = x14 * x17
x19 = x10 + x18 * x8 - x9
x20 = x2 * x8
x21 = x20 * x6
x22 = 2 * Gp
x23 = G * u
x24 = 6 * x23
x25 = -x14
x26 = sinh(x16)
x27 = x25 * x26
x28 = -x22 + x24 - x27
x29 = u ^ 4
x30 = x29 * x6
x31 = 8 * xi
x32 = 2 * u
x33 = Sp * u
x34 = S * x8
x35 = xi ^ 2
x36 = 4 * x35
x37 = Bp * x29
x38 = 3 * B
x39 = u ^ 5
x40 = S * x39
x41 = 6 * xi
x42 = G * Gp
x43 = x29 * x42
x44 = x4 * xi
x45 = Sp * x8
x46 = B ^ 2
x47 = 9 * x39
x48 = G ^ 2
x49 = 18 * x48
x50 = S ^ 2
x51 = x39 * x50
x52 = u ^ 7
x53 = Bp * x52
x54 = B * S
x55 = B * Bp * x41
x56 = u ^ 6
x57 = x56 * xi
x58 = 12 * xi
x59 = S * x42 * x52
x60 = x39 * x42
x61 = u ^ 8
x62 = 18 * x46
x63 = S * x61
x64 = 36 * x48
x65 = x42 * x63
x66 = u ^ 10 * x50
x67 = Bp * x38
x68 = x35 * x56
x69 = S * u ^ 9 * xi
x70 = 6 * x42
x71 = 9 * x46
x72 = u ^ 11 * x50
x73 = x35 * x52
x74 = x27 * x30
x75 = G * x74
x76 = x10 * x5
x77 = 8 * x4
x78 = 4 * x3 + x77 + 2
x79 = x76 + x78
x80 = 5 * x3 + 9 * x4 + 3
x81 = x76 + x80
x82 = -Bp * x79 + x13 * x81
x83 = u * x5
x84 = 6 * x5
x85 = x17 * x25
x86 = x15 * x84 * x85
x87 = -x76
x88 = x80 + x87
x89 = x78 + x87
x90 = x0 * x5
x91 = x32 * x6
x92 = x91 * xi_y
x93 = x6 * x8
x94 = 6 * x93
x95 = 4 * x83
x96 = x2 * xi_x
x97 = x95 * x96
x98 = -2 * x4 + x45 + 1
x99 = 4 * x98
x100 = x26 * x6
x101 = u * x100
x102 = 8 * x5
x103 = St * x20
x104 = x101 * xi_y
x105 = x102 * x45
x106 = 4 * x98 ^ 2
x107 = Bh * x26
x108 = -Gp
x109 = G * x12
x110 = x108 + x109
x111 = 2 * x30
x112 = x110 * x111
x113 = x9 * xi_y
x114 = 3 * x100
x115 = 4 * x45
x116 = x115 * x98
x117 = x112 * xi_y
x118 = 6 * u
x119 = Gpp + x118 * (x108 + 2 * x23)
x120 = cosh(x15) ^ 2
x121 = x120 * x91
x122 = Bpt * x2
x123 = x121 * x96
x124 = 6 * x21
x125 = Bt * x120
x126 = Spp - 4 * x33 + 6 * x34
x127 = x30 * xi_y
x128 = x25 ^ 2
x129 = x128 * x26
x130 = x127 * x22
x131 = Gh * x111
x132 = 6 * x120 * x6
x133 = x9 * x96
x134 = Bpp + x118 * (B * x32 + x11)
x135 = 4 * x120 * x14 * x5
x136 = St * x2 * x29
x137 = x20 * xi_x
x138 = x135 * x137
x139 = Sh * x29
x140 = x28 * x5
x141 = 2 * x140
x142 = 2 * x8
x143 = x140 * x142 * xi_y
x144 = x111 * x120
x145 = Bp * x144 * x25
x146 = Sp * x29
x147 = x146 * x96
x148 = x146 * xi_y
x149 = x2 * x30
x150 = 2 * x149
x151 = x128 * x144
x152 = -Gp + x109 - x27
x153 = Gt * x150
x154 = x30 * x96
x155 = x154 * x22
x156 = -x98
x157 = x112 * x96
x158 = x22 - x24 - x27
x159 = 4 * xi
x160 = 3 * x42
x161 = 9 * x48
x162 = x95 * xi_y
x163 = Sh * x8
x164 = x91 * x96
x165 = x120 * x92
x166 = Bh * x120
x167 = x101 * x96
x168 = Bt * x26
x169 = x135 * x8 * xi_y
x170 = x110 + x27
x171 = 2 * x158 * x5
x172 = x137 * x171
AA[1,1] = x2 * x7
AA[1,2] = 0
BB[1,1] = x21 * (x19 + 8)
BB[1,2] = -x28 * x30
CC[1,1] = x20 * (15 * B * x29 * x35 + 36 * B * x40 + 27 * B * x50 * x61 + 9 * B * x8 - Bp * x0 * x36 - Bp * x32 - Bp * x40 * x58 - 10 * S * x37 + u * x36 + 24 * x1 * xi + x17 * x82 * x83 - x31 * x45 + x31 - 12 * x33 + 48 * x34 - x37 * x38 + x39 * x49 - x39 * x55 - x41 * x9 - 6 * x43 + 32 * x44 + x46 * x47 + x49 * x72 + x49 * x73 - 8 * x50 * x53 + 12 * x51 - 6 * x53 * x54 + 42 * x54 * x57 - x55 * x63 + x57 * x62 + x57 * x64 - x58 * x60 - x58 * x65 - 12 * x59 + x62 * x63 + x62 * x69 + x63 * x64 + x64 * x69 - x66 * x67 - x66 * x70 - x67 * x68 - x68 * x70 + x71 * x72 + x71 * x73 - 6 * x75)
CC[1,2] = -x90 * (Gp * x1 * x84 - 4 * Gp * (2 * x3 + 4 * x4 + 1) + x24 * x88 + x26 * (-Bp * x89 + x13 * x88) - x86)
SS[1] = -(Bh * x112 + Bp * x117 - Bp * x127 * x27 + Bph * x101 + Bpp * x104 - Bpp * x123 - Gh * x94 - Gp * x94 * xi_y + Gph * x91 + Gpp * x92 - Spp * x97 - Spt * x2 * x95 + x102 * x103 + x103 * x99 - x104 * x134 + x105 * x96 - x106 * x96 + x107 * x14 * x30 - 3 * x107 * x93 - x113 * x114 + x116 * x96 + x117 * x14 + x117 * x85 - x119 * x92 - x121 * x122 + x123 * x134 + x124 * x125 + x125 * x14 * x150 + x126 * x97 + x127 * x129 + x130 * x85 - x131 * x18 + x132 * x133 + x135 * x136 + x135 * x147 + x138 * x156 + x138 - x139 * x141 - x141 * x148 - x143 * x156 - x143 - x145 * x96 + x151 * x96 + x152 * x153 + x152 * x155 + x152 * x157 + x96 * x99)
AA[2,1] = 0
AA[2,2] = x7
BB[2,1] = x149 * x158
BB[2,2] = -x93 * (x19 - 8)
CC[2,1] = x2 * x90 * (2 * Gp * x79 - x24 * x81 + x26 * x82 - x86)
CC[2,2] = x142 * (x120 * x83 * (Bp * x89 - x13 * x88) - x159 * x45 + x159 - x160 * x66 - x160 * x68 + x161 * x72 + x161 * x73 + x32 * x35 - 6 * x33 + 24 * x34 - x41 * x60 - x41 * x65 - 3 * x43 + 16 * x44 + x47 * x48 + x49 * x57 + x49 * x63 + x49 * x69 + 6 * x51 - 6 * x59 + 3 * x75)
SS[2] = -(-Bp * x157 - Bp * x74 * x96 + Bph * x121 + Bpp * x165 - Bpp * x167 - Bt * x110 * x150 - Gp * x124 * xi_x + Gpp * x164 + Gpt * x2 * x91 - Gt * x124 - Sph * x95 - Spp * x162 - x101 * x122 + x102 * x163 + x105 * xi_y - x106 * xi_y + x111 * x14 * x166 - x113 * x132 + x114 * x133 + x116 * xi_y + x117 * x170 - x119 * x164 + x126 * x162 + x129 * x154 + x130 * x170 + x131 * x170 - x134 * x165 + x134 * x167 - x135 * x139 - x135 * x148 + x136 * x171 + x14 * x149 * x168 - x14 * x157 - x145 * xi_y + x147 * x171 + x151 * xi_y + x153 * x18 + x155 * x18 - x156 * x169 + x156 * x172 + x157 * x18 + x163 * x99 - x166 * x94 + 3 * x168 * x21 - x169 + x172 + xi_y * (x115 - x77 + 4))


    nothing
end


function Sd_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        potential,  u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,        G      ,        S      ,    Fx     ,    Fy     ,
        Bp    ,        Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,        Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,        G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,        G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,        Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,        Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,        G_xx   ,        S_xx   ,
        B_yy  ,        G_yy   ,        S_yy   ,
                        G_xy   ,        S_xy  ,   St
    ) = vars

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Sp")
    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Sp")
    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")

x0 = u ^ 3
x1 = S * x0
x2 = u * xi + x1 + 1
x3 = B * x0
x4 = exp(x3)
x5 = 8 * x4
x6 = 8 * xi
x7 = u ^ 2
x8 = S * x7
x9 = 1 / u
x10 = x8 + x9 + xi
x11 = x10 ^ 3
x12 = x4 * x6
x13 = Sp * x7
x14 = -2 * x1 + x13 + 1
x15 = 4 * x14
x16 = x10 ^ 2
x17 = x16 * x4
x18 = x15 * x17
x19 = Fy * u
x20 = x19 + xi_y
x21 = Sh + Sp * x20
x22 = x21 * x7
x23 = x14 * x20
x24 = x22 - x23 + xi_y
x25 = G * x0
x26 = cosh(x25)
x27 = 4 * x26
x28 = 8 * xi_y
x29 = sinh(x25)
x30 = Fx * u
x31 = x30 + xi_x
x32 = Sp * x31 + St
x33 = -x14 * x31 + x32 * x7 + xi_x
x34 = x29 * x33
x35 = Fxh + Fxp * x20
x36 = Fyp * x31 + Fyt
x37 = Fx * x7
x38 = -Fxp * u + x37
x39 = -x38
x40 = Fy * x7
x41 = -Fyp * u + x40
x42 = -x41
x43 = 2 * x0
x44 = 2 * x7
x45 = Fy * x43 - Fyp * x44 + Fypp * u
x46 = Fx * x43 - Fxp * x44 + Fxpp * u
x47 = Gh + Gp * x20
x48 = x0 * x47
x49 = u ^ 4
x50 = 3 * x49
x51 = -G * x50 + Gp * x0
x52 = x20 * x51
x53 = x48 - x52
x54 = Gp * x31 + Gt
x55 = x0 * x54
x56 = 2 * x55
x57 = 2 * x30 + 2 * xi_x
x58 = 2 * x17
x59 = x20 ^ 2
x60 = 2 * x19 + 2 * xi_y
x61 = Spp * x20
x62 = Sph + x61
x63 = 4 * x7
x64 = x21 * x43 - x62 * x7
x65 = -x64
x66 = 4 * x20
x67 = Fyh + Fyp * x20
x68 = u * x67 - x20 * x42 + xi_yy
x69 = 6 * x49
x70 = S * x69 - 4 * Sp * x0 + Spp * x7
x71 = Bh + Bp * x20
x72 = -B * x50 + Bp * x0
x73 = x20 * x72
x74 = x0 * x71 - x73
x75 = 4 * x24
x76 = 4 * x33
x77 = x31 * x51
x78 = x55 - x77
x79 = Spp * x31 + Spt
x80 = x32 * x43 - x7 * x79
x81 = -x31 * x70 - x80
x82 = u * x35 + u * x36 - x20 * x39 - x31 * x42 + 2 * xi_xy
x83 = Gpp * x20
x84 = Gph + x83
x85 = -x0 * x84 + x47 * x50
x86 = -x85
x87 = Gpp * x31 + Gpt
x88 = 12 * u ^ 5
x89 = G * x88 - Gp * x69 + Gpp * x0
x90 = -x0 * x87 + x50 * x54
x91 = -x31 * x89 - x90
x92 = Bp * x31 + Bt
x93 = x0 * x92
x94 = -x31 * x72 + x93
x95 = -G * x69 + 2 * Gp * x0
x96 = -B * x69 + 2 * Bp * x0
x97 = 2 * u
x98 = Bph + Bpp * x20
x99 = -x0 * x98 + x50 * x71
x100 = B * x88 - Bp * x69 + Bpp * x0
x101 = x31 ^ 2
x102 = 4 * x31
x103 = Fxp * x31 + Fxt
x104 = u * x103 - x31 * x39 + xi_xx
x105 = Bpp * x31 + Bpt
x106 = -x0 * x105 + x50 * x92
ABCS[1] = 0
ABCS[2] = -x2 ^ 3 * x5
ABCS[3] = x2 ^ 2 * x4 * (8 * Sp * u - x6 - 24 * x8)
ABCS[4] = -(-12 * x10 ^ 4 * x4 + x10 * x4 * (x26 * (-x53 * x76 - x75 * x78) + x29 * (x65 * (8 * x30 + 8 * xi_x) - 8 * x7 * (Sc + x20 * x79 + x31 * x61 + x31 * x62) + x81 * (8 * x19 + x28) + x82 * (-8 * x1 + 4 * x13 + 4) - 8 * xi_xy)) + x10 * (x24 * x29 * (4 * x48 - x51 * (4 * x19 + 4 * xi_y)) - x26 * (x15 * x68 - x63 * (Spp * x59 + Ss + x60 * x62) + x65 * x66 + x66 * (-x20 * x70 - x64) + x74 * x75 - 4 * xi_yy)) + x11 * x12 + x11 * x5 * x9 + x12 * x14 * x16 * x9 + x16 * (x26 * (x42 ^ 2 - x43 * (Bpp * x59 + Bs + x60 * x98) + x44 * x67 + x45 * x60 + 2 * x53 ^ 2 - x60 * x99 + x60 * (-x100 * x20 - x99) + x68 * x96 + 2 * x74 ^ 2 + x74 * (2 * Fyp * u - 2 * x40) - x97 * (Fyph + Fypp * x20)) + x29 * (2 * x0 * (Gpp * x59 + Gs + x60 * x84) - x60 * x86 - x60 * (-x20 * x89 - x85) - x68 * x95 - (2 * x48 - 2 * x52) * (2 * x0 * x71 - x20 * x96 - x41))) + x18 * xi ^ 2 + x18 / x7 - x24 ^ 2 * x27 + x26 * x58 * (-x43 * (Gc + x20 * x87 + x31 * x83 + x31 * x84) + x51 * x82 + x53 * (-x38 - x94) + x57 * x86 + x60 * x91 + x78 * (x0 * x71 - x41 - x73)) + x29 * x58 * (u * (Fxph + Fxpp * x20) + u * (Fypp * x31 + Fypt) - x20 * x46 - x31 * x45 - x35 * x7 - x36 * x7 - x39 * x42 - x53 * (-x51 * x57 + x56)) + x34 * x4 * (8 * x22 - 8 * x23 + x28) + (x10 * (x26 * (x102 * x80 - x102 * x81 - x104 * x15 + x63 * (Sb - Spp * x101 + x57 * x79) + x76 * x94 + 4 * xi_xx) + x34 * (-x51 * (4 * x30 + 4 * xi_x) + 4 * x55)) + x16 * (x26 * (x103 * x44 - x104 * x96 + x106 * x57 + x39 ^ 2 + x43 * (B2 - Bpp * x101 + x105 * x57) + x46 * x57 - x57 * (-x100 * x31 - x106) + 2 * x78 ^ 2 + 2 * x94 ^ 2 - x94 * (2 * Fxp * u - 2 * x37) - x97 * (Fxpp * x31 + Fxpt)) + x29 * (2 * x0 * (Gb - Gpp * x101 + x57 * x87) - x104 * x95 + x57 * x90 - x57 * x91 + (x56 - 2 * x77) * (-x31 * x96 + x38 + 2 * x93))) - x27 * x33 ^ 2) * exp(2 * x3))
    

    nothing
end




# this is another coupled equation, for Bd and Gd. the notation used is
#
# ( A11 d_uu Bd + A12 d_uu Gd + B d_u Bd + B2 d_u Gd + C11 Bd + C12 Gd ) = -S1
# ( A21 d_uu Bd + A22 d_uu Gd + B21 d_u Bd + B22 d_u Gd + C21 Bd + C22 Gd ) = -S2

function BdGd_eq_coeff!(AA::Matrix, BB::Matrix, CC::Matrix, SS::Vector, vars::Tuple, ::Inner)
    (
        u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
        B     ,       G      ,        S      ,    Fx     ,    Fy     ,  Sd,
        Bp    ,       Gp     ,        Sp     ,    Fxp    ,    Fyp    ,
        Bpp   ,       Gpp    ,        Spp    ,    Fxpp   ,    Fypp   ,
        B_x   ,       G_x    ,        S_x    ,    Fx_x   ,    Fy_x   ,
        B_y   ,       G_y    ,        S_y    ,    Fx_y   ,    Fy_y   ,
        Bp_x  ,       Gp_x   ,        Sp_x   ,    Fxp_x  ,    Fyp_x  ,
        Bp_y  ,       Gp_y   ,        Sp_y   ,    Fxp_y  ,    Fyp_y  ,
        B_xx  ,       G_xx   ,        S_xx   ,
        B_yy  ,       G_yy   ,        S_yy   ,
                       G_xy   ,        S_xy
    ) = vars

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")

   

   x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = 8 * x2
x4 = S * x0
x5 = u * xi
x6 = x5 + 1
x7 = x4 + x6
x8 = x7 ^ 3
x9 = u * x8
x10 = x7 ^ 2
x11 = 3 * u
x12 = G * x11
x13 = u ^ 2
x14 = G * x0
x15 = x13 * tanh(x14)
x16 = 2 * x5
x17 = -Sp * x13 + x16 + 4 * x4 + 1
x18 = x8 * (-B * x11 + Bp)
x19 = sech(x14)
x20 = -Fyp
x21 = Fy * u
x22 = x20 + x21
x23 = Sh * u
x24 = 2 * x13
x25 = S * x24
x26 = x25 * xi_y
x27 = 2 * x4
x28 = x27 - 1
x29 = Fy * x28 + x23 + x26
x30 = 8 * x13
x31 = x10 * x2 * (2 * Sd * x0 + x6 ^ 2) / u
x32 = x16 + x27 + 2
x33 = u * xi_y
x34 = 5 * x13
x35 = -2 * Fyph + u * (Fy ^ 2 * x34 + 2 * Fyh + Fyp ^ 2 - 4 * Fyp * x21) - 2 * x33 * (-x20 - 2 * x21)
x36 = exp(2 * x1)
x37 = -Fxp
x38 = Fx * u
x39 = x37 + x38
x40 = St * u
x41 = x25 * xi_x
x42 = x39 * (Fx * x28 + x40 + x41)
x43 = 2 * x38
x44 = u * xi_x
x45 = x7 * (-2 * Fxpt + u * (Fx ^ 2 * x34 + Fxp ^ 2 - 4 * Fxp * x38 + 2 * Fxt) - 2 * x44 * (-x37 - x43))
x46 = Fxh * u
x47 = 3 * G
x48 = Fyt * u
x49 = Fx * Fyp
x50 = x13 * x47
x51 = Fxp * Fy
x52 = 4 * x2
x53 = sinh(x14)
x54 = 4 * x13
x55 = cosh(x14)
x56 = u ^ 4
x57 = Fx * x56
x58 = Bp * x0
x59 = 3 * B
x60 = x56 * x59
x61 = u ^ 5 * x59
x62 = 3 * x1
AA[1,1] = 0
AA[1,2] = 0
BB[1,1] = -x3 * x9
BB[1,2] = 0
CC[1,1] = -x10 * x3 * (x15 * x7 * (-Gp + x12) + x17)
CC[1,2] = x15 * x18 * x3
SS[1] = -(x0 * x19 * x52 * x7 * (-Fxh * Gp + Fxp * Gh - Fyp * Gt + Fyt * Gp - Gh * x38 + Gt * x21 + x33 * (-Fx * Gp + 3 * Fxp * G) + x46 * x47 - x47 * x48 - x49 * x50 + x50 * x51 + xi_x * (-Fyp * x12 + Gp * x21)) - x19 * x22 * x29 * x30 + x19 * x32 * x35 + x19 * x36 * (x30 * x42 - 2 * x45) - x31 * (12 * B * u - 4 * Bp))
AA[2,1] = 0
AA[2,2] = 0
BB[2,1] = 0
BB[2,2] = -x52 * x9
CC[2,1] = -x18 * x2 * x24 * sinh(2 * x14)
CC[2,2] = -x10 * x17 * x52
SS[2] = -(-x13 * x52 * x55 * (Fx * Sh * x13 - Fxp * x23 + 4 * Fy * S * x57 + Fy * St * x13 - Fy * x43 - Fyp * x40 + x22 * x41 + x26 * x39 - x27 * x49 - x27 * x51 + x49 + x51) + x2 * x32 * x55 * (-Bh * Fxp * x0 + Bh * x57 - Bt * Fy * x56 + Bt * Fyp * x0 + 5 * Fx * Fy * x0 + Fxh * x58 - Fxh * x60 + Fxp * Fyp * u - Fxph - Fypt - Fyt * x58 + Fyt * x60 - x24 * x49 - x24 * x51 + x33 * (Bp * Fx * x0 - Fxp * x62 - Fxp + x43) - x44 * (-Fyp * x62 + Fyp + x21 * (Bp * x13 - 2)) + x46 + x48 + x49 * x61 - x51 * x61) + x22 * x29 * x53 * x54 + x31 * (-6 * G * u + 2 * Gp) - x35 * x53 * x7 + x36 * x53 * (x42 * x54 - x45))
    
    nothing
end



function A_eq_coeff!(ABCS::Vector, vars::Tuple, ::Inner)
    (
        potential, u, xi, xi_x, xi_y, xi_xx, xi_yy, xi_xy,
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

    @tilde_inner("B")
    @tilde_inner("G")
    @tilde_inner("S")
    @tilde_inner("Fx")
    @tilde_inner("Fy")

    @hat_inner("B")
    @hat_inner("G")
    @hat_inner("S")
    @hat_inner("Fx")
    @hat_inner("Fy")

    @bar_inner("B")
    @bar_inner("G")
    @bar_inner("S")

    @star_inner("B")
    @star_inner("G")
    @star_inner("S")

    @tilde_inner("Fxp")
    @tilde_inner("Fyp")

    @hat_inner("Fxp")
    @hat_inner("Fyp")

    @cross_inner("G")
    @cross_inner("S")


    x0 = u ^ 3
x1 = B * x0
x2 = exp(x1)
x3 = S * x0
x4 = u * xi + 1
x5 = x3 + x4
x6 = x5 ^ 4
x7 = x2 * x6
x8 = 2 * u
x9 = 1 / u
x10 = 4 * x9
x11 = u ^ 2
x12 = S * x11
x13 = x12 + x9 + xi
x14 = Sh * u
x15 = 2 * x12
x16 = 2 * x3
x17 = x16 - 1
x18 = Fy * x17 + x14 + x15 * xi_y
x19 = G * x0
x20 = cosh(x19)
x21 = 4 * x20
x22 = x11 * x21
x23 = 3 * u
x24 = G * x23
x25 = G * x11
x26 = 3 * x25
x27 = Fy * x26 + Gh + x24 * xi_y
x28 = 4 * x27
x29 = sinh(x19)
x30 = u ^ 4
x31 = x29 * x30
x32 = Sh * x30
x33 = Sp * x11
x34 = Fy ^ 2
x35 = x11 * x34
x36 = Bh * x0
x37 = Fy * x36
x38 = 4 * x0
x39 = Fy * Sh
x40 = Sp * x30
x41 = u ^ 6
x42 = 3 * B
x43 = 2 * S * x41
x44 = u ^ 5
x45 = B * x44
x46 = x34 * x45
x47 = 8 * S * x44
x48 = 2 * x0
x49 = S * x8 - Sp
x50 = -x49
x51 = u ^ 8
x52 = x34 * x51
x53 = B * S
x54 = 6 * x53
x55 = xi_y ^ 2
x56 = 3 * x45
x57 = 4 * x14
x58 = 2 * x30
x59 = S * x58
x60 = 4 * Spp
x61 = u * x60
x62 = 3 * x1
x63 = 12 * x41
x64 = x53 * x63
x65 = 14 * x3
x66 = x33 - x65 + 1
x67 = u * xi_y
x68 = x5 ^ 2
x69 = Gp * u
x70 = 6 * Gh
x71 = Fy * x11
x72 = -Gp
x73 = -x24 - x72
x74 = Gp * x0
x75 = 6 * Bh
x76 = Fy * x44
x77 = G * x76
x78 = G * x30
x79 = 15 * x78
x80 = 2 * Gpp
x81 = 18 * G
x82 = B * u ^ 7 * x81
x83 = 4 * Gpp
x84 = Fy * x83
x85 = x1 - 1
x86 = 4 * x1
x87 = G * u
x88 = 9 * x87
x89 = Fy * u
x90 = x29 * x8
x91 = 6 * x1
x92 = 2 * x11
x93 = Bp * x92
x94 = Fy * x63
x95 = Bp * x58
x96 = 4 * Bpp
x97 = -B * x23 + Bp
x98 = x8 * x97
x99 = B ^ 2
x100 = 18 * x99
x101 = G ^ 2
x102 = 18 * x101
x103 = 18 * x41
x104 = x101 * x103 + x103 * x99
x105 = St * u
x106 = Fx * x17 + x105 + x15 * xi_x
x107 = Fx * x26 + Gt
x108 = x107 + x24 * xi_x
x109 = St * x30
x110 = Bt * x0
x111 = Fx * x110
x112 = Fx * x38
x113 = Fx ^ 2
x114 = -x11 * x113
x115 = B * Fx
x116 = x113 * x51
x117 = xi_x ^ 2
x118 = 6 * x1 + 6
x119 = u * xi_x
x120 = 6 * Fx
x121 = 6 * Gt
x122 = x115 * x44
x123 = Fx * u
x124 = Fx * G
x125 = x123 * xi_y
x126 = Fx * Fy
x127 = x0 * x126
x128 = Fx * xi_y
x129 = Fx * x11
x130 = x129 * xi_y
x131 = 3 * Bh
x132 = 3 * Bt
ABCS[1] = x7 * x8
ABCS[2] = 8 * x7
ABCS[3] = x10 * x7
ABCS[4] = -(4 * x13 ^ 4 * x2 + x13 * (4 * u * x20 * (Bh * Fy * x43 + Bh * x32 - Fyh * x16 + Fyh * x33 + Fyh - Spp * x34 * x48 - Ss * u + u * x50 * xi_yy + x34 * x40 - x34 * x47 + x35 - x37 - x38 * x39 + x39 * x41 * x42 - 3 * x46 + x52 * x54 - x55 * x8 * (-S * x56 + Spp + 3 * x12) + x67 * (Bh * x59 + Fy * (-x61 - x62 + x64 + x66) + x32 * x42 - x57)) - x18 * x28 * x31) + x18 ^ 2 * x22 + x2 * (-8 * x106 * x11 * x18 * x29 + x13 ^ 2 * (x20 * x48 * (Fxh * x26 - Fxh * x69 - Fy * Gt * x56 + Fyt * x26 - Fyt * x69 + 2 * Gc + Gh * x110 + 3 * Gh * x122 - 2 * Gp * x127 - Gp * x130 - Gt * x30 * x42 * xi_y - Gt * x36 + x121 * x67 + x121 * x71 - x124 * x131 * x44 + x125 * x83 + 30 * x126 * x78 + 27 * x128 * x19 + x129 * x70 + x129 * x84 + x132 * x77 + x132 * x78 * xi_y + xi_x * (u * (27 * Fy * x25 + Gh * x62 - Gp * x89 - x131 * x19 + x70 + x84) + xi_y * (24 * x25 + x83)) - xi_xy * (2 * Gp - 6 * x87)) + x29 * (x108 * x28 * x41 + x92 * (-Fxp + x123) * (Fyp - x89)) + (Sd * x48 + x4 ^ 2) * (8 * S * x0 - 4 * x33 - 4) / x11) - x6 * x8 * (-Bd * x20 ^ 2 * x97 - Gd * x73) + (u * x29 * (16 * Fx * S * x76 + Fxh * x16 - Fxh * x33 - Fxh + Fy * St * x38 + Fyt * x16 - Fyt * x33 - Fyt + 14 * S * x128 * x30 + Sc * x8 + Sh * x112 - Sp * x0 * x128 - Sp * x126 * x58 + 4 * St * x11 * xi_y + x119 * (Fy * (-x33 + x61 + x65 - 1) + x57 + xi_y * (12 * x12 + x60)) - x125 - x126 * x92 + x127 * x60 + x130 * x60 - x50 * x8 * xi_xy) + x20 * x30 * (x106 * x27 + x108 * x18)) * (x10 + 4 * x12 + 4 * xi)) + x68 * (-x20 * (B * Bh * x94 + Bh ^ 2 * x58 - Bs * x8 - Fyh * x91 + Fyh * x93 - Fyp ^ 2 + 2 * Fyp * x89 + G * Gh * x94 + Gh ^ 2 * x58 - x0 * x34 * x96 + x100 * x52 + x102 * x52 + x34 * x95 - x35 - 12 * x37 - 30 * x46 - x55 * (24 * B * x0 + 4 * Bpp * u - x104) + x92 * xi_y * (-27 * B * x71 + Bp * x89 - Fy * x96 + x100 * x76 + x102 * x76 + x19 * x70 + x75 * x85) + x98 * xi_yy) + x90 * (Fy * x45 * x70 - Fyh * x26 + Fyh * x69 + 2 * Gh * x36 - Gs + x34 * x74 - x34 * x79 + x34 * x82 - x35 * x80 - x55 * (12 * x25 - x45 * x81 + x80) + x67 * (x19 * x75 + x70 * x85 - x84 + x89 * (Gp + x88 * (x86 - 3))) - x70 * x71 + x73 * xi_yy + x75 * x77)) + (x106 ^ 2 * x22 + x13 * (-u * x21 * (Bt * Fx * x43 + Bt * x109 + Fxt * x16 - Fxt * x33 - Fxt + Sb * u + St * x112 + 3 * St * x115 * x41 + u * x49 * xi_xx - x111 - x113 * x40 + x113 * x47 - x113 * x56 + x114 + x116 * x54 + x117 * x118 * x3 + x119 * (Bt * x59 - Fx * (x62 - x64 + x66) + 4 * x105 + x109 * x42)) - 4 * x106 * x108 * x31) + x68 * (-x20 * (B2 * x8 + Bt ^ 2 * x58 + Bt * x115 * x63 - Fxp ^ 2 + 2 * Fxp * x123 + Fxt * x91 - Fxt * x93 + Gt ^ 2 * x58 + Gt * x124 * x63 + x100 * x116 + x102 * x116 + 12 * x111 + 30 * x113 * x45 - x113 * x95 + x114 + x117 * (24 * x1 + x104) + x92 * xi_x * (Bt * x118 + u * (-Bp * Fx + 3 * u * (G * x107 * x8 + x0 * x120 * x99 + 9 * x115))) - x98 * xi_xx) - x90 * (Bt * G * x120 * x44 + Fxt * x26 - Fxt * x69 + Gb + Gt * x11 * x120 + 2 * Gt * x110 - x113 * x74 + x113 * x79 + x113 * x82 + 6 * x117 * x25 * (x62 + 2) + x119 * (6 * Bt * x19 + Gt * x118 + x123 * (x72 + x88 * (x86 + 3))) + x121 * x122 - x73 * xi_xx))) * exp(2 * x1))

    nothing
end
