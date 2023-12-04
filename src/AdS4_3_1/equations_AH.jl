
function AH_eq_coeff(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   ,  G   ,  S    , Fx    , Fy    , Sd ,
        B1p  , Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp , Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x , G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
	B1_y , G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
        B1p_x, Gp_x,  Sp_x , Fxp_x , Fyp_x ,
        B1p_y,  Gp_y,  Sp_y , Fxp_y , Fyp_y ,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

    @tilde_outer("B1p")
    @tilde_outer("Gp")
    @tilde_outer("Sp")
    @tilde_outer("Fxp")
    @tilde_outer("Fyp")

    @hat_outer("B1p")
    @hat_outer("Gp")
    @hat_outer("Sp")
    @hat_outer("Fxp")
    @hat_outer("Fyp")

   x0 = cosh(G)
x1 = S * x0
x2 = exp(B)
x3 = x2 / 2
x4 = exp(-B) / 2
x5 = sinh(G)
x6 = S * x5
x7 = x6 / 2
x8 = Sp * x5
x9 = sigma0_y * x8
x10 = Gp * x1
x11 = Fy + xi_y
x12 = Gp * x11
x13 = x1 / 2
x14 = Gh + x12
x15 = Bt * x1
x16 = Fxp * x1
x17 = Gt * x6
x18 = Sp * x0
x19 = 2 * x18
x20 = Bp * S
x21 = x0 * x20
x22 = 2 * Fx
x23 = Gp * x6
x24 = -Sp + x20
x25 = 2 * sigma0_x
x26 = 2 * xi_x
x27 = Fyp * x1
x28 = Gh * x6
x29 = 2 * xi_y
x30 = 2 * Fy
x31 = x23 * x30
x32 = x23 * x29
x33 = x2 * x8
x34 = x10 * x2
x35 = Sp + x20
x36 = x8 / 2
x37 = sigma0_x * x7
x38 = Spp * x5
x39 = sigma0_x * sigma0_y
x40 = Fxp * x11
x41 = Fx + xi_x
x42 = Fyp * x41
x43 = x11 * x38
x44 = sigma0_y * x41
x45 = Gp ^ 2
x46 = x12 * x18
x47 = Gp * x41
x48 = Gpp * x11
x49 = Gpp * x41
x50 = Fxh + x40 + xi_xy
x51 = Fyt + x42 + xi_xy
x52 = x14 / 2
x53 = Gt + x47
x54 = x53 / 2
x55 = x13 * (Gph + x48)
x56 = x13 * (Gpt + x49)
x57 = x18 * x52
x58 = x18 * x54
x59 = x10 / 2
x60 = Bpt * S
x61 = x0 * x60
x62 = Bt * x18
x63 = Fxpp * S
x64 = x0 * x63
x65 = Gpt * x6
x66 = Gt * x8
x67 = Fxp * x18
x68 = Spp * x0
x69 = Fx * x26
x70 = Bt * x21
x71 = Bp * x17
x72 = Bt * x23
x73 = Gt * x10
x74 = Fx ^ 2
x75 = xi_x ^ 2
x76 = Bpp * S
x77 = x0 * x76
x78 = Gpp * x6
x79 = 2 * Fxp
x80 = x21 * x79
x81 = x23 * x79
x82 = Bp ^ 2 * S
x83 = x0 * x82
x84 = S * x45
x85 = x0 * x84
x86 = 2 * Bp
x87 = x23 * x86
x88 = Gp * x86
x89 = -Spp + x82 + x84
x90 = x0 * (x76 + x89) + x6 * (Gpp + x88)
x91 = Gp * S
x92 = Gpp * S
x93 = 4 * Gp
x94 = 2 * Spp
x95 = Fypp * S
x96 = x0 * x95
x97 = Gh * x8
x98 = Gph * x6
x99 = Bh * x21
x100 = Gh * x10
x101 = Fy ^ 2
x102 = xi_y ^ 2
x103 = Fy * x29
x104 = -Gpp + x88
x105 = -x76 + x89
x106 = 2 * Fyp
axx = -x1 * x3
ayy = -x1 * x4
axy = x6
bx = Fyp * x7 + sigma0_y * x10 - x11 * x8 + x12 * x13 + x13 * x14 - x3 * (-Fx * x19 + x15 + x16 + x17 + x21 * x22 + x22 * x23 + x25 * (x0 * x24 + x23) + x26 * (-x18 + x21 + x23)) - x9
by = x4 * (Bh * x1 + Fxp * x2 * x6 + Fy * x19 + Gt * x1 * x2 + 2 * sigma0_y * (x0 * x35 - x23) + x18 * x29 + x2 * x25 * (x10 - x8) + x21 * x29 + x21 * x30 - x22 * x33 + x22 * x34 - x26 * x33 + x26 * x34 - x27 - x28 - x31 - x32)
cc = -Fxp * x9 / 2 + Fxpp * sigma0_y * x7 - Fyp * sigma0_x * x36 + Fypp * x37 + Gp * sigma0_x * x27 + Gp * sigma0_y * x16 + Gpp * x1 * x39 + S ^ 2 * Sdp + 2 * S * Sd * Sp + sigma0_x * x13 * x48 + sigma0_x * x23 * x52 - sigma0_x * x43 - sigma0_x * x46 / 2 + sigma0_x * x55 + sigma0_x * x57 + sigma0_xy * x10 + sigma0_xy * x8 + sigma0_y * x13 * x49 - sigma0_y * x18 * x47 / 2 + sigma0_y * x23 * x54 + sigma0_y * x56 + sigma0_y * x58 + x11 * x37 * x45 + x11 * x56 + x11 * x58 + x12 * x53 * x7 + x14 * x47 * x7 + x16 * x52 + x27 * x54 - x3 * (4 * Bp * Fx * x23 * xi_x + Fx * x61 + Fx * x62 + Fx * x64 + Fx * x65 + Fx * x66 - Fx * x67 + Fx * x70 + Fx * x71 + Fx * x72 + Fx * x73 + Fx * x80 + Fx * x81 + Fxp * x15 + Fxp * x17 + Fxpt * x1 + Fxt * x18 + Fxt * x21 + Fxt * x23 + sigma0_x ^ 2 * x90 + sigma0_x * (x0 * (Bt * Sp - Fx * x94 - Fxp * Sp + Gt * x91 + x20 * (Bt + x79) + x22 * x76 + x22 * x82 + x22 * x84 + x60 + x63) + x26 * x90 + x5 * (Bt * x91 + Gpt * S + Gt * Sp + x20 * (Fx * x93 + Gt) + x22 * x92 + x79 * x91)) + sigma0_xx * (x0 * x35 + x23) + x18 * xi_xx + x21 * xi_xx + x23 * xi_xx + x61 * xi_x + x62 * xi_x + x64 * xi_x + x65 * xi_x + x66 * xi_x - x67 * xi_x - x68 * x69 - x68 * x74 - x68 * x75 + x69 * x77 + x69 * x78 + x69 * x83 + x69 * x85 + x70 * xi_x + x71 * xi_x + x72 * xi_x + x73 * xi_x + x74 * x77 + x74 * x78 + x74 * x83 + x74 * x85 + x74 * x87 + x75 * x77 + x75 * x78 + x75 * x83 + x75 * x85 + x75 * x87 + x80 * xi_x + x81 * xi_x) + x36 * x50 + x36 * x51 - x38 * x39 - x38 * x44 + x39 * x45 * x6 + x4 * (Bh * Fy * Gp * S * x5 + Bh * Fy * Sp * x0 + Bh * Fyp * S * x0 + Bh * Gp * S * x5 * xi_y + Bh * Sp * x0 * xi_y + 2 * Bp * Fy * Fyp * S * x0 + Bp * Fy * Gh * S * x5 + 4 * Bp * Fy * Gp * S * x5 * xi_y + Bp * Fyh * S * x0 + 2 * Bp * Fyp * S * x0 * xi_y + Bp * Gh * S * x5 * xi_y + 2 * Bp * Gp * S * x101 * x5 + 2 * Bp * Gp * S * x102 * x5 + Bp * S * x0 * xi_yy + Bph * Fy * S * x0 + Bph * S * x0 * xi_y + 2 * Bpp * Fy * S * x0 * xi_y + Bpp * S * x0 * x101 + Bpp * S * x0 * x102 + Fy * Fyp * Sp * x0 + 2 * Fy * Spp * x0 * xi_y - Fy * x100 - Fy * x96 - Fy * x97 - Fy * x98 - Fy * x99 - Fyh * x18 - Fyh * x23 + Fyp * Sp * x0 * xi_y - Fyp * x28 - Fyp * x31 - Fyp * x32 - Fyph * x1 + Spp * x0 * x101 + Spp * x0 * x102 + sigma0_y ^ 2 * (-x0 * x105 + x104 * x6) + sigma0_y * (x0 * (-Bh * x24 + Bph * S + Fy * x94 + Fyp * Sp - Gh * x91 + x106 * x20 + x30 * x76 - x30 * x82 - x30 * x84 - x95) - x29 * (x0 * x105 - x104 * x6) + x5 * (Bh * Gp * S + Bp * S * (Fy * x93 + Gh) - Gh * Sp - Gph * S - x106 * x91 - x30 * x92)) + sigma0_yy * (x0 * x24 - x23) - x100 * xi_y - x101 * x78 - x101 * x83 - x101 * x85 - x102 * x78 - x102 * x83 - x102 * x85 - x103 * x78 - x103 * x83 - x103 * x85 - x18 * xi_yy - x23 * xi_yy - x96 * xi_y - x97 * xi_y - x98 * xi_y - x99 * xi_y) - x40 * x8 - x41 * x43 - x41 * x46 + x41 * x55 + x41 * x57 - x42 * x8 + x44 * x45 * x7 + x50 * x59 + x51 * x59 + x7 * (Fxph + Fxpp * x11) + x7 * (Fypp * x41 + Fypt)



    return axx, ayy, axy, bx, by, cc
end

function AH_eq_res(vars::Tuple, ::Outer)
    (
        sigma0, sigma0_x, sigma0_y, sigma0_xx, sigma0_yy, sigma0_xy,
        xi    , xi_x    , xi_y    , xi_xx    , xi_yy, xi_xy,
        B1   ,  G   ,  S    , Fx    , Fy    , Sd ,
        B1p  ,  Gp  ,  Sp   , Fxp   , Fyp   , Sdp,
        B1pp ,  Gpp ,  Spp  , Fxpp  , Fypp  ,
        B1_x ,  G_x ,  S_x  , Fx_x  , Fy_x  , Sd_x,
        B1_y ,  G_y ,  S_y  , Fx_y  , Fy_y  , Sd_y,
    ) = vars

    @tilde_outer("B1")
    @tilde_outer("G")
    @tilde_outer("S")
    @tilde_outer("Fx")
    @tilde_outer("Fy")
    @tilde_outer("Sd")

    @hat_outer("B1")
    @hat_outer("G")
    @hat_outer("S")
    @hat_outer("Fx")
    @hat_outer("Fy")
    @hat_outer("Sd")

   x0 = sinh(G)
x1 = S * x0
x2 = x1 / 2
x3 = Sp * x0
x4 = sigma0_x * x3
x5 = cosh(G)
x6 = S * x5
x7 = sigma0_x * x6
x8 = Fy + xi_y
x9 = Fx + xi_x
x10 = x3 * x9
x11 = Gp * x8
x12 = Gp * x9
x13 = sigma0_y * x6
x14 = Gh + x11
x15 = x14 / 2
x16 = Gt + x12
x17 = x16 / 2
x18 = x6 * x9
x19 = sigma0_y ^ 2
x20 = Sp * x5
x21 = Bp * x6
x22 = Gp * x1
x23 = Bp * x8
x24 = x6 * (Bh + x23)
x25 = x1 * x14
x26 = sigma0_x ^ 2
x27 = Bp * x9
x28 = Bt + x27
x29 = x1 * x16

Fxp * sigma0_y * x2 + Fyp * sigma0_x * x2 + Gp * sigma0_y * x7 + S ^ 2 * Sd + sigma0_xy * x1 - sigma0_y * x10 - sigma0_y * x4 - x10 * x8 + x11 * x7 / 2 + x12 * x13 / 2 + x13 * x17 + x15 * x18 + x15 * x7 + x17 * x6 * x8 + x2 * (Fxh + Fxp * x8 + xi_xy) + x2 * (Fyp * x9 + Fyt + xi_xy) - x4 * x8 + (-Fxp * x7 + 2 * Sp * sigma0_x * x5 * x9 + Sp * x26 * x5 + Sp * x5 * x9 ^ 2 - sigma0_x * x1 * x12 - sigma0_x * x29 - sigma0_xx * x6 - x18 * x28 - x21 * x26 - x22 * x26 - x27 * x7 - x28 * x7 - x29 * x9 - x6 * (Fxp * x9 + Fxt + xi_xx)) * exp(B) / 2 + (-Fyp * x13 - sigma0_y * x1 * x11 + 2 * sigma0_y * x20 * x8 + sigma0_y * x24 - sigma0_y * x25 - sigma0_yy * x6 + x13 * x23 + x19 * x20 + x19 * x21 - x19 * x22 + x20 * x8 ^ 2 + x24 * x8 - x25 * x8 - x6 * (Fyh + Fyp * x8 + xi_yy)) * exp(-B) / 2

end
