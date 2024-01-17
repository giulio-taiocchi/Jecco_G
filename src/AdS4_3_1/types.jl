
"""
Extend this type for different `Potential` choices
"""
abstract type Potential end

"""
Extend this type for different `InitialData` choices
"""
abstract type InitialData end

"""
Extend this type for different `Diagnostics` operations
"""
abstract type Diagnostics end

"""
Extend this type for different `GaugeCondition`s
"""
abstract type GaugeCondition end

Base.@kwdef struct ConstantGauge <: GaugeCondition
    # order of the FD operator for solving the AH equation
    fd_order :: Int = 2
end

Base.@kwdef struct Advect_xi{T} <: GaugeCondition
    xi_vx    :: T   = 0.0
    xi_vy    :: T   = 0.0
    # order of the FD operator for solving the AH equation
    fd_order :: Int = 2
end

Base.@kwdef struct ConstantAH{T} <: GaugeCondition
    u_AH     :: T   = 1.0
    kappa    :: T   = 1.0
    # order of the FD operator for solving the xi_t and AH equations
    fd_order :: Int = 2
end

"""
Parameters for the Apparent Horizon Finder
"""
Base.@kwdef struct AHF
    itmax     :: Int      = 200
    # modified from 10^-12
    epsilon   :: Float64  = 1e-10
end

"""
Extend this type for different `EvolutionEquations`. Needs members `ahf` and
`gaugecondition`.
"""
abstract type EvolutionEquations end

Base.@kwdef struct EvolTest0{TG<:GaugeCondition} <: EvolutionEquations
    gaugecondition :: TG  = ConstantGauge()
    ahf            :: AHF = AHF()
end

# Took out T and TP<:Potential, from the argument, since here we don't have any potential. How does this change run_model?
Base.@kwdef struct AffineNull{TG<:GaugeCondition} <: EvolutionEquations
    #potential      :: TP  = ZeroPotential()
    gaugecondition :: TG  = ConstantAH()
    ahf            :: AHF = AHF()
end

"""
Parameters for the time evolution
"""
Base.@kwdef struct Integration{T,Tdt,S}
    dt              :: Tdt  = :auto
    tmax            :: T
    ODE_method      :: S    = AB3()
    adaptive        :: Bool = false
    # relative tolerance for adaptive integrators
    reltol          :: Float64 = 1e-6
    filter_poststep :: Bool = true
end

"""
Parameters for Input/Output
"""
Base.@kwdef struct InOut
    # negative values suppress output
    out_boundary_every          :: Int  = -1
    out_bulk_every              :: Int  = -1
    out_gauge_every             :: Int  = -1
    out_bulkconstrained_every   :: Int  = -1

    out_boundary_every_t        :: Float64  = -1.0
    out_bulk_every_t            :: Float64  = -1.0
    out_gauge_every_t           :: Float64  = -1.0
    out_bulkconstrained_every_t :: Float64  = -1.0

    checkpoint_every_walltime_hours :: Float64 = -1.0

    # stop and checkpoint upon reaching this walltime
    max_walltime       :: Float64 = 1.e20

    # trigger termination
    termination_from_file :: Bool    = true
    check_file_every      :: Int     = 10
    termination_file      :: String  = "TERMINATE"

    # name of script
    _parfile           :: String  = splitext(basename(Base.source_path()))[1]

    # use name of script by default
    out_dir            :: String  = _parfile
    checkpoint_dir     :: String  = _parfile

    recover            :: Symbol  = :auto
    recover_dir        :: String  = _parfile

    # be very careful with this option! it will remove the whole folder contents
    # if set to true! use only for fast debugging runs
    remove_existing    :: Bool  = false
end


abstract type AbstractVars{T} <: AbstractVector{T} end

struct BulkEvolved{T} <: AbstractVars{T}
    B  :: Array{T,3}
    G   :: Array{T,3}
end

struct BulkConstrained{T} <: AbstractVars{T}
    S    :: Array{T,3}
    Fx   :: Array{T,3}
    Fy   :: Array{T,3}
    Bd  :: Array{T,3}
    Gd   :: Array{T,3}
    Sd   :: Array{T,3}
    A    :: Array{T,3}
end

struct Bulk{T} <: AbstractVars{T}
    B   :: Array{T,3}
    G    :: Array{T,3}
    S    :: Array{T,3}
    Fx   :: Array{T,3}
    Fy   :: Array{T,3}
    Bd  :: Array{T,3}
    Gd   :: Array{T,3}
    Sd   :: Array{T,3}
    A    :: Array{T,3}
end

struct BulkDeriv{T}
    Du_B   :: Array{T,3}
    Du_G    :: Array{T,3}
    Du_S    :: Array{T,3}
    Du_Fx   :: Array{T,3}
    Du_Fy   :: Array{T,3}
    Du_Sd   :: Array{T,3}
    Du_Bd  :: Array{T,3}
    Du_Gd   :: Array{T,3}
    Du_A    :: Array{T,3}
    Duu_B  :: Array{T,3}
    Duu_G   :: Array{T,3}
    Duu_S   :: Array{T,3}
    Duu_Fx  :: Array{T,3}
    Duu_Fy  :: Array{T,3}
    Duu_A   :: Array{T,3}
end

struct Boundary{T} <: AbstractVars{T}
    a3  :: Array{T,3}
    fx1 :: Array{T,3}
    fy1 :: Array{T,3}
end

struct Gauge{T} <: AbstractVars{T}
    xi  :: Array{T,3}
end

@inline varlist(::BulkEvolved)     = (:B, :G)
@inline varlist(::BulkConstrained) = (:S, :Fx, :Fy, :Bd,  :Gd, :Sd, :A)
@inline varlist(::Bulk)            = (:B,  :G, :S, :Fx, :Fy, :Bd, 
                                      :Gd,:Sd, :A)
@inline varlist(::Boundary)        = (:a3, :fx1, :fy1)
@inline varlist(::Gauge)           = (:xi,)


"""
    BulkEvolved{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are evolved in time: B, G
"""
function BulkEvolved{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    B  = Array{T}(undef, Nu, Nx, Ny)
    G   = Array{T}(undef, Nu, Nx, Ny)
    BulkEvolved{T}(B,  G)
end

"""
    BulkConstrained{T}(undef, Nu, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the bulk variables
that are constrained (not evolved in time): S, Fx, Fy, Bd, Gd, Sd, A
"""
function BulkConstrained{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    S    = Array{T}(undef, Nu, Nx, Ny)
    Fx   = Array{T}(undef, Nu, Nx, Ny)
    Fy   = Array{T}(undef, Nu, Nx, Ny)
    Bd  = Array{T}(undef, Nu, Nx, Ny)
    Gd   = Array{T}(undef, Nu, Nx, Ny)
    Sd   = Array{T}(undef, Nu, Nx, Ny)
    A    = Array{T}(undef, Nu, Nx, Ny)
    BulkConstrained{T}(S, Fx, Fy, Bd,Gd, Sd, A)
end

"""
    Bulk(bulkevol::BulkEvolved, bulkconstrain::BulkConstrained)

Construct a container to hold all the bulk variables where the evolved variables
point to the given bulkevol struct and the constrained ones point to the bulkconstrain struct
"""
function Bulk(bulkevol::BulkEvolved{T}, bulkconstrain::BulkConstrained{T}) where {T}
    B    = bulkevol.B
    G     = bulkevol.G
    S     = bulkconstrain.S
    Fx    = bulkconstrain.Fx
    Fy    = bulkconstrain.Fy
    Bd   = bulkconstrain.Bd
    Gd    = bulkconstrain.Gd
    Sd    = bulkconstrain.Sd
    A     = bulkconstrain.A
    Bulk{T}(B,  G,  S, Fx, Fy, Bd, Gd,Sd, A)
end

function BulkDeriv{T}(::UndefInitializer, Nu::Int, Nx::Int, Ny::Int) where {T<:Real}
    Du_B    = Array{T}(undef, Nu, Nx, Ny)
    Du_G     = Array{T}(undef, Nu, Nx, Ny)
    Du_S     = Array{T}(undef, Nu, Nx, Ny)
    Du_Fx    = Array{T}(undef, Nu, Nx, Ny)
    Du_Fy    = Array{T}(undef, Nu, Nx, Ny)
    Du_Sd    = Array{T}(undef, Nu, Nx, Ny)
    Du_Bd   = Array{T}(undef, Nu, Nx, Ny)
    Du_Gd    = Array{T}(undef, Nu, Nx, Ny)
    Du_A     = Array{T}(undef, Nu, Nx, Ny)
    Duu_B   = Array{T}(undef, Nu, Nx, Ny)
    Duu_G    = Array{T}(undef, Nu, Nx, Ny)
    Duu_S    = Array{T}(undef, Nu, Nx, Ny)
    Duu_Fx   = Array{T}(undef, Nu, Nx, Ny)
    Duu_Fy   = Array{T}(undef, Nu, Nx, Ny)
    Duu_A    = Array{T}(undef, Nu, Nx, Ny)
    BulkDeriv{T}(Du_B,  Du_G, Du_S, Du_Fx, Du_Fy, Du_Sd, Du_Bd,
                Du_Gd, Du_A,
                 Duu_B,  Duu_G, Duu_S, Duu_Fx, Duu_Fy, Duu_A)
end

"""
    Boundary{T}(undef, Nx, Ny)

Construct a container of uninitialized Arrays to hold all the boundary
variables: a3, fx1, fy1

These variables are automatically defined on a `(1,Nx,Ny)` grid, rather than
a `(Nx,Ny)` one, so that the same Dx and Dy differential operators defined for
the bulk quantities can also straightforwardly apply on them. Remember that the
axis along which the operator applies is defined on the operator itself. So, by
defining things this way, the Dx operator (which acts along the 2nd index) will
also do the correct thing when acting on a3, fx1, fy1.
"""
function Boundary{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    a3  = Array{T}(undef, 1, Nx, Ny)
    fx1 = Array{T}(undef, 1, Nx, Ny)
    fy1 = Array{T}(undef, 1, Nx, Ny)
    Boundary{T}(a3, fx1, fy1)
end

"""
    Gauge{T}(undef, Nx, Ny)

Construct a container with an uninitialized Array to hold the gauge variable xi

As in the Boundary struct, xi is automatically defined on a `(1,Nx,Ny)` grid,
rather than a `(Nx,Ny)` one, so that the same Dx and Dy differential operators
defined for the bulk quantities can also straightforwardly apply on it.
"""
function Gauge{T}(::UndefInitializer, Nx::Int, Ny::Int) where {T<:Real}
    xi  = Array{T}(undef, 1, Nx, Ny)
    Gauge{T}(xi)
end


@inline getB(ff::BulkEvolved)       = ff.B
@inline getG(ff::BulkEvolved)        = ff.G

@inline getS(ff::BulkConstrained)    = ff.S
@inline getFx(ff::BulkConstrained)   = ff.Fx
@inline getFy(ff::BulkConstrained)   = ff.Fy
@inline getBd(ff::BulkConstrained)  = ff.Bd
@inline getGd(ff::BulkConstrained)   = ff.Gd
@inline getSd(ff::BulkConstrained)   = ff.Sd
@inline getA(ff::BulkConstrained)    = ff.A

@inline getB(ff::Bulk)              = ff.B
@inline getG(ff::Bulk)               = ff.G
@inline getS(ff::Bulk)               = ff.S
@inline getFx(ff::Bulk)              = ff.Fx
@inline getFy(ff::Bulk)              = ff.Fy
@inline getBd(ff::Bulk)             = ff.Bd
@inline getGd(ff::Bulk)              = ff.Gd
@inline getSd(ff::Bulk)              = ff.Sd
@inline getA(ff::Bulk)               = ff.A

@inline geta3(ff::Boundary)          = ff.a3
@inline getfx1(ff::Boundary)         = ff.fx1
@inline getfy1(ff::Boundary)         = ff.fy1

@inline getxi(ff::Gauge)             = ff.xi


@inline function Base.length(ff::AbstractVars)
    vars = varlist(ff)
    sum_l = 0
    for x in vars
        f = getproperty(ff,x)   # f will point to each variable in the ff struct
        sum_l += length(f)
    end
    sum_l
end
@inline Base.size(ff::AbstractVars) = (length(ff),)

# indexing. this is just a linear indexing through all the arrays. adapted from
# RecursiveArrayTools
@inline Base.firstindex(ff::AbstractVars) = 1
@inline Base.lastindex(ff::AbstractVars) = length(ff)

@inline function Base.getindex(ff::AbstractVars, i::Int)
    vars = varlist(ff)
    @inbounds for x in vars
        f  = getproperty(ff,x)
        i -= length(f)
        if i <= 0
            return f[length(f)+i]
        end
    end
end
@inline function Base.setindex!(ff::AbstractVars, v, i::Int)
    vars = varlist(ff)
    @inbounds for x in vars
        f  = getproperty(ff,x)
        i -= length(f)
        if i <= 0
            f[length(f)+i] = v
            break
        end
    end
end

Base.similar(ff::BulkEvolved{T}) where{T} =
    BulkEvolved{T}(similar(ff.B), similar(ff.G))

Base.similar(ff::Boundary{T}) where{T} =
    Boundary{T}(similar(ff.a3), similar(ff.fx1), similar(ff.fy1))

Base.similar(ff::Gauge{T}) where{T} = Gauge{T}(similar(ff.xi))


"""
    unpack(ff::AbstractVars)

Return `Array` with all the different fields in the structure
"""
function unpack(ff::AbstractVars)
    vars = varlist(ff)
    [getproperty(ff,x) for x in vars]
end


"""
    BulkPartition{N,A} <: AbstractPartition{N,A}

Container to store (bulk) quantities that may be spread across different grid
partitions. This is to be thought as a `Tuple` (or an `Array`) of `Bulk` objects.
"""
struct BulkPartition{N,A} <: AbstractPartition{N,A}
    x :: NTuple{N,A}
end
BulkPartition(x...) = BulkPartition(tuple(x...))

Base.similar(bulks::BulkPartition{N,A}) where {N,A} = BulkPartition{N,A}(similar.(bulks.x))


function BulkEvolved(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkEvolved(bulks[i]), N)
    BulkPartition(f)
end

function BulkConstrained(bulks::BulkPartition{N}) where{N}
    f = ntuple(i -> BulkConstrained(bulks[i]), N)
    BulkPartition(f)
end


struct EvolVars{T,N,A} <: FlattenedVector{T,N,A}
    x :: NTuple{N,A}
end
EvolVars(x::NTuple{N,A}) where {N,A} = EvolVars{eltype(A),N,A}(x)
EvolVars(x...) = EvolVars(tuple(x...))

function EvolVars(ff::AbstractVector{A}) where{A<:AbstractArray}
    x = Tuple(ff)
    EvolVars{eltype(A),length(x),A}(x)
end

"""
    EvolVars(boundary::Boundary, gauge::Gauge, bulkevols::NTuple)

Build a container to store all the evolved quantities as elements of an
`NTuple`. The idea is to treat them as a single column vector for the point of
view of the time evolution routine.
"""
function EvolVars(boundary::Boundary, gauge::Gauge,
                  bulkevols::BulkPartition{Nsys}) where {Nsys}
    f1 = unpack(boundary)
    f2 = unpack(gauge)

    # this technique allocates a bit of memory and is not fully type-stable, but
    # i think it's not a problem since we only call this function once
    f3 = [unpack(bulkevol) for bulkevol in bulkevols]
    EvolVars([f1; f2; f3...])
end

Base.similar(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(similar.(ff.x))
Base.copy(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(copy.(ff.x))
Base.zero(ff::EvolVars{T,N,S}) where {T,N,S} = EvolVars{T,N,S}(zero.(ff.x))


#modified: before was div(N-4, 4). I sent the second 4 to 2 because of the number of bulkevolved functions!
# the first 4 should be related to the boundary variable plus the gauge
@inline getudomains(::EvolVars{T,N}) where {T,N} = div(N-4, 2)

@inline geta3(ff::EvolVars)   = ff.x[1]
@inline getfx1(ff::EvolVars)  = ff.x[2]
@inline getfy1(ff::EvolVars)  = ff.x[3]
@inline getxi(ff::EvolVars)   = ff.x[4]

function getB(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    Nvalue = Nsys *2+4
    #println("The value of N is $Nvalue")
    @assert i > 0
    @assert i <= Nsys
    #println("i is $i")
    #changed from ff.x[5 + (i-1)*4]
    ff.x[5 + (i-1)*2]
end


function getG(ff::EvolVars, i::Int)
    Nsys = getudomains(ff)
    @assert i > 0
    @assert i <= Nsys
    #changed from ff.x[7 + (i-1)*4]
    ff.x[6 + (i-1)*2]
end

@inline getboundary(ff::EvolVars) = Boundary(geta3(ff), getfx1(ff), getfy1(ff))
@inline getgauge(ff::EvolVars)    = Gauge(getxi(ff))

@inline getbulkevolved(ff::EvolVars, i::Int) =
    BulkEvolved(getB(ff,i), getG(ff,i))

function getbulkevolvedpartition(ff::EvolVars)
    Nsys = getudomains(ff)
    f = ntuple(i -> getbulkevolved(ff,i), Nsys)
    BulkPartition(f)
end


struct BulkHorizon{T}
    B_uAH      :: Array{T,3}
    G_uAH       :: Array{T,3}
    S_uAH       :: Array{T,3}
    Fx_uAH      :: Array{T,3}
    Fy_uAH      :: Array{T,3}
    Sd_uAH      :: Array{T,3}
    Bd_uAH     :: Array{T,3}
    Gd_uAH      :: Array{T,3}
    A_uAH       :: Array{T,3}

    Du_B_uAH   :: Array{T,3}
    Du_G_uAH    :: Array{T,3}
    Du_S_uAH    :: Array{T,3}
    Du_Fx_uAH   :: Array{T,3}
    Du_Fy_uAH   :: Array{T,3}
    Du_Sd_uAH   :: Array{T,3}
    Du_Bd_uAH  :: Array{T,3}
    Du_Gd_uAH   :: Array{T,3}
    Du_A_uAH    :: Array{T,3}

    Duu_B_uAH  :: Array{T,3}
    Duu_G_uAH   :: Array{T,3}
    Duu_S_uAH   :: Array{T,3}
    Duu_Fx_uAH  :: Array{T,3}
    Duu_Fy_uAH  :: Array{T,3}
    Duu_A_uAH   :: Array{T,3}
end
function BulkHorizon{T}(Nx::Int, Ny::Int) where {T<:Real}
    B_uAH      = Array{T}(undef, 1, Nx, Ny)
    G_uAH       = Array{T}(undef, 1, Nx, Ny)
    S_uAH       = Array{T}(undef, 1, Nx, Ny)
    Fx_uAH      = Array{T}(undef, 1, Nx, Ny)
    Fy_uAH      = Array{T}(undef, 1, Nx, Ny)
    Sd_uAH      = Array{T}(undef, 1, Nx, Ny)
    Bd_uAH     = Array{T}(undef, 1, Nx, Ny)
    Gd_uAH      = Array{T}(undef, 1, Nx, Ny)
    A_uAH       = Array{T}(undef, 1, Nx, Ny)

    Du_B_uAH   = Array{T}(undef, 1, Nx, Ny)
    Du_G_uAH    = Array{T}(undef, 1, Nx, Ny)
    Du_S_uAH    = Array{T}(undef, 1, Nx, Ny)
    Du_Fx_uAH   = Array{T}(undef, 1, Nx, Ny)
    Du_Fy_uAH   = Array{T}(undef, 1, Nx, Ny)
    Du_Sd_uAH   = Array{T}(undef, 1, Nx, Ny)
    Du_Bd_uAH  = Array{T}(undef, 1, Nx, Ny)
    Du_Gd_uAH   = Array{T}(undef, 1, Nx, Ny)
    Du_A_uAH    = Array{T}(undef, 1, Nx, Ny)

    Duu_B_uAH  = Array{T}(undef, 1, Nx, Ny)
    Duu_G_uAH   = Array{T}(undef, 1, Nx, Ny)
    Duu_S_uAH   = Array{T}(undef, 1, Nx, Ny)
    Duu_Fx_uAH  = Array{T}(undef, 1, Nx, Ny)
    Duu_Fy_uAH  = Array{T}(undef, 1, Nx, Ny)
    Duu_A_uAH   = Array{T}(undef, 1, Nx, Ny)

    BulkHorizon{T}(B_uAH, G_uAH, S_uAH, Fx_uAH, Fy_uAH,
                   Sd_uAH, Bd_uAH, Gd_uAH, A_uAH, Du_B_uAH,
                   Du_G_uAH, Du_S_uAH, Du_Fx_uAH,
                   Du_Fy_uAH, Du_Sd_uAH, Du_Bd_uAH, Du_Gd_uAH,
                   Du_A_uAH, Duu_B_uAH, Duu_G_uAH, Duu_S_uAH,
                   Duu_Fx_uAH, Duu_Fy_uAH, Duu_A_uAH)
end

struct HorizonCache{T1,T2,D}
    bulkhorizon :: BulkHorizon{T1}
    axx         :: Vector{T2}
    ayy         :: Vector{T2}
    axy         :: Vector{T2}
    bx          :: Vector{T2}
    by          :: Vector{T2}
    cc          :: Vector{T2}
    b_vec       :: Vector{T2}
    Dx_2D       :: D
    Dy_2D       :: D
    Dxx_2D      :: D
    Dyy_2D      :: D
    Dxy_2D      :: D
    _Dx_2D      :: D
    _Dy_2D      :: D
    _Dxx_2D     :: D
    _Dyy_2D     :: D
    _Dxy_2D     :: D
end
