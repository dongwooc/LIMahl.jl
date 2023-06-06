"""
    functions to find the real-space line-intensity power spectrum
"""

struct LineModel{F<:Function}
    νrest::Unitful.Quantity
    νobs::Unitful.Quantity
    cosmo_model::Cosmology.AbstractCosmology
    cosmo_ns::Float64
    cosmo_σ8::Float64
    cosmo_Ωb::Float64
    cosmo_Tcmb::Float64
    Mhmin_Msun::Float64
    Mhmax_Msun::Float64
    nM::Integer
    Mh_range::Vector{Float64}
    z::Float64
    hmf_σ2::Vector{Float64}
    hmf_R::Vector{Float64}
    ρm::Float64
    Pk::F
    hmf_dlnσ2dlnR::Vector{Float64}
    hmf_lnν::Vector{Float64}
    hmf_dndlnM::Vector{Float64}
    hmf_dndM::Vector{Float64}
    bofM_func::Function
    bofM::Vector{Float64}
    LofM_func::Function
    LofM::Vector{Unitful.Quantity}
    CLT::Unitful.Quantity
    Tmean_µK::Float64
    bmean::Float64
    Pshot_µK2_Mpc3::Float64
    function LineModel(νrest,νobs,cosmo_model,cosmo_ns,cosmo_σ8,cosmo_Ωb,cosmo_Tcmb,Mhmin,Mhmax,nM,LofM,bofM_func=tinker10bM)
        Mhmin_Msun,Mhmax_Msun,Mh_range,z,hmf_σ2,R,ρm,Pk,hmf_dlnσ2dlnR,lnν,hmf_dndlnM,hmf_dndM,bofM_func,bofM_eval,LofM,LofM_eval,CLT,Tmean,bmean,Pshot = LineModel_calc(νrest,νobs,cosmo_model,cosmo_ns,cosmo_σ8,cosmo_Ωb,cosmo_Tcmb,Mhmin,Mhmax,nM,LofM,bofM_func)
        F = typeof(Pk)
        new{F}(νrest,νobs,cosmo_model,cosmo_ns,cosmo_σ8,cosmo_Ωb,cosmo_Tcmb,Mhmin_Msun,Mhmax_Msun,nM,Mh_range,z,hmf_σ2,R,ρm,Pk,hmf_dlnσ2dlnR,lnν,hmf_dndlnM,hmf_dndM,bofM_func,bofM_eval,LofM,LofM_eval,CLT,Tmean,bmean,Pshot)
    end
end

function LineModel_calc(νrest::Unitful.Quantity,νobs::Unitful.Quantity,cosmo_model::Cosmology.AbstractCosmology,cosmo_ns::Float64,cosmo_σ8::Real,cosmo_Ωb::Real,cosmo_Tcmb::Real,Mhmin::Unitful.Quantity,Mhmax::Unitful.Quantity,nM::Integer,LofM::Function,bofM_func::Function)
    Mhmin_Msun = Unitful.uconvert(Unitful.NoUnits, Mhmin/UnitfulAstro.Msun)
    Mhmax_Msun = Unitful.uconvert(Unitful.NoUnits, Mhmax/UnitfulAstro.Msun)
    Mh_range = exp.(LinRange(log(Mhmin_Msun),log(Mhmax_Msun),nM))
    z = Unitful.uconvert(Unitful.NoUnits,νrest/νobs-1)
    ρm = cosmo_model.Ω_m*2.775e11*cosmo_model.h^2 # Msol/Mpc^3
    R = (Mh_range*3/(4*pi*ρm)).^(1/3) # Mpc
    #Pk, hmf_σ2, hmf_dlnσ2dlnR, lnν, hmf_dndlnM = pk_and_hmf_from_sigma8(z,cosmo_model.h,cosmo_Ωb,cosmo_model.Ω_m,cosmo_model.Ω_Λ,cosmo_ns,cosmo_σ8,cosmo_Tcmb,R)
    #Pk = pk_from_sigma8(z,cosmo_model.h,cosmo_Ωb,cosmo_model.Ω_m,cosmo_model.Ω_Λ,cosmo_ns,cosmo_σ8,cosmo_Tcmb)
    #hmf_σ2 = MatterPower.sigma2.(Pk,R)
    #hmf_dlnσ2dlnR = MatterPower.dsigma2dR.(Pk,R) .* R ./ hmf_σ2
    #lnν = 2 * log(1.686) .- log.(hmf_σ2)
    #hmf_dndlnM = -hmf_dlnσ2dlnR .* HaloMF.tinker08MF.(lnν, z, 200) ./ (Mh_range*3/ρm)
    #hmf_dndM = hmf_dndlnM ./ Mh_range
    hmf_σ2 = Vector{Float64}(undef,nM)
    lnν = Vector{Float64}(undef,nM)
    hmf_dlnσ2dlnR = Vector{Float64}(undef,nM)
    hmf_dndlnM = Vector{Float64}(undef,nM)
    hmf_dndM = Vector{Float64}(undef,nM)
    Pk = pk_and_hmf_from_sigma8!(z,cosmo_model.h,cosmo_Ωb,cosmo_model.Ω_m,cosmo_model.Ω_Λ,cosmo_ns,cosmo_σ8,cosmo_Tcmb,R, hmf_σ2, lnν, hmf_dlnσ2dlnR, hmf_dndlnM, hmf_dndM)
    LofM_eval = LofM(Mh_range.*UnitfulAstro.Msun)
    bofM_eval = bofM_func.(lnν)
    CLT = Unitful.uconvert(Unitful.µK*UnitfulAstro.Mpc^3/UnitfulAstro.Lsun,3.12e31*UnitfulAstro.Mpc^2*Unitful.µK*Unitful.km/Unitful.s^4/UnitfulAstro.Lsun*(1+z)^2/(νrest^3*Cosmology.H(cosmo_model,z)))
    Tmean = Unitful.uconvert(Unitful.NoUnits,CLT / UnitfulAstro.Mpc^3 * integrate(Mh_range,hmf_dndM.*LofM_eval) / Unitful.µK)
    bmean = Unitful.uconvert(Unitful.NoUnits,CLT / UnitfulAstro.Mpc^3 * integrate(Mh_range,hmf_dndM.*LofM_eval.*bofM_eval) / Unitful.µK)/Tmean
    Pshot = Unitful.uconvert(Unitful.NoUnits,CLT^2 / UnitfulAstro.Mpc^6 * integrate(Mh_range,hmf_dndM.*LofM_eval.^2) / Unitful.µK^2)
    Mhmin_Msun,Mhmax_Msun,Mh_range,z,hmf_σ2,R,ρm,Pk,hmf_dlnσ2dlnR,lnν,hmf_dndlnM,hmf_dndM,bofM_func,bofM_eval,LofM,LofM_eval,CLT,Tmean,bmean,Pshot
end

function Pline_realspace(l::LineModel,k::Unitful.Quantity)
    (l.Pshot_µK2_Mpc3+(l.Tmean_µK*l.bmean)^2*l.Pk(Unitful.uconvert(Unitful.NoUnits,k*UnitfulAstro.Mpc))) * Unitful.µK^2 * UnitfulAstro.Mpc^3
end
