"""
    functions to find the real-space line-intensity power spectrum
"""

struct LineModel{T<:Real,S<:Unitful.AbstractQuantity,M<:Unitful.AbstractQuantity,C<:Cosmology.AbstractCosmology,I<:Integer,F<:Function}
    νrest::S
    νobs::S
    cosmo_model::C
    cosmo_ns::T
    cosmo_σ8::T
    cosmo_Ωb::T
    cosmo_Tcmb::T
    Mhmin_Msun::T
    Mhmax_Msun::T
    nM::I
    Mh_range::Any
    z::T
    hmf_σ2::Any
    hmf_R::Any
    ρm::T
    Pk::Any
    hmf_dlnσ2dlnR::Any
    hmf_lnν::Any
    hmf_dndlnM::Any
    hmf_dndM::Any
    LofM_func::F
    LofM::Any
    CLT::Any
    Tmean::Any
    function LineModel(νrest::S,νobs::S,cosmo_model::C,cosmo_ns::T,cosmo_σ8::T,cosmo_Ωb::T,cosmo_Tcmb::T,Mhmin::M,Mhmax::M,nM::I,LofM::F) where {T,S,M,C,I,F}
        # want to check for Quantity dimension; not sure how to do this yet
        Mhmin_Msun = Unitful.uconvert(Unitful.NoUnits, Mhmin/UnitfulAstro.Msun)
        Mhmax_Msun = Unitful.uconvert(Unitful.NoUnits, Mhmax/UnitfulAstro.Msun)
        Mh_range = exp.(LinRange(log(Mhmin_Msun),log(Mhmax_Msun),nM))
        z = Unitful.uconvert(Unitful.NoUnits,νrest/νobs-1)
        ρm = cosmo_model.Ω_m*2.775e11*cosmo_model.h^2 # Msol/Mpc^3
        R = (Mh_range*3/(4*pi*ρm)).^(1/3) # Mpc
        Pk = pk_from_sigma8(z,cosmo_model.h,cosmo_Ωb,cosmo_model.Ω_m,cosmo_model.Ω_Λ,cosmo_ns,cosmo_σ8,cosmo_Tcmb)
        hmf_σ2 = MatterPower.sigma2.(Pk,R)
        hmf_dlnσ2dlnR = MatterPower.dsigma2dR.(Pk,R) .* R ./ hmf_σ2
        lnν = 2 * log(1.686) .- log.(hmf_σ2)
        hmf_dndlnM = -hmf_dlnσ2dlnR .* HaloMF.tinker08MF.(lnν, z, 200) ./ (Mh_range*3/ρm)
        hmf_dndM = hmf_dndlnM ./ Mh_range
        LofM_eval = LofM(Mh_range.*UnitfulAstro.Msun)
        CLT = Unitful.uconvert(Unitful.µK*UnitfulAstro.Mpc^3/UnitfulAstro.Lsun,3.12e31*UnitfulAstro.Mpc^2*Unitful.µK*Unitful.km/Unitful.s^4/UnitfulAstro.Lsun*(1+z)^2/(νrest^3*Cosmology.H(cosmo_model,z)))
        Tmean = Unitful.uconvert(Unitful.µK,CLT / UnitfulAstro.Mpc^3 * integrate(Mh_range,hmf_dndM.*LofM_eval))
        new{T,S,M,C,I,F}(νrest,νobs,cosmo_model,cosmo_ns,cosmo_σ8,cosmo_Ωb,cosmo_Tcmb,Mhmin_Msun,Mhmax_Msun,nM,Mh_range,z,hmf_σ2,R,ρm,Pk,hmf_dlnσ2dlnR,lnν,hmf_dndlnM,hmf_dndM,LofM,LofM_eval,CLT,Tmean)
    end
end
