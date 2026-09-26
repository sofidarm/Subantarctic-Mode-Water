
hlp="""

This script will create and plot an isotherm depth (8degC) in the 
annual mean OCCA estimate using the `isosurface` function (MeshArrays.jl)

It :

- Requires `v0.5.16` (or up) of `MeshArrays.jl` that defines `layerfraction`, 
  `coldlayer`, `hotlayer`
- Downloads data on the fly if needed (Climatology.jl)
- Formats the data (360x160x50) as a MeshArray that `isosurface` etc 
  know what to do with
- `isosurface` returns a MeshArray (one level) too that `heatmap` (Makie.jl) 
  then knows what to do with
- `layerfraction`, `coldlayer`, `hotlayer` return a MeshArray (50 levels) of 
  values between 0 and 1 based on `θ>=T` for example in the case of `hotlayer`
- Should work on the other grids (ECCO4 or IAP) that MeshArrays.jl support, 
  after modifying `GridSpec` and `rd` accordingly. 
"""

if true
    using Pkg;
    Pkg.update()
    Pkg.activate(temp=true)
    Pkg.add(["MeshArrays","Climatology","CairoMakie","NetCDF"])
end

using MeshArrays, Climatology, CairoMakie, NetCDF

get_occa_variable_if_needed("DDtheta")

function rd(filename, varname; nr=50, month=1)
    fil = NetCDF.open(filename, varname)
    siz = size(fil)
    tmp = fil[:,:,1:n,month]
    tmp[findall(tmp.<-1e22)] .= 0.0
    return tmp
end

γ=MeshArrays.GridSpec(ID=:onedegree)
Γ=MeshArrays.GridLoad(γ;option="full")
msk=1.0 .*(Γ.hFacC.>0); replace!(msk,0=>NaN);

pth=ScratchSpaces.OCCA; n=50; fil=joinpath(pth,"DDtheta.0406clim.nc")
θ=read(rd(fil,"theta",month=1),MeshArray(γ,Float32,n))

d=isosurface(θ,8.0,Γ)
heatmap(d)

vol=Γ.RAC*Γ.DRF*Γ.hFacC;
m=hotlayer(θ,8.0,Γ);
m=coldlayer(θ,8.0,Γ);
m=layerfraction(θ,-Inf,8.0,Γ);
m=layerfraction(θ,7.0,8.0,Γ);
heatmap(m[:,10])

println(hlp)
