
"""
- this script will create and plot an isotherm depth (8degC) in the annual mean OCCA estimate using the `isosurface` function (MeshArrays.jl)
- it downloads data on the fly if needed (Climatology.jl)
- it formats the data (360x160x50) as a MeshArray that `isosurface` knows what to do with
- `isosurface` returns a MeshArray too that `heatmap` (Makie.jl) then knows what to do with
- the same code would work on the ECCO4 or IAP grids that MeshArrays.jl also supports
"""

using MeshArrays, Climatology, CairoMakie, NetCDF

γ=MeshArrays.GridSpec(ID=:onedegree)
Γ=MeshArrays.GridLoad(γ;option="full")

msk=1.0 .*(Γ.hFacC.>0); replace!(msk,0=>NaN);

get_occa_variable_if_needed("DDtheta")

function rd(filename, varname,n)
   fil = NetCDF.open(filename, varname)
   siz = size(fil)
   tmp = zeros(siz[1:2]...,n)
   [tmp .+= fil[:,:,1:n,t] for t=1:12]
   tmp ./= 12.0
   tmp[findall(tmp.<-1e22)] .= 0.0
   return tmp
end

pth=ScratchSpaces.OCCA; n=50; fil=joinpath(pth,"DDtheta.0406clim.nc")
θ=msk*read(rd(fil,"theta",n),MeshArray(γ,Float32,n))

d=isosurface(θ,8.0,Γ)
heatmap(d)

