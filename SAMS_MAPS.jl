using NCDatasets
using Statistics
using CairoMakie
using GeoMakie
using Plots
using Pkg
#using GMT
using GibbsSeaWater
using MeshArrays


##______________________ 1. INsert the theta & salt OCCA Files ________________________________________
Tfile    = "/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/DDtheta.0406clim.nc"
Sfile    = "/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/DDsalt.0406clim.nc"
dsT      = NCDataset(Tfile);dsS = NCDataset(Sfile)
println(dsT);println(dsS)
lon      = dsT["Longitude_t"][:];lat     = dsT["Latitude_t"][:]
depth    = dsT["Depth_c"][:]
depth_w  = dsT["Depth_w"][:];println(depth_w)
##_______________________ 2. SELECT THE Southern Ocean ______________________________________________
ilat     = findall(lat .<= -30)
iz       = findall(depth .<= 1000) ##selects cells whose centres are shallower than 1000 m; it does not impose an exact 1000-m boundary.
lon_SO   = lon ##southern Ocean area size(lon_SO) to see the size
##The f0 is julia indicating numbers are stored as float32
lat_SO   = lat[ilat]
depth_SO = depth[iz];println(depth_SO)##print what depths we extracted
T        = dsT["theta"][:, ilat, iz, :];size(T) ##lon, lat, depth,month
S        = dsS["salt"][:, ilat, iz, :];size(S)
dz       = diff(depth_w) ##difference between the edges of each cell

###=============================== 3.SAMW CRITERIA : LI ET AL. 2021=================================================
#Tmin  = 4.0; Tmax = 15.0 --- CONSERVATIVE TEMPERATURE
#Smin  = 34.0 ; Smax = 35.8 ----> absolute salinity 
###-------------------------------- 3a) CONVERT CONSERVATIVE TEMP to POTENTIAL TEMPERATURE ----------------------------------
CTmin = 4.0;CTmax = 15.0 
SAmin = 34.4;SAmax = 35.6# Absolute Salinity (g/kg)
# OCCA theta is potential temperature
Tmin  = GibbsSeaWater.gsw_pt_from_ct(SAmin, CTmin)
Tmax  = GibbsSeaWater.gsw_pt_from_ct(SAmax, CTmax)
println("Equivalent potential-temperature limits:")
println("Tmin = ", Tmin, " °C")
println("Tmax = ", Tmax, " °C")

##-------------------------------- 3b) CONVERT SALINITY ---------------------------------------------------
# OCCA salt is Practical Salinity.
# SA -> SP depends on:# pressure, longitude and latitude.
# Therefore Smin and Smax must be calculated at every
# lon × lat × depth point.

nlon       = length(lon_SO);
nlat       = length(lat_SO);
nz         = length(depth_SO)
Smin_local = fill(NaN, nlon, nlat, nz);Smax_local = fill(NaN, nlon, nlat, nz)
#We are converting the published SAMW salinity limits into the equivalent limits for OCCA at every grid point
### ONCE THE values are converted from TEOS to psu we then ask whether they are inside the published salinity limits
for k in 1:nz
    for j in 1:nlat
        for i in 1:nlon
            p = GibbsSeaWater.gsw_p_from_z(-depth_SO[k],lat_SO[j]) # Convert OCCA depth to pressure. GSW uses z < 0 below sea level.
            # Lower SA boundary -> local SP boundary
            Smin_local[i,j,k] =GibbsSeaWater.gsw_sp_from_sa(SAmin,p,lon_SO[i],lat_SO[j])
            # Upper SA boundary -> local SP boundary
	    Smax_local[i,j,k] =GibbsSeaWater.gsw_sp_from_sa(SAmax,p,lon_SO[i],lat_SO[j])
        end
    end
end

##-----5. Inspect the converted criteria ----------------------------
println("Potential temperature range = ",Tmin, " to ", Tmax, " °C")
println("Smin varies from ", minimum(Smin_local), " to ", maximum(Smin_local))
println("Smax varies from ", minimum(Smax_local), " to ", maximum(Smax_local))


## In order to apply the T & S limits to the 4D T & S of occa the new T& S have to be broadcasted to 4D ##
##-----6. SAMW mask -------------------------------------------------
Smin4 = reshape(Smin_local, nlon, nlat, nz, 1)
Smax4 = reshape(Smax_local, nlon, nlat, nz, 1)
samw  = coalesce.( (T .>= Tmin) .& (T .<= Tmax) .& (S .>= Smin4) .& (S .<= Smax4), false)
println(size(samw))

##__________________________________  4. THICKNESS CALCULATION _________________________________________________________________________________
dz_SO     = dz[iz];size(dz_SO)
##----Now i need to broadcast the dz_SO because Julia does not understand what to do with the size of this variable like xarray does----
dz4       = reshape(dz_SO, 1, 1, :, 1) ## reshape(array, new_dimensions...)
samw_dz   = samw .* dz4## now i multiply my smaw layers with the reshaped thickness ##
thickness = sum(samw_dz, dims=3) ##sum over the depth dimension. That 1 remains because Julia keeps the dimension you summed over.
size(thickness)
thickness = dropdims(thickness, dims=3) ##now drop the depth dimension 

##---------PLOT SAMW THICKNESS: JAN-JUNE -------------------- ##
levels = 0:50:1000
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
Z      = thickness[:, :, mon]
title  = "SAMW ($(round(Tmin, digits=1))–$(round(Tmax, digits=1)) °C, $(round(minimum(Smin_local), digits=1))–$(round(maximum(Smax_local), digits=1)) psu), $(months[mon])"

fig = CairoMakie.Figure(size = (2400, 650))
cf  = nothing
for mon in 1:6
    Z     = thickness[:, :, mon]
    ax    = CairoMakie.Axis(fig[1, mon],xlabel = "Longitude (°E)",  ylabel = mon == 1 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :lighttest)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6) # Land
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end
# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title, fontsize = 35) # Overall title
fig
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/SAMW_thickness_Jan_Jun.png", fig)


##_________________________________  5. SHALLOWEST & DEEPEST SAMW DEPTH: JAN_JUN  _______________________________________________________
nlon, nlat, nz, nmon = size(samw) ## the dimensions of the array that I will create has to have the same dimensions as samw
shallow_depth        = fill(NaN, nlon, nlat, nmon) #Create an empty array for the shallowest depth.We no longer need a depth dimension because depth will be the value stored in the array.
deep_depth           = fill(NaN, nlon, nlat, nmon) #Create an empty array for the deepest depth

for mon in 1:nmon ##loop on the months
    for j in 1:nlat ##loop on the latitudes
        for i in 1:nlon ##loop on the longitudes
            kk = findall(samw[i, j, :, mon]) ##finds all the indices where the sawm is True (we made a mask above)
            if !isempty(kk) ##continue only if the kk is NOT empty
		    shallow_depth[i, j, mon] = depth_SO[first(kk)] ##finds the first index (and then value) of True
		    deep_depth[i, j, mon]    = depth_SO[last(kk)] ##finds the last index (and then value) of False
            end
        end
    end
end

##-----------PLOT THE SHALLOWEST DEPTH of SAMW -----------------------------------------------------
levels = 0:50:1000
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
#title  = " Shallowest depth: SAMW ($(round(Tmin, digits=1))–$(round(Tmax, digits=1)) °C, $(round(minimum(Smin_local), digits=1))–$(round(maximum(Smax_local), digits=1)) psu))"
title  = "Deepest depth: SAMW ($(round(Tmin, digits=1))–$(round(Tmax, digits=1)) °C, $(round(minimum(Smin_local), digits=1))–$(round(maximum(Smax_local), digits=1)) psu)"

fig = CairoMakie.Figure(size = (2400, 650))
cf  = nothing
for mon in 1:6
    #Z  = shallow_depth[:, :, mon]
    Z     = deep_depth[:,:,mon]
    ax    = CairoMakie.Axis(fig[1, mon],xlabel = "Longitude (°E)",  ylabel = mon == 1 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :Paired_10)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6) ##Land
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end

# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title, fontsize = 35)
fig
#CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Shallowest_Depth_SAMW_Jan_Jun.png", fig)
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Deepest_Depth_SAMW_Jan_Jun.png", fig)

##____________________-____________ 5. LI ET AL 2021 CRITERION:  JULY _ DECEMBER _________________________________________________________
# SAMW THICKNESS OF JUL-DEC 
levels = 0:50:1000
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing
for mon in 7:12
    col   = mon - 6
    Z     = thickness[:, :, mon]
    ax    = CairoMakie.Axis(fig[1, col], xlabel = "Longitude (°E)", ylabel = col == 1 ? "Latitude (°N)" : "", title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z,levels = levels, colormap = :lighttest)
    CairoMakie.poly!(ax, GeoMakie.land(), color = :gray80, strokecolor = :black, strokewidth = 0.6)  # Land
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    # Keep latitude labels only on first panel
    if col != 1
        CairoMakie.hideydecorations!(ax, grid = false)
    end
end
# Horizontal colorbar spanning all six columns
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false,label = "Thickness (m)", ticks = 0:100:1000, height = 40,width = 1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title,fontsize = 35)
fig
CairoMakie.save( "/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/SAMW_thickness_Jul_Dec.png", fig)

##________________________________________- 6. SHALLOWEST AND DEEPEST SAMW JUL-DEC ____________________________________________
nlon, nlat, nz, nmon = size(samw) ## the dimensions of the array that I will create has to have the same dimensions as samw
shallow_depth        = fill(NaN, nlon, nlat, nmon) #Create an empty array for the shallowest depth.We no longer need a depth dimension because depth will be the value stored in the array.
deep_depth           = fill(NaN, nlon, nlat, nmon) #Crate an empty array for the deepest depth
for mon in 1:nmon ##loop on the months
    for j in 1:nlat ##loop on the latitudes
        for i in 1:nlon ##loop on the longitudes
            kk = findall(samw[i, j, :, mon]) ##finds all the indices where the sawm is True (we made a mask above)
            if !isempty(kk) ##continue only if the kk is NOT empty
                    shallow_depth[i, j, mon] = depth_SO[first(kk)] ##finds the first index (and then value) of True
                    deep_depth[i, j, mon]    = depth_SO[last(kk)] ##finds the last index (and then value) of False
            end
        end
    end
end
##--------------------------PLOT THE SHALLOWEST DEPTH of SAMW ---------------------------------------------------------------------
#levels = 0:50:1000
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
#title  = " Shallowest depth: SAMW ($(round(Tmin, digits=1))–$(round(Tmax, digits=1)) °C, $(round(minimum(Smin_local), digits=1))–$(round(maximum(Smax_local), digits=1)) psu))"
title  = "Deepest depth: SAMW ($(round(Tmin, digits=1))–$(round(Tmax, digits=1)) °C, $(round(minimum(Smin_local), digits=1))–$(round(maximum(Smax_local), digits=1)) psu)"
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing
for mon in 7:12
#    Z  = shallow_depth[:, :, mon]
    Z     = deep_depth[:,:,mon]
    ax    = CairoMakie.Axis(fig[1, mon-6],xlabel = "Longitude (°E)",  ylabel = mon ==7 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :Paired_10)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6)
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end
# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title, fontsize = 35)
fig
#CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Shallowest_Depth_SAMW_Jul_Dec.png", fig)
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Deepest_Depth_SAMW_Jul_Dec.png", fig)



#`````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````````
##=================================== 6. SAMW CRITERION : T & S Borreguero et Rintoul_2010  =================================================
#..............................................................................................................................................

Tmin  = 8; Tmax = 9.5
Smin  = 34.58 ; Smax = 34.68 
samw  = coalesce.((T .>= Tmin) .& (T .<= Tmax) .& (S .>= Smin) .& (S .<= Smax),false)
#true # satisfies T-S criterion; false  # does not satisfy it;missing   # T or S is missing. coalesce.(condition, false): Wherever condition is missing, replace it with false.
size(samw) ##lon,lat,depth,time
dz_SO = dz[iz];size(dz_SO)## THickenss of the 0-1000m layer ##

##----Now i need to broadcast the dz_SO because Julia does not understand what to do with the size of this variable like xarray does----
dz4       = reshape(dz_SO, 1, 1, :, 1) ## reshape(array, new_dimensions...)
samw_dz   = samw .* dz4## now i multiply my smaw layers with the reshaped thickness ##
thickness = sum(samw_dz, dims=3) #sum over the depth dimension.#That 1 remains because Julia keeps the dimension you summed over.
size(thickness)
thickness = dropdims(thickness, dims=3) ##now drop the depth dimension 
levels    = 0:50:700

##__________________________________________ 7. MAke the SAMW THICKNESS: JAN-JUNE _________________________________________________________________________
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
Z      = thickness[:, :, mon]
title  = "SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu), $(months[mon])"
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing
for mon in 1:6
    Z     = thickness[:, :, mon]
    ax    = CairoMakie.Axis(fig[1, mon],xlabel = "Longitude (°E)",  ylabel = mon == 1 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :lighttest)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6)
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end

# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],"SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)", fontsize = 35)
fig
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/SAMW_thickness_Jan_Jun_BORREGUERO.png", fig)   


#_________________________________________ 8. SHALLOWEST & DEEPEST SAMW   ____________________________________________________
nlon, nlat, nz, nmon = size(samw) ## the dimensions of the array that I will create has to have the same dimensions as samw
shallow_depth        = fill(NaN, nlon, nlat, nmon) #Create an empty array for the shallowest depth.We no longer need a depth dimension because depth will be the value stored in the array.
deep_depth           = fill(NaN, nlon, nlat, nmon) #Create an empty array for the deepest depth

for mon in 1:nmon ##loop on the months
    for j in 1:nlat ##loop on the latitudes
        for i in 1:nlon ##loop on the longitudes
            kk = findall(samw[i, j, :, mon]) ##finds all the indices where the sawm is True (we made a mask above)
            if !isempty(kk) ##continue only if the kk is NOT empty
                    shallow_depth[i, j, mon] = depth_SO[first(kk)] ##finds the first index (and then value) of True
                    deep_depth[i, j, mon]    = depth_SO[last(kk)] ##finds the last index (and then value) of False
            end
        end
    end
end


##-----------PLOT THE SHALLOWEST DEPTH of SAMW -----------------------------------------------------
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
#title  = "Shallowest depth of SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)"
title  = "Deepest depth of SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)"
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing
for mon in 1:6
    #Z  = shallow_depth[:, :, mon]
    Z     = deep_depth[:,:,mon]
    ax    = CairoMakie.Axis(fig[1, mon],xlabel = "Longitude (°E)",  ylabel = mon == 1 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :Paired_10)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6)
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end
# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title, fontsize = 35)
fig
#CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Shallowest_Depth_SAMW_Jan_Jun_BORREGUERO.png", fig)
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Deepest_Depth_SAMW_Jan_Jun__BORREGUERO.png", fig)

##--------------------------------------------------------------------------------------------------------------------------------------------
###______________________________________________________ 9. JULY _ DECEMBER _____________________________________________________________________________________
# SAMW THICKNESS OF JUL-DEC 
levels = 0:50:700
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing

for mon in 7:12
    col   = mon - 6
    Z     = thickness[:, :, mon]
    ax    = CairoMakie.Axis(fig[1, col], xlabel = "Longitude (°E)", ylabel = col == 1 ? "Latitude (°N)" : "", title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z,levels = levels, colormap = :lighttest)
    CairoMakie.poly!(ax, GeoMakie.land(), color = :gray80, strokecolor = :black, strokewidth = 0.6)
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    # Keep latitude labels only on first panel
    if col != 1
        CairoMakie.hideydecorations!(ax, grid = false)
    end
end
# Horizontal colorbar spanning all six columns
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false,label = "Thickness (m)", ticks = 0:100:1000, height = 40,width = 1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],"SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)",fontsize = 35)
fig
CairoMakie.save( "/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/SAMW_thickness_Jul_Dec_BORREGUERO.png", fig)

##___________________________________________ 10. SHALLOWEST AND DEEPEST SAMW JUL-DEC _______________________________________________________

nlon, nlat, nz, nmon = size(samw) ## the dimensions of the array that I will create has to have the same dimensions as samw
shallow_depth        = fill(NaN, nlon, nlat, nmon) #Create an empty array for the shallowest depth.We no longer need a depth dimension because depth will be the value stored in the array.
deep_depth           = fill(NaN, nlon, nlat, nmon) #Crate an empty array for the deepest depth
for mon in 1:nmon ##loop on the months
    for j in 1:nlat ##loop on the latitudes
        for i in 1:nlon ##loop on the longitudes
            kk = findall(samw[i, j, :, mon]) ##finds all the indices where the sawm is True (we made a mask above)
            if !isempty(kk) ##continue only if the kk is NOT empty
                    shallow_depth[i, j, mon] = depth_SO[first(kk)] ##finds the first index (and then value) of True
                    deep_depth[i, j, mon]    = depth_SO[last(kk)] ##finds the last index (and then value) of False
            end
        end
    end
end

##-----------PLOT THE SHALLOWEST DEPTH of SAMW -----------------------------------------------------
months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
title  = "Shallowest depth of SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)"
#title  = "Deepest depth of SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu)"
fig    = CairoMakie.Figure(size = (2400, 650))
cf     = nothing
for mon in 7:12
    Z     = shallow_depth[:, :, mon]
#    Z     = deep_depth[:,:,mon]
    ax    = CairoMakie.Axis(fig[1, mon-6],xlabel = "Longitude (°E)",  ylabel = mon ==7 ? "Latitude (°N)" : "",title = months[mon], titlesize = 34,limits = (0, 360, -80, -30))
    cf    = CairoMakie.contourf!(ax,lon_SO,lat_SO,Z;levels = levels,colormap = :Paired_10)
    CairoMakie.poly!(ax,GeoMakie.land(); color = :gray80,strokecolor = :black, strokewidth = 0.6)
    # Shifted land copy for 0–360°
    land2 = CairoMakie.poly!(ax,GeoMakie.land();color = :gray80, strokecolor = :black, strokewidth = 0.6)
    Makie.translate!(land2, 360, 0, 0)
    if mon != 1 ##this is to keep the latitutde label only in the first plot 
       CairoMakie.hideydecorations!(ax, grid = false)
    end
end

# Horizontal colorbar below all six plots
CairoMakie.Colorbar(fig[2, 1:6],cf,vertical = false, label = "Thickness (m)",ticks = 0:100:1000,height=40,width=1200,ticklabelsize = 20,labelsize = 32)
CairoMakie.Label(fig[0, 1:6],title, fontsize = 35)
fig
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Shallowest_Depth_SAMW_Jul_Dec_BORREGUERO.png", fig)
#CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/Deepest_Depth_SAMW_Jul_Dec_BORREGUERO.png", fig)


##======================================================================================================================================================
##==============================================  11.USING ECCO DATA ===================================================================================
##+=====================================================================================================================================================


###________________________________TIME MEAN SAMW + STD _________________________________________________________________________________________________
##----------- FIRST Identify which tiles have cells south of 30S ---------
gridfile = "/home/user/Desktop/hdd2/ECCO_DATA/Downloads/ECCO_L4_GEOMETRY_LLC0090GRID_V4R4/GRID_GEOMETRY_ECCO_V4r4_native_llc0090.nc"
grid     = NCDataset(gridfile)
XC       = grid["XC"][:,:,:];println(size(XC)) ##longitude
YC       = grid["YC"][:,:,:];println(size(YC)) ##latitude

for t in 1:13 ##inspect min lon and lat of each tile to find which tiles contain latitutdes southern than -30S
    println("Tile ", t, ": ", round(minimum(YC[:, :, t]), digits=2), "° to ", round(maximum(YC[:, :, t]), digits=2), "°")
end
SO_tiles = [1, 2, 4, 5, 9, 10, 12, 13]

##---------------- NOW SELECT THE DEPTHS ------------------------
Z   = grid["Z"][:]# # where T/S values are located
Zp1 = grid["Zp1"][:] ## upper and lower boundary of cells
drF = grid["drF"][:] ## how thick those levels are
println("Z size   = ", size(Z)) ;println("Zp1 size = ", size(Zp1));
println("drF size = ", size(drF))

##------------- FIND THE INDICES OF THE CElls to extract up to 1000m
iz  = findall(Z .>= -1000.0)
println("Selected indices: ", iz)
println("Number of levels: ", length(iz))
println("Last selected center: ", Z[last(iz)], " m")
println("Bottom of last selected cell: ", Zp1[last(iz)+1], " m")
println("Next center: ", Z[last(iz)+1], " m")
println("Bottom of next cell: ", Zp1[last(iz)+2], " m")

##--------------------CELL 29 is the cell for 1000m depth -----------------------------------------
iz     = 1:29
Z_SO   = Z[iz];drF_SO = drF[iz]

println("Number of levels: ", length(iz))
println("Deepest cell centre: ", Z_SO[end], " m")
println("Bottom boundary: ", Zp1[30], " m")
println("Total nominal depth: ", sum(skipmissing(drF_SO)), " m")


### ___________ EXTRACT 0-1000m from the SO_tiles ________________________________
##check dimensions of a single file
datafile = "/home/user/Desktop/hdd2/ECCO_DATA/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_TEMP_SALINITY_LLC0090GRID_DAILY_V4R4/OCEAN_TEMPERATURE_SALINITY_day_mean_2017-12-26_ECCO_V4r4_native_llc0090.nc"
ds       = NCDataset(datafile)
println(dimnames(ds["THETA"]));println(size(ds["THETA"]))
println(dimnames(ds["SALT"]));println(size(ds["SALT"]))

##------------- Southern Ocean  TILES & DEPTH OF SAMW-----------------------------
SO_tiles = [1, 2, 4, 5, 9, 10, 12, 13]
iz       = 1:29 ##these are the indices of the cells of up to 1000m
T        = ds["THETA"][:, :, SO_tiles, iz, 1];
XC_SO    = XC[:, :, SO_tiles]
S        = ds["SALT"][:, :, SO_tiles, iz, 1];
YC_SO    = YC[:, :, SO_tiles]
println("T size = ", size(T));println("S size = ", size(S))

##---------- extract tiles from the geographical tiles ---------------------------
XC_SO    = XC[:, :, SO_tiles]
YC_SO    = YC[:, :, SO_tiles]
Z_SO     = Z[iz]
drF_SO   = drF[iz]

println(size(XC_SO))   # (90, 90, 8)
println(size(YC_SO))   # (90, 90, 8)
println(size(Z_SO))    # (29,)
println(size(drF_SO))  # (29,)

##-------------- CONVERT THE Li et al., 2021 criteria to celsius and psu ------------------------------------------------
SAmin = 34.4   # g/kg
SAmax = 35.6   # g/kg
CTmin = 4.0    # °C
CTmax = 15.0   # °C

nlon  = size(XC_SO, 1);nlat  = size(XC_SO, 2)
ntile = size(XC_SO, 3);nz    = length(Z_SO)
##Smin_local and Smax_local are not ECCO salinities. We have not touched S yet.
Smin_local = fill(NaN, nlon, nlat, ntile, nz)
Smax_local = fill(NaN, nlon, nlat, ntile, nz)

for k in 1:nz
    for t in 1:ntile
        for j in 1:nlat
            for i in 1:nlon
                lon = XC_SO[i,j,t]
                lat = YC_SO[i,j,t]
                # Skip invalid/land coordinates if present
                if isfinite(lon) && isfinite(lat) #This checks that the longitude and latitude at this grid point are valid numbers
                    p     = GibbsSeaWater.gsw_p_from_z(Z_SO[k], lat) # converts the depth of level k into sea pressure.
                    spmin = GibbsSeaWater.gsw_sp_from_sa(SAmin, p, lon, lat) #At this particular longitude, latitude and pressure, what Practical Salinity (SP) would produce the Absolute Salinity SA = 34.4 g/kg
                    spmax = GibbsSeaWater.gsw_sp_from_sa(SAmax, p, lon, lat) #t this location and pressure, what SP would produce the SA = 35.6 g/k
                    if spmin < 1e10 && spmax < 1e10 ## Keep only valid GSW results because the GSW conversion fails close to Antarctica
                        Smin_local[i,j,t,k] = spmin
                        Smax_local[i,j,t,k] = spmax
		    end
                end
            end
        end
    end
end
println(size(Smin_local));println(size(Smax_local))
println("NaNs in Smin = ", count(isnan, Smin_local)) ##i cannot use extrema command because it DOES NOT IGNORE NaNs which I have put inside
println("NaNs in Smax = ", count(isnan, Smax_local))


###------------- APPLY THE Li et al. 2021 criterion on temperatures -------------------------------------------
CTmin = 4.0;CTmax = 15.0
Tmin  = GibbsSeaWater.gsw_pt_from_ct(SAmin, CTmin);println("Tmin = ", Tmin)
Tmax  = GibbsSeaWater.gsw_pt_from_ct(SAmax, CTmax);println("Tmax = ", Tmax)

###______________________________________________- NOW APPLY THE SAMW MASK __________________________________________________--
samw  = coalesce.( (T .>= Tmin) .& (T .<= Tmax) .&  (S .>= Smin_local) .& (S .<= Smax_local), false) ##the mask is True only when both conditions are satisfied
println(size(samw)) ##check that the mask works well
println("SAMW cells = ", count(samw))

##----------- Calculate the thickness of cells of the Southern Ocean  ---------------
drF_SO    = drF[iz] ##thickness of the 29 levels of the SO limit of 1000m
dz4       = reshape(drF_SO, 1, 1, 1, :)
samw_dz   = samw .* dz4
thickness = sum(samw_dz, dims=4)
thickness = dropdims(thickness, dims=4) #SAMW thickness in metres for each horizontal ECCO grid cell on this particular day.

##------------------ CREATE THE MONTHLY MEAN SAMW THICKNESS --------------------------
datadir  = "/home/user/Desktop/hdd2/ECCO_DATA/Downloads/ECCO_V4r4_PODAAC/ECCO_L4_TEMP_SALINITY_LLC0090GRID_DAILY_V4R4"
files    = sort(filter(f -> endswith(f, ".nc"), readdir(datadir, join=true)))
println("Number of daily files = ", length(files))
##------------------  initialize one accunulator for each month  ------------------------------------------
years    = 1992:2017
nyears   = length(years)   # 26
ym_sum   = zeros(Float64, 90, 90, length(SO_tiles), 12, nyears) # Sum of daily thicknesses for each month of each year
ym_count = zeros(Int, 12, nyears) ## Number of days in each month of each year

##For every day from 1992–2017 → identify SAMW → calculate its thickness → add that daily thickness to the appropriate calendar month 
#→ eventually divide by the number of days to obtain the January–December climatology.

for (n, file) in enumerate(files) #Go through every daily ECCO file
    m    = match(r"(\d{4})-(\d{2})-(\d{2})", basename(file)) # Find the date inside the filename
    m    === nothing && continue #Skip anything that doesn't contain a date
    yr   = parse(Int, m.captures[1]) #Extract the year and months and convert the strings into numbers
    mon  = parse(Int, m.captures[2])

    if yr < 1992 || yr > 2017 ## the file has to be between 1992 to 2017  || means or.
        continue
    end
    yi   = yr - 1992 + 1 # Convert year into array index:# 1992 -> 1, 1993 -> 2, ..., 2017 -> 26
    ds   = NCDataset(file) #Open that day's NetCDF file
    T    = ds["THETA"][:, :, SO_tiles, iz, 1] #We are deliberately processing one day at a time so we don't try to load 26 years of 3-D T/S data into RAM.
    S    = ds["SALT"][:, :, SO_tiles, iz, 1] #Read temperature & sal for the region/depths we selected

    samw = coalesce.((T .>= Tmin) .& (T .<= Tmax) .& (S .>= Smin_local) .& (S .<= Smax_local), false) # Identify SAMW cells. true   → this cell satisfies the SAMW criteria. false  → it doesn't
    daily_thickness = dropdims( sum(samw .* dz4, dims=4), dims=4) #Convert the SAMW mask into thickness
    ym_sum[:, :, :, mon, yi] .+= daily_thickness # # Add this day's thickness to its specific month AND year
    ym_count[mon, yi] += 1  # Count days in this specific month AND year
    close(ds) # e don't want Julia keeping thousands of NetCDF files open simultaneously.
    if n % 100 == 0
        println("Processed ", n, " / ", length(files))
    end
end

###--------- Calculate the monthly mean thickness for every year --------------------
using Statistics
ym_mean = similar(ym_sum)
for yi in 1:nyears
    for mon in 1:12
        ym_mean[:, :, :, mon, yi] =
            ym_sum[:, :, :, mon, yi] ./ ym_count[mon, yi]
    end
end

# 1992–2017 climatological monthly mean
monthly_mean = dropdims( mean(ym_mean, dims=5),dims=5)
# Interannual standard deviation of the monthly means
monthly_std  = dropdims( std(ym_mean, dims=5), dims=5)

###-------------------------------------------------- MAKE THE PLOT -------------------------------------------------------------------------
levels       = 0:50:1000
months       = ["January", "February", "March","April", "May", "June"]
lon_plot     = vec(XC_SO);lat_plot = vec(YC_SO) # Longitude and latitude are identical for every month
cmax         = maximum(monthly_mean[:, :, :, 1:6]) # Common colour scale across January-June
fig_mean_1_6 = CairoMakie.Figure(size = (1400, 750))
for mon in 1:6
    # First 3 months on row 1, next 3 on row 2
    row    = div(mon - 1, 3) + 1
    col    = mod(mon - 1, 3) + 1
    ax     = GeoMakie.GeoAxis(fig_mean_1_6[row, col], dest = "+proj=eqearth", title = months[mon])
    h_plot = vec(monthly_mean[:, :, :, mon]) # Monthly mean SAMW thickness for this month
    # Keep only finite coordinates and thickness values
    good   = isfinite.(lon_plot) .& isfinite.(lat_plot) .& isfinite.(h_plot)
    # Plot native ECCO cells
    CairoMakie.scatter!(ax,lon_plot[good],lat_plot[good];color = h_plot[good], colormap = :lighttest, colorrange = (minimum(levels), maximum(levels)), markersize = 5, marker = :rect, strokewidth = 0)
    # Grey land and black coastline
    CairoMakie.poly!( ax, GeoMakie.land(); color = :gray80, strokecolor = :black, strokewidth = 0.6)
end
# One common colorbar for all six months
CairoMakie.Colorbar( fig_mean_1_6[:, 4], limits = (minimum(levels), maximum(levels)), colormap = :lighttest, ticks = collect(levels), label = "Mean SAMW thickness (m)")
fig_mean_1_6
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/MONTHLY_MEAN_SAMW_THICKNESS_ECCO.png", fig)

###---------------------------------- MAKE THE STD PLOT --------------------------------------------------------
months       = ["January", "February", "March","April", "May", "June"]
lon_plot     = vec(XC_SO);lat_plot = vec(YC_SO) # Longitude and latitude are identical for every month
cmax         = maximum(monthly_mean[:, :, :, 1:6]) # Common colour scale across January-June
fig_mean_1_6 = CairoMakie.Figure(size = (1400, 750))
for mon in 1:6
    row    = div(mon - 1, 3) + 1  # First 3 months on row 1, next 3 on row 2
    col    = mod(mon - 1, 3) + 1
    ax     = GeoMakie.GeoAxis(fig_mean_1_6[row, col], dest = "+proj=eqearth", title = months[mon])
    h_plot = vec(monthly_std[:, :, :, mon])  # Monthly mean SAMW thickness for this month
    good   = isfinite.(lon_plot) .& isfinite.(lat_plot) .& isfinite.(h_plot) # Keep only finite coordinates and thickness values
    CairoMakie.scatter!(ax,lon_plot[good],lat_plot[good];color = h_plot[good], colormap = :lighttest, colorrange = (minimum(levels), maximum(levels)), markersize = 5, marker = :rect, strokewidth = 0) # Plot native ECCO cells
    CairoMakie.poly!( ax, GeoMakie.land(); color = :gray80, strokecolor = :black, strokewidth = 0.6)# Grey land and black coastline
end

# One common colorbar for all six months
CairoMakie.Colorbar( fig_mean_1_6[:, 4], limits = (minimum(levels), maximum(levels)), colormap = :lighttest, ticks = collect(levels), label = "Monthly STD of SAMW thickness (m)")
fig_mean_1_6
CairoMakie.save("/home/user/Desktop/hdd2/ECCO_DATA/OCCA_DATA/MONTHLY_STD_SAMW_THICKNESS_ECCO.png", fig)

























# Flatten all 8 native tiles
lon_plot = vec(XC_SO)
lat_plot = vec(YC_SO)
h_plot   = vec(thickness)

# Keep finite values
good = isfinite.(lon_plot) .&
       isfinite.(lat_plot) .&
       isfinite.(h_plot)

cmax = maximum(h_plot[good])

fig = CairoMakie.Figure(size = (1100, 650))

ax = GeoMakie.GeoAxis(
    fig[1,1],
    dest = "+proj=eqearth",
    title = "Daily SAMW thickness"
)

p = CairoMakie.scatter!(
    ax,
    lon_plot[good],
    lat_plot[good];
    color = h_plot[good],
    colormap = :viridis,
    colorrange = (0, cmax),
    markersize = 5,
    marker = :rect,
    strokewidth = 0
)

# Grey land + black coastline
CairoMakie.poly!(
    ax,
    GeoMakie.land();
    color = :gray80,
    strokecolor = :black,
    strokewidth = 0.6
)

CairoMakie.Colorbar(
    fig[1,2],
    p,
    label = "SAMW thickness (m)"
)

fig















##********************----------------------NOTES!***************************************************************************************<
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

###---------- CHeck how many gaps i have inside the SAMW water mass ----------------------------------
##_______________________________________________________________________________________________________
nlon, nlat, nz, nmon = size(samw)
ngaps = fill(NaN, nlon, nlat, nmon)

for mon in 1:nmon
    for j in 1:nlat
        for i in 1:nlon
            profile = samw[i, j, :, mon]
            kk      = findall(profile) ##find all the indices where the samw is TRUE
            if !isempty(kk)
                kfirst = first(kk) #first index of True
                klast  = last(kk) #last index of true
                inside = profile[kfirst:klast] ## keep only the indices where samw are consecutively True
#		For this longitude, latitude and month, go down through the SAMW depth range and count how many times we move from a SAMW cell (true) into a non-SAMW cell (false). Store that number in ngaps.
                ngaps[i, j, mon] = count(k -> inside[k-1] && !inside[k],2:length(inside)) ##looks for every transition from True-> False
		##is the previous cell (inside[k-1] true and (&&) the current one Flase (!)? ). The 2nd arguments shows WHICH values of k to test
            end
        end
    end
end


##Check the minimm and maximum number of SAMW gaps ##
vals = ngaps[.!isnan.(ngaps)] #Take all the values from ngaps that are not NaN, and put them into a new 1-D array called vals
println("Minimum number of gaps = ", minimum(vals))
println("Maximum number of gaps = ", maximum(vals))
println("Unique values = ", sort(unique(vals)))


##----------------------------- CHECK THE GAP THICKNESS ------------------------------------------------------------------
###______________________________________________________________________________________________________________________
gap_thickness = fill(NaN, nlon, nlat, nmon)

for mon in 1:nmon
    for j in 1:nlat
        for i in 1:nlon
            profile = samw[i, j, :, mon] ##select the samw of a specific lon and lat 
            kk      = findall(profile) ##find all the indices where samw EXISTS
            if !isempty(kk) ##if the water column is not empty 
                kfirst = first(kk)
                klast  = last(kk)
                # Start with zero gap thickness
                gap_thickness[i, j, mon] = 0.0
                # Look only between shallowest and deepest SAMW cells
                for k in kfirst:klast #look in the water column where you have identified the samw first and alst True
                    if !profile[k]#If this particular depth cell is not SAMW, do what follows.
                        gap_thickness[i, j, mon] += dz_SO[k] #Whenever we encounter a non-SAMW depth cell between the shallowest and deepest SAMW cells, add that cell's physical thickness to the total gap thickness.
                    end
                end
            end
        end
    end
end








###******************** NOTES ______________________________________________________________
###*_*_*_*_*_*_ *_*_*_*_*_*_ MAKE the ONE MONTH thickness contourf *_*_*_*_*_*_*_*_*_*_*_*_*_*_*_*___
###------------------------------------------------------------------------------------------------------
#levels = 0:50:1000
#months = ["Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"]
#mon    = 1
#Z      = thickness[:, :, mon]
#title  = "SAMW ($(Tmin)–$(Tmax) °C, $(Smin)–$(Smax) psu), $(months[mon])"
#pastel_colors = Makie.cgrad(:pastel, 20; categorical = true, rev = true)
#fig = CairoMakie.Figure(size = (1100, 500))
#ax  = CairoMakie.Axis(fig[1, 1],xlabel = "Longitude (°E)", ylabel = "Latitude (°N)", title = title,limits = (0, 360, -80, -30))
#cf  = CairoMakie.contourf!(ax,lon_SO,lat_SO, Z,levels = 0:50:1000, colormap = :lighttest)
#--- Make tha land grey and the coastline black------------------------
#This plots the land exactly as GeoMakie stores it: from −180° to +180°. Since my plot only shows 0° to 360°, I mainly see the land between 0° and 180° correctly.
#land1 = CairoMakie.poly!( ax, GeoMakie.land(); color = :gray80,  strokecolor = :black, strokewidth = 0.8)
#This creates a second identical copy of all the land. At this point, land1 and land2 are sitting directly on top of each other.
#land2 = CairoMakie.poly!( ax, GeoMakie.land(); color = :gray80, strokecolor = :black, strokewidth = 0.8)
#Makie.translate!(land2, 360, 0, 0)#This moves the second copy 360° to the right.
#CairoMakie.Colorbar(fig[1, 2],cf,label = "Thickness (m)", ticks = 0:100:1000)
#fig
#CairoMakie.save("SAMW_thickness_January.png", fig)

