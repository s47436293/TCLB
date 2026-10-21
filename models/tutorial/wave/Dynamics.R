AddField(name="u", stencil2d=1)
AddField(name="v", stencil2d=1)

AddQuantity(name="U")

AddSetting(name="Speed")
AddSetting(name="Viscosity")
AddSetting(name="Value", zonal=TRUE)

AddNodeType(name="Dirichlet", group="BOUNDARY")