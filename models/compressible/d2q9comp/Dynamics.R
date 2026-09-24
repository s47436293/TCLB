
#Off Equillibrium compoenent for D2Q9 Lattice

AddDensity(name="g[0]", dx=0, dy=0, group="g")
AddDensity(name="g[1]", dx=1, dy=0, group="g")
AddDensity(name="g[2]", dx=0, dy=1, group="g")
AddDensity(name="g[3]", dx=-1,dy=0, group="g")
AddDensity(name="g[4]", dx=0, dy=-1,group="g")
AddDensity(name="g[5]", dx=1, dy=1, group="g")
AddDensity(name="g[6]", dx=-1,dy=1, group="g")
AddDensity(name="g[7]", dx=-1,dy=-1, group="g")
AddDensity(name="g[8]", dx=1, dy=-1, group="g")


#Fields

AddField(name="rho",   dx=c(-1,1), dy=c(-1,1), group="macro",   comment="Density")
AddField(name="ux", dx=c(-1,1), dy=c(-1,1), group="macro", comment="Velocity x")
AddField(name="uy", dx=c(-1,1), dy=c(-1,1), group="macro", comment="Velocity y")
AddField(name="theta", dx=c(-1,1), dy=c(-1,1), group="macro",   comment="Lattice temperature")
AddField(name="s",     dx=c(-2,2), dy=c(-2,2), group="entropy", comment="Entropy")
AddField(name="a1xx",  dx=c(0,0), dy=c(0,0), group="tensor",  comment="Second-order off-equilibrium moment")
AddField(name="a1xy", dx=c(0,0), dy=c(0,0), group="tensor", comment="Second-order off-equilibrium moment")
AddField(name="a1yy", dx=c(0,0), dy=c(0,0), group="tensor", comment="Second-order off-equilibrium moment")
AddField(name="Exx",  dx=c(0,0), dy=c(0,0), group="psi", comment = "Lattice Closure and polyatomic inconsistency correction sum" )
AddField(name="Eyy",  dx=c(0,0), dy=c(0,0), group="psi", comment = "Lattice Closure and polyatomic inconsistency correction sum")

#Stages

AddStage("Collision", "CalcRelaxation",
    load = DensityAll$group == "g",,
    save = Fields$group %in% c("g", "entropy", "macro"))

AddStage("Moments", "CalcMoments",
    load = DensityAll$group %in% c("g"),
    save = Fields$name %in% c("rho", "ux", "uy"),
    can.overwrite = TRUE)

AddStage("Temperature", "CalcTheta",
    load = FALSE,
    save = Fields$name == "theta",
    can.overwrite = TRUE)


AddStage("PsiEQ", "CalcEqnPsi",
    load = FALSE,
    save = Fields$group %in% c("psi"))

AddStage("A1", "CalcA1",
    load = DensityAll$group %in% c("g"),
    save = Fields$group == "tensor")

AddStage("InitMacro", "SetMacro",
    save = Fields$group %in% c("macro", "entropy"))

AddStage("InitCollide", "CalcInitCollision",
    load = FALSE,
    save = Fields$group %in% c("g", "tensor", "psi"))

AddAction("Init",      c("InitMacro", "InitCollide"))
AddAction("Iteration", c("Collision", "Moments", "Temperature", "PsiEQ", "A1"))

#Physical properties
AddSetting(name="Viscosity", default=0.0, unit="kg/m/s", comment="dynamic viscosity mu")
AddSetting(name="gamma_g",   default=1.4,  comment="heat capacity ratio")
AddSetting(name="Pr",        default=0.71, comment="Prandtl number")
AddSetting(name="T_r", default = 1.0, unit = "K", comment ="Reference temperature, lattice theta = T/T_r" )


#Scheme
AddSetting(name="sigma",   default=0.4, comment="a1 blend: sigma*a1_pr + (1-sigma)*a1_fd")
AddSetting(name="D_dim",   default=2.0, comment="spatial dimension in kappa=(D+2)/D")
AddSetting(name="chi",     default=0.33333333333333,  comment="MUSCL parameter")

#Initial condition - zonal: set per named zone painted in the XML
AddSetting(name="InitRho",   default=1.0, zonal=T, unit="kg/m3", comment="initial density")
AddSetting(name="T_inf", default=1.0, zonal=T, unit="K",     comment="initial temperature, theta = T_inf/T_r")
AddSetting(name="InitUx", default=0.0, zonal=T, unit="m/s", comment="initial velocity, x")
AddSetting(name="InitUy", default=0.0, zonal=T, unit="m/s", comment="initial velocity, y")


#Boundary Condition Nodes
AddNodeType(name="WOutflow", group="XBOUNDARY")
AddNodeType(name="EOutflow", group="XBOUNDARY")
AddNodeType(name="WBuffer",  group="XBOUNDARY")
AddNodeType(name="EBuffer",  group="XBOUNDARY")

AddNodeType(name="SOutflow", group="YBOUNDARY")
AddNodeType(name="NOutflow", group="YBOUNDARY")
AddNodeType(name="SBuffer",  group="YBOUNDARY")
AddNodeType(name="NBuffer",  group="YBOUNDARY")
AddNodeType(name="Wall",     group="YBOUNDARY")
#Outputs

AddQuantity(name="Rho")
AddQuantity(name="U", vector=TRUE)
AddQuantity(name="Theta")
AddQuantity(name="Entropy")
AddQuantity(name="Pressure")
