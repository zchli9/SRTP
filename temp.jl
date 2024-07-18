using WaterLily,StaticArrays
# sdf an map for a moving circle in y-direction
function sdf(x,t)
    norm2(SA[x[1]-192,mod(x[2]-384,384)-192])-32
end
function map(x,t)
    x.-SA[0.,t/2]
end
# make a body
body = AutoBody(sdf, map)
# y-periodic boundary conditions
Simulation((512,384), (1,0), 32; body, perdir=(2,))