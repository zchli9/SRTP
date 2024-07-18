# WaterLily例子解析

## 1. Flow over a circle

```julia
using WaterLily # 引入WaterLily库
# 定义函数，用于创建一个圆形障碍物并初始化流体模拟环境。函数接受两个位置参数n和m（可能用于定义模拟的空间分辨率或网格大小），以及两个关键字参数Re（雷诺数，默认值为250）和U（流速，默认值为1）
function circle(n,m;Re=250,U=1)
    radius, center = m/8, m/2 # 计算圆形障碍物的半径和中心位置
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius) # 使用AutoBody定义障碍物的几何形状，该处使用了一个匿名函数来计算向量x到圆心距离减去半径的值，这里的t猜测可能是时间，目的是提供将来可能需要的灵活性，以便在CFD中引入时间依赖的行为
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body) # Simulation函数创建一个模拟环境，其中包括网格大小(n,m)、流速(U,0)、障碍物半径以及通过雷诺数计算得到的动力粘度ν
end
circ = circle(3*2^6,2^7) # circ代表初始化后的模拟环境，也即circle中的Simulation
t_end = 10 # 设置模拟的结束时间
sim_step!(circ,t_end) # sim_step!函数用于执行模拟，对模拟环境进行更新
using Plots # 引入Plots库
contour(circ.flow.p') # contour函数接受一个参数，即流场的压力分布circ.flow.p'，'表示转置，最终生成压力等高线图
```

## 2. 3D Taylor Green Vortex

```julia
using WaterLily
# TGV函数用于初始化一个流体动力学模拟，特别是泰勒-格林涡流的初始速度场，pow代表涡流的大小、Re代表雷诺数、T代表数据类型、mem代表内存结构
function TGV(; pow=6, Re=1e5, T=Float64, mem=Array)
    L = 2^pow; U = 1; ν = U*L/Re # 定义涡流大小、速度和粘性
    # 泰勒-格林涡流初始速度场函数，xyz=[x,y,z]表示位置矢量
    function uλ(i,xyz)
        x,y,z = @. (xyz-1.5)*π/L # 将坐标点xyz缩放并平移，以便涡流位于模拟的中心
        # 根据i的值（代表不同的速度分量），使用Taylor Green Vortex公式计算相应的速度分量
        i==1 && return -U*sin(x)*cos(y)*cos(z) # u_x
        i==2 && return  U*cos(x)*sin(y)*cos(z) # u_y
        return 0.                              # u_z
    end
    # 模拟初始化，包含模拟的空间尺寸(L, L, L)、初始位置(0, 0, 0)
    return Simulation((L, L, L), (0, 0, 0), L; U, uλ, ν, T, mem)
end
import CUDA # 利用NVIDIA的CUDA技术，在支持CUDA的GPU上进行加速计算
@assert CUDA.functional() # 确保CUDA环境可用；如果CUDA.functional()返回false，则程序会在这里抛出错误并停止执行
vortex = TGV(T=Float32,mem=CUDA.CuArray) # mem=CUDA.CuArray表明使用CUDA的CuArray作为数据存储结构，以便在GPU上进行高效的并行计算
sim_step!(vortex,1) # 1可能为时间步长或迭代次数
```

## 3. Moving bodies

```julia
using StaticArrays # StaticArrays提供的静态数组在进行小规模矩阵运算时，比标准的动态数组有更好的性能
# hover函数用于模拟一个物体在流体中的悬停并进行振荡和旋转运动，其中L为位置参数、Re为雷诺数、U为特征速度、amp为振幅、ϵ是一个用于控制计算精度的小的正数、thk为厚度
function hover(L=2^5;Re=250,U=1,amp=π/4,ϵ=0.5,thk=2ϵ+√2)
    # 定义线段SDF(Signed Distance Function，符号距离函数)
    function sdf(x,t)
        y = x .- SA[0,clamp(x[2],-L/2,L/2)]
        √sum(abs2,y)-thk/2
    end
    # 描述物体随时间变化的振荡运动和旋转，其中R为旋转矩阵、x为位置向量、α为旋转角度
    function map(x,t)
        α = amp*cos(t*U/L); R = SA[cos(α) sin(α); -sin(α) cos(α)]
        R * (x - SA[3L-L*sin(t*U/L),4L])
    end
    Simulation((6L,6L),(0,0),L;U,ν=U*L/Re,body=AutoBody(sdf,map),ϵ) # 创建一个模拟实例，包含模拟区域的大小、速度U、粘度ν、通过sdf和map定义的物体
end
```

## 4. Circle inside an oscillating flow

```julia
# circle函数用于在二维域中心创建一个圆形物体，并对其施加一个随时间变化的流体力。函数接受两个位置参数n和m
function circle(n,m;Re=250,U=1)
    # 在二维域中心定义一个圆
    radius = m/8 # 计算圆的半径，大小相对于域的尺寸可调
    body = AutoBody((x,t)->√sum(abs2, x .- (n/2,m/2)) - radius) # 创建一个圆形物体body，(n/2,m/2)为域中心
    # 定义一个随时间变化的体力g以及一个周期性方向perdir
    accelScale, timeScale = U^2/2radius, radius/U # 计算加速度比例accelScale和时间比例timeScale
    g(i,t) = i==1 ? -2accelScale*sin(t/timeScale) : 0 # 条件运算符，g根据索引i和时间t来决定力的大小和方向。如果i为1，则表示在x方向上施加一个正弦变化的力-2accelScale*sin(t/timeScale)，否则为0
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body, g, perdir=(1,)) # 创建一个模拟，传入域的尺寸(n,m)、初始速度(U,0)、半径、粘度、物体、力的函数、方向，其中perdir=(1,)是一个元组，指定了周期性边界条件的方向
end
```

## 5. Accelerating reference frame

```julia
using WaterLily
# 定义随时间变化的速度边界条件，模拟加速参考系中的流动
Ut(i,t::T;a0=0.5) where T = i==1 ? convert(T, a0*t) : zero(T) # 函数的返回值类型为T，这里使用了Julia的类型参数特性，where T表示T是一个泛型类型。当i为1时，返回a0*t的值并将其转换为T类型，否则返回T类型的零值
# 将Ut传递给创建模拟的函数Simulation
sim = Simulation((256,256), Ut, 32) # (256,256)元组表示模拟环境的尺寸，Ut定义了速度边界条件，32可能表示模拟的某个特定参数，例如时间步长或空间分辨率
```

## 6. Periodic and convective boundary conditions

```julia
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
```
