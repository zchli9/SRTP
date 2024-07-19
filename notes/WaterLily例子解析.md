# WaterLily例子解析

## 1. Flow over a circle

**功能**：模拟流体流过圆形障碍物的状态

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

**功能**：模拟3D泰勒-格林涡流

**Taylor Green Vortex**：流体动力学中常用的一个问题，用来研究湍流和数值方法的性能。它是一个经典的、理想化的流体动力学问题，具有明确的解析解，因此在验证数值模拟方法和湍流模型时非常有用。Taylor Green Vortex的初始条件在一个立方体域内定义，流体在初始时刻表现出特定的**涡流**结构。随着时间的推移，这些涡流会发生非线性相互作用，最终导致**湍流**的发展。
在一个立方体域$[0,L]^3$中，3D Taylor Green Vortex的初始速度场可以表示为
$$
u(x,y,z,0)=U_0sin(\frac{x}{L})cos(\frac{y}{L})cos(\frac{z}{L})\\
v(x,y,z,0)=-U_0cos(\frac{x}{L})sin(\frac{y}{L})cos(\frac{z}{L})\\
w(x,y,z,0)=0
$$
其中，$U_0$是初始速度幅值，$L$是立方体域的边长，$u$、$v$、$w$分别是速度场在$x$、$y$、$z$方向的分量。

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

**功能**：模拟移动和旋转的刚体在流体中的状态

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

**功能**：模拟圆形物体在振荡流体中的状态

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

**功能**：模拟流体在加速参考系中的流动

```julia
using WaterLily
# 定义随时间变化的速度边界条件，模拟加速参考系中的流动
Ut(i,t::T;a0=0.5) where T = i==1 ? convert(T, a0*t) : zero(T) # 函数的返回值类型为T，这里使用了Julia的类型参数特性，where T表示T是一个泛型类型。当i为1时，返回a0*t的值并将其转换为T类型，否则返回T类型的零值
# 将Ut传递给创建模拟的函数Simulation
sim = Simulation((256,256), Ut, 32) # (256,256)元组表示模拟环境的尺寸，Ut定义了速度边界条件，32可能表示模拟的某个特定参数，例如时间步长或空间分辨率
```

## 6. Periodic and convective boundary conditions

**功能**：模拟周期性边界条件和对流边界条件

**Periodic boundary conditions(周期性边界条件)**：假定在某个方向上，流动或物理量在边界上的值是重复的。这意味着在一个边界上的值与相对的边界上的值相同。这种条件通常用于模拟无限或周期性的问题，例如**流体在无限长的管道中的流动**。设想一个长方形区域，如果在左边界和右边界上施加周期性边界条件，那么流体在左边界流出的情况和右边界流入的情况将是相同的。

**Convective boundary conditions(对流边界条件)**：适用于描述开放边界上的流体流动。通常用于模拟边界处流体的出流或入流情况。这种边界条件考虑了流体沿着流动方向的变化，可以避免在边界处的非物理反射。在管道出口/入口处应用对流边界条件，可以**描述流体从管道中顺利排出/排入的过程**，而不会对流场产生反射影响。

- 创建周期性边界条件模拟

```julia
using WaterLily,StaticArrays
# 描述一个在y方向上移动的圆形物体
function sdf(x,t) # sdf函数返回的值代表点x到圆形物体边界的距离
    norm2(SA[x[1]-192,mod(x[2]-384,384)-192])-32 # norm2用于计算欧式距离的平方，通过SA创建了一个静态数组来优化性能，mod(x[2]-384,384)-192表示对x[2]-384取除以384的余数再减去192，使得结果在-192到191之间，可能是为了确保索引在SA数组的有效范围内
end
# map函数将原位置x沿y方向根据时间t进行平移，返回一个新的位置；结合上面的mod运算，可以实现周期性边界条件
function map(x,t)
    x.-SA[0.,t/2]
end
# 根据sdf和map创建一个物体
body = AutoBody(sdf, map)
# y方向的周期性边界条件
Simulation((512,384), (1,0), 32; body, perdir=(2,)) # 指定模拟的空间尺寸(512,384)、速度(1,0)、y方向的周期性边界条件perdir=(2,)
```

- 创建对流边界条件模拟

```julia
# 将标志exitBC=true传递给Simulation函数以启用对流边界条件
using WaterLily
body = AutoBody(sdf, map)
Simulation((512,384), u_BC=(1,0), L=32; body, exitBC=true)
```

## 7. Writing to a VTK file

**功能**：将模拟数据写入VTK文件，一般用于可视化。

- 最简单的写入器样式

```julia
# 使用WaterLily和WriteVTK库创建一个仿真sim，并将结果写入VTK文件中；VTK(Visualization Toolkit)是一种广泛用于科学数据可视化的文件格式
using WaterLily,WriteVTK # WriteVTK提供了将数据写入VTK格式文件的功能
sim = make_sim(...) # 创建一个仿真实例
# make a writer
writer = vtkwriter("simple_writer") # 通过vtkwriter函数创建一个VTK文件写入器实例simple_writer，用于将仿真数据写入文件。vtkwriter函数接受文件名或路径作为参数，用于指定输出文件的位置和名称
write!(writer,sim) # 将sim中的仿真数据写入创建的VTK文件中，!表示该函数会修改其一个或多个参数，这里表示向文件中写入数据
close(writer) # 关闭文件
```

- 采用自定义属性字典写入文件

```julia
using WaterLily,WriteVTK
# 创建一个带有一些属性(如速度、压力、体积和λ₂等)的写入器，需要输出到CPU阵列来保存文件(|> Array)，定义四个函数velocity、pressure、_body和lamda
velocity(a::Simulation) = a.flow.u |> Array; # velocity函数接受一个Simulation类型的对象作为参数，并返回一个数组，用于提取仿真对象中的特定数据，下同
pressure(a::Simulation) = a.flow.p |> Array;
_body(a::Simulation) = (measure_sdf!(a.flow.σ, a.body, WaterLily.time(a)); a.flow.σ |> Array;) # 通过measure_sdf!函数计算并提取体积分数σ
lamda(a::Simulation) = (@inside a.flow.σ[I] = WaterLily.λ₂(I, a.flow.u); a.flow.σ |> Array;)
# 创建字典custom_attrib，将字符串键映射到上述定义的函数。这个字典定义了将要写入VTK文件的自定义属性及其对应的数据提取函数
custom_attrib = Dict(
    "Velocity" => velocity,
    "Pressure" => pressure,
    "Body" => _body,
    "Lambda" => lamda
)
writer = vtkWriter("advanced_writer"; attrib=custom_attrib) # attrib=custom_attrib参数将上一步创建的custom_attrib字典传递给写入器，这样写入器就知道需要写入那些属性以及如何从仿真对象中提取这些属性的数据
...
close(writer)
```

vtkwriter函数会自动将数据写入`.pvd`文件，Paraview可以读取该文件。

```julia
# vtkwriter函数的原型
custom_vtk_function(a::Simulation) = ... |> Array
```

## 8. Restarting from a VTK file

**功能**：从VTK文件重新启动仿真。

```julia
using WaterLily,ReadVTK # ReadVTK包用于从.pvd文件中读取模拟数据
sim = make_sim(...)
writer = restart_sim!(sim; fname="file_restart.pvd") # 重新启动模拟，fname代表指定重启文件的名称
write!(writer, sim) # 将数据附加到用于重新启动的文件中
close(writer)
```

## 9.一个2D参数化壳体的例子
```julia
using WaterLily

# 定义一个函数来生成具有特定曲率的边界
function curved_square(n, m; curvature_top=-0.05, curvature_bottom=0.05, curvature_left_right=0.05)
    # 网格尺寸
    Lx, Ly = m/2, n/2
    
    # 定义 SDF 函数
    sdf(x, t) = begin
        # 初始化距离为无穷大
        dist = Inf
        
        # 顶部曲线
        curve_top(x) = curvature_top * ((x-n/2)^2-(Ly/2)^2)  + Ly / 2 + m/2
        dist_top =  x[2] - curve_top(x[1]) 
        
        # 底部曲线
        curve_bottom(x) = curvature_bottom * ((x-n/2)^2-(Ly/2)^2) - Ly / 2 + m/2
        dist_bottom = curve_bottom(x[1]) - x[2]  
        
        # 左侧曲线
        curve_left(y) = -curvature_left_right * ((y-m/2)^2-(Lx/2)^2) - Lx / 2 + n/2
        dist_left = curve_left(x[2]) - x[1]
        
        # 右侧曲线
        curve_right(y) = curvature_left_right * ((y-m/2)^2-(Lx/2)^2) + Lx / 2 + n/2
        dist_right = x[1] - curve_right(x[2]) 
        
        # 返回最距离作为 SDF
        max(dist_top, dist_bottom, dist_left, dist_right)
    end
    
    # 创建 AutoBody 对象
    body = AutoBody(sdf)
    
    return body
end
function sim_generate(n, m; curvature_top=-0.01, curvature_bottom=0.01, curvature_left_right=0.01)
    body = curved_square(n, m; curvature_top=curvature_top, curvature_bottom=curvature_bottom, curvature_left_right=curvature_left_right)
    sim = Simulation((n, m), (1,0), m; ν=1e-3, body=body)
    return sim
end
include("TwoD_plots.jl")
sim_gif!(sim_generate(64,64), duration=10, clims=(-5, 5), plotbody=true)
```
使用SDF函数模拟曲线，组合成几何体，这里使用了和叶学长论文一致的三个曲率参数来表示，可以调整参数试试效果。
