# WaterLily笔记

## Julia

### 技巧

- 下载与安装
  - Julia安装包下载地址为：<https://julialang.org/downloads>
  - 安装教程：<https://www.runoob.com/julia/julia-environment.html>
- IDE选择
  - Julia语言适配的IDE为Juno，但其官网显示今后不再更新，推荐使用vscode
- vscode环境配置
  - 教程参考：<https://zhuanlan.zhihu.com/p/163809924>
- vscode运行Julia代码
  - 点击右上角三角形下的子选项Julia: Execute active File in REPL
  - 教程参考：<https://blog.csdn.net/hfy1237/article/details/124500193>
- 在vscode终端中输入`julia`即可进入REPL界面

### 语法

- 输出`println()`
- 不需要指定变量类型，它可以自动推断等号右侧的对象类型
- 命名规范
  - 变量名采用小写
  - 使用下划线`_`来分变量名中的单词
  - 类型`Type`和模块`Module`的名称使用大写字母开头，并且用大写字母而不是用下划线分隔单词
  - 函数（function）和宏（macro）的名字使用小写，不使用下划线
  - 会对输入参数进行更改的函数要使用`!`结尾。这些函数有时叫做"mutating"或"in-place"函数，因为它们在被调用后会修改他们的输入参数的内容而不仅仅只是返回一个值
- 注释
  - 单行注释`#`
  - 多行注释`#=`和`=#`

## WaterLily

- 安装package
  - 教程参考：<https://blog.csdn.net/qq_43630810/article/details/130613174>

  ```julia
    using Pkg
    Pkg.add("WaterLily")
  ```

- 模拟流体流过圆形障碍物的状态的例子

```julia
using WaterLily
function circle(n,m;Re=250,U=1)
    radius, center = m/8, m/2
    body = AutoBody((x,t)->√sum(abs2, x .- center) - radius)
    Simulation((n,m), (U,0), radius; ν=U*radius/Re, body)
end
circ = circle(3*2^6,2^7)
t_end = 10
sim_step!(circ,t_end)
using Plots
contour(circ.flow.p')
```

将上述代码用REPL运行，即可得到一个等高线图

### 定义列表

```julia
# WaterLily.jl
function Simulation(dims::NTuple{N}, u_BC, L::Number; Δt=0.25, ν=0., g=nothing, U=nothing, ϵ=1, perdir=(), uλ=nothing, exitBC=false, body::AbstractBody=NoBody(), T=Float32, mem=Array)
"""
Simulation用于创建和管理流体动力学模拟的整个过程。其中，
dims代表模拟域的维度；
u_BC代表模拟域速度边界条件，可以是一个元组或一个随时间变化的函数；
L代表模拟的长度尺度；
U代表模拟的速度尺度；
Δt代表初始时间步长；
ν代表缩放粘度；
g代表域加速度；
ϵ代表BDIM内核宽度；
perdir代表在(i,)方向的域周期性边界条件；
exitBC代表在i=1方向的对流出口边界条件；
uλ代表生成初始速度场的函数；
body代表浸没的几何体；
T代表数组元素类型；
mem代表内存位置，可以是Array、CuArray、ROCm分别对应于CPU、NVIDIA GPU或AMD GPU设备。
"""
time(sim::Simulation)
"""
返回模拟的当前累积时间。
"""
sim_time(sim::Simulation)
"""
返回模拟的当前无量纲时间。
"""
sim_step!(sim::Simulation,t_end;remeasure=true,max_steps=typemax(Int),verbose=false)
sim_step!(sim::Simulation;remeasure=true)
"""
将模拟sim集成到无量纲时间t_end。如果remeasure=true，那么物体将会在每一个时间步长被重新测量。它能够被设置成false来为静态几何体加速模拟。
"""
measure!(sim::Simulation,t=sum(sim.flow.Δt))
"""
测量一个动态物体来更新流场和泊松方程系数。
"""
# WriteVTK
vtkWriter
write!
default_attrib
pvd_collection
# ReadVTK
restart_sim!
check_nthreads(::Val{1})
__init__()
```
