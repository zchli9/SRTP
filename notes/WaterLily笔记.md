# WaterLily笔记

## Julia

### 技巧

- 下载与安装
  - Julia安装包下载地址为：https://julialang.org/downloads
  - 安装教程：https://www.runoob.com/julia/julia-environment.html
- IDE选择
  - Julia语言适配的IDE为Juno，但其官网显示今后不再更新，推荐使用vscode
- vscode环境配置
  - 教程参考：https://zhuanlan.zhihu.com/p/163809924
- vscode运行Julia代码
  - 点击右上角三角形下的子选项Julia: Execute active File in REPL
  - 教程参考：https://blog.csdn.net/hfy1237/article/details/124500193
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
  - 教程参考：https://blog.csdn.net/qq_43630810/article/details/130613174

  ```julia
    using Pkg
    Pkg.add("WaterLily")
  ```

- 一个水流流过圆圈的例子

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