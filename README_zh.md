# DynamicAlphaHull

[English documentation](README.md)

`DynamicAlphaHull` 用于根据清洗后的 WGS 84 物种出现点构建地理分布范围。它基于
Delaunay 三角剖分构造 alpha-hull，自动选择满足点覆盖率与分区数要求的 alpha，并返回
可直接用于空间分析的 `terra::SpatVector` 多边形。

它适用于希望用比最小凸包更好地保留凹陷、孔隙和离散分布的范围多边形，同时不想维护
旧脚本或在运行时联网下载海岸数据的场景。

## 项目背景

本包将两套较早的 R 工作流现代化：`alphahull` 的 Delaunay、alpha-shape 与
alpha-hull 几何算法，以及 `rangeBuilder` 根据点覆盖率和多边形分区数自动选择 alpha
的范围构建逻辑。

这里不是简单封装旧包，而是将需要的功能重组为独立的现代 R 包：

- 使用 `terra` 统一完成矢量几何、投影、缓冲与叠加；包内不使用 `sf`；
- 不依赖 `alphahull`、`rangeBuilder`、`rnaturalearth` 或网络服务；
- 将预融合的 Natural Earth 1:50m 陆地数据作为 lazy internal data 内置，海岸裁切可完全离线进行；
- 保留 alpha-hull 的几何语义与动态选择逻辑；Delaunay 三角剖分与 `alphahull` 一样取自
  `interp`，使输出与参考实现逐位一致。

本包专注于从已清洗的出现记录构建 alpha-hull 范围；它不是通用 GIS 框架，也不负责
分类学、坐标质量或异常点清洗。

## 一致性与性能

几何算法按参考实现 `alphahull` / `rangeBuilder` 移植，保证输出完全一致。尤其是
`delvor()` 使用 `interp::tri.mesh` —— 与 `alphahull` 完全相同的单精度三角剖分。这是
有意的：三角剖分的边顺序决定 alpha-hull 的圆弧顺序，因此必须使用同一三角剖分才能保证
多边形逐位一致。基于 `deldir` 的三角剖分虽然得到相同的边集，但顺序不同，最终多边形也
不同，故 `delvor()` 不再使用它。

其余性能优化在不改变结果的前提下保留：

- `ahull()` 先按圆心距离筛除不可能相交的圆弧对，只让可能相交的组合进入原有裁剪状态机；
- `ashape()` 用向量化和 `data.table` 聚合替换逐边 `rank()` 与临时表；
- `getDynamicAlphaHull()` 每轮 alpha 都重新三角化，与 `rangeBuilder` 完全一致，使固定
  种子下消耗相同的 `interp` jitter 序列；
- `getDynamicAlphaHull()` 的近似重复点距离矩阵只构建一次，之后增量维护：丢弃一个点不会改变
  其余点之间的距离，因此每轮只需删掉被丢弃点的行列，不再整体重建 O(n²) 矩阵。最近点对的
  选择逻辑不变，丢弃的点序列逐位相同；
- `ah2terra()` 移植 `rangeBuilder::ah2sf` 的圆弧重排与环闭合算法；`dummycoor()` 与
  `alphahull` 一样使用 `interp::in.convex.hull`。

下面把本包的 `getDynamicAlphaHull()` 与它取代的传统 `alphahull` + `rangeBuilder`
流程直接对比：同一批真实出现记录、相同参数、单 R 进程（无并行 worker），每次调用前
`set.seed(1)`，使两个实现消耗完全相同的 RNG 流并选出相同的 alpha 序列。海岸裁切排除
在外（`clipToCoast = "no"`），与参考管线的调用方式一致。用 `microbenchmark` 实测
（R 4.6，Windows）：

| 物种 | 点数 | alphahull + rangeBuilder (s) | DynamicAlphaHull (s) | 提速 | alpha |
| --- | ---: | ---: | ---: | ---: | ---: |
| Rubus idaeus | 3,995 | 153.97 | 5.62 | 27x | alpha15 |
| Potentilla nivea | 1,654 | 24.72 | 5.33 | 4.6x | alpha30 |
| Neillia incisa | 1,146 | 3.43 | 0.43 | 7.9x | alpha4 |
| Ceanothus fendleri | 801 | 2.35 | 0.17 | 14x | alpha2 |

Rubus 与 Potentilla 较重：自适应搜索每轮 alpha 都在 `interp` 里重新三角化（分别 15 次和
30 次），且序列受与参考实现一致的 RNG 对等性约束固定，这限制了最大物种还能提升的上限。
两个实现为每个物种选出相同的 alpha（表中可见）；多边形等价性由下面的 2,296 物种验证覆盖。

在 2,296 个真实物种与 rangeBuilder 基线的对比中，2,292 个（99.8%）逐位一致，面积差异
中位数 0.001%。实际耗时随点构型、alpha 序列和机器而变化。

## 椭球几何与四个不一致的物种

本包用 `terra` 的 GeographicLib 椭球度量距离与面积，并在与 alpha-hull 构建相同的平面
空间里做点是否在多边形内的判定。而参考实现 `rangeBuilder` 的基线是在 `sf` 的 S2
球面几何开启状态下运行的——S2 把地球当成正球体。椭球更精确（地球在赤道隆起），因此两者
有出入时，本包的结果更精确。

2,296 个物种中有 4 个与冻结的基线不一致，原因是库层面的差异，而非算法错误：

- `Pouzolzia guatemalana var. nivea` —— 覆盖率检查里 S2 球面 vs 平面的点在多边形内判定；
- `Brosimum rubescens` 与 `Ficus lutea` —— 删除近似重复点时 S2/LWGEOM 椭球测地线 vs
  本包的大圆距离，在近乎重合的记录上有亚毫米级差异；
- `Potentilla elegans` —— hull 跨越反子午线，其自交在 `sf` 与 `terra` 各自链接的 GEOS
  版本里被判为相反的有效性。

此外 `Ficus lutea` 还会让参考实现直接崩溃（未处理的 `try-error` 传给了 `sf`），而本包
会返回最小凸包。

## 安装

从 GitHub 安装开发版本：

```r
install.packages("remotes")
remotes::install_github("wyx619/DynamicAlphaHull")
```

运行时依赖 `data.table`、`interp`、`purrr` 与 `terra`，缺失的 R 依赖会自动安装。
某些平台安装 `terra` 时还需要 GDAL/PROJ 等系统组件。

## 快速开始

最小输入为含经度、纬度的 data frame 或矩阵。下面构建范围并保留陆地部分：

```r
library(DynamicAlphaHull)

occurrences <- data.frame(
  Longitude = c(116.0, 116.5, 116.4, 116.1),
  Latitude = c(39.0, 39.0, 39.4, 39.3)
)

range <- getDynamicAlphaHull(
  occurrences,
  buff = 1000,
  clipToCoast = "terrestrial"
)

range$alpha
#> [1] "alpha3"

plot(range$hull, col = "lightblue", border = "steelblue")
points(occurrences$Longitude, occurrences$Latitude, pch = 16)
```

本包已为输出对象注册 base R 的 `plot()` 方法。执行
`library(DynamicAlphaHull)` 后可直接使用 `plot(range$hull)`，无需为了绘图再
`library(terra)`。

返回列表包含：

- `range$hull`：WGS 84 坐标系下的 `terra::SpatVector` 多边形；
- `range$alpha`：实际选用的 alpha 标签，如 `"alpha0.07"`；若使用最小凸包降级结果，
  则为 `"alphaMCH"`。

## 动态 alpha 选择

`getDynamicAlphaHull()` 从 `initialAlpha` 开始，每次增加 `alphaIncrement`，直到：

1. 产生有效多边形；
2. 多边形分区数不超过 `partCount`；
3. 范围与至少 `fraction` 比例的输入点相交。

若到达 `alphaCap` 仍不能满足条件，函数返回带缓冲的最小凸包（MCH）；共线输入也直接
使用这一降级路径。`verbose = TRUE` 会打印每次尝试的 alpha。

```r
set.seed(20260811)
occurrences <- data.frame(
  Longitude = runif(500, 115.8, 116.7),
  Latitude = runif(500, 38.8, 39.6)
)

range <- getDynamicAlphaHull(
  occurrences,
  fraction = 0.98,
  partCount = 1,
  buff = 1000,
  initialAlpha = 0.005,
  alphaIncrement = 0.005,
  alphaCap = 0.08,
  clipToCoast = "no",
  verbose = TRUE
)

range$alpha
#> [1] "alpha0.07"
```

### 参数与单位

| 参数 | 含义 |
| --- | --- |
| `fraction` | 必须被范围覆盖的点比例，取值 `(0, 1]`。 |
| `partCount` | 允许的最大不相连多边形数。 |
| `buff` | 最终范围缓冲距离，单位为米；在 Equal Earth（EPSG:8857）中执行。 |
| `initialAlpha`、`alphaIncrement`、`alphaCap` | alpha 搜索序列；alpha 直接在输入经纬度平面上计算，单位为“度”。 |
| `clipToCoast` | `"no"` 不裁切；`"terrestrial"` 仅保留陆地；`"aquatic"` 仅保留海洋。 |

alpha 不是米制或测地距离。该包沿用经纬度平面上的 alpha-hull 定义，因而对大范围、
跨越 180° 经线或接近极地的数据应谨慎解释；区域尺度数据应按点间距试验 alpha。

### 自定义列名和矩阵输入

当经纬度不在前两列时，用 `coordHeaders` 指定；若输入恰有两列，则直接按第一、二列
解释：

```r
records <- data.frame(
  species = "example_species",
  lon = c(116.0, 116.5, 116.4, 116.1),
  lat = c(39.0, 39.0, 39.4, 39.3)
)

range <- getDynamicAlphaHull(records, coordHeaders = c("lon", "lat"), clipToCoast = "no")
matrixRange <- getDynamicAlphaHull(as.matrix(records[, c("lon", "lat")]))
```

缺失值、非有限值和完全重复坐标会被移除；少于三个有效唯一坐标会报错。

## 离线海岸裁切

`clipToCoast = "terrestrial"` 或 `"aquatic"` 使用内置 Natural Earth 1:50m 陆地
图层。该图层已预先融合为 lazy internal data，只在当前 R 会话首次使用时恢复并缓存；
它不访问网络、不下载数据、也不写出外部文件：

```r
land <- loadWorldMap()
plot(land, col = "grey85", border = NA)
plot(range$hull, add = TRUE, border = "firebrick", lwd = 2)
```

## 底层几何接口

大多数使用者只需要 `getDynamicAlphaHull()`。如需固定 alpha、检查中间对象或对多个
alpha 复用同一 Delaunay 网格，可使用：

```r
coordinates <- as.matrix(occurrences[, c("Longitude", "Latitude")])
mesh <- delvor(coordinates)
shape <- ashape(mesh, alpha = 0.05)
hull <- ahull(mesh, alpha = 0.05)
polygon <- ah2terra(hull)

plot(polygon)
```

将已有 `delvor` 对象传给 `ashape()` 或 `ahull()` 会复用三角剖分，适合试验多个
alpha。

## 项目结构

```text
DynamicAlphaHull/
├── R/
│   ├── getDynamicAlphaHull.R  # 动态范围构建与离线海岸裁切
│   ├── delvor.R               # Delaunay 网格
│   ├── ashape.R               # alpha-shape 边筛选
│   ├── ahull.R                # alpha-hull 圆弧生成与裁剪
│   ├── complement.R           # 补集几何
│   ├── ah2terra.R             # 圆弧 hull 转 terra 多边形
│   └── geometry.R             # 圆、旋转、弧与多边形基础几何
├── data/ne_50m_land.rda       # lazy internal data：已融合的离线陆地图层
├── tests/testthat/             # 几何、Delaunay 与动态范围测试
├── man/                        # roxygen2 生成的帮助页面
├── DESCRIPTION                 # 包元数据和依赖
└── NAMESPACE                   # roxygen2 生成的导出声明
```

所有公开函数以及注册的 `plot.SpatVector()` 方法都有完整 roxygen2 文档和可直接运行的
最小示例。包级 roxygen 导入声明会生成 `NAMESPACE`，因此实现代码直接调用已导入的依赖，
不再散落 `包名::函数()` 写法。测试覆盖 Delaunay 构建、基础几何、动态范围、离线裁切和
退化输入；开发时先更新文档，再运行：

```r
roxygen2::roxygenise(".")
testthat::test_local(".")
```

## 依赖与边界

运行时仅依赖 `data.table`、`interp`、`purrr` 与 `terra`。`interp` 提供底层
Delaunay 三角剖分（与 `alphahull` 相同），`terra` 提供全部空间对象和 GIS 运算；本包
在此之上实现 alpha-shape、alpha-hull 与动态范围选择。

请在调用前完成分类学、坐标精度、海陆一致性和异常点清洗。范围多边形只是输入记录及
参数的几何概括，并不自动等同于物种真实分布或适生范围。
