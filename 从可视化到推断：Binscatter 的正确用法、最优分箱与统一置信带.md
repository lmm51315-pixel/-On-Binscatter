# 从可视化到推断：Binscatter 的正确用法、最优分箱与统一置信带

> **姓名：** 马俊豪 (中山大学)，郝泓程（中山大学）
>
> **E-mail:**  <majh28@mail2.sysu.edu.cn>

---

> Source: Cattaneo M D, et al.（2024）. On Binscatter. American Economic Review, 114(5): 1488-1514. —[Link](https://www.aeaweb.org/articles?id=10.1257/aer.20221576)
>
> 完成重现资料: [On Binscatter](https://github.com/lmm51315-pixel/-On-Binscatter)
> 其中放置了重现该文的数据、复现内容、原文和附录等相关参考资料

---

**目录**

[toc]

---


## 1. 文章内容介绍

随着经济学实证研究步入大数据时代，面对百万量级的观测样本，传统的散点图往往因数据过度重叠而失效，难以揭示变量间的潜在关系。分箱散点图（Binscatter）作为一种有效的降维可视化工具应运而生。其基本思想是将解释变量 $x$ 的值域划分为若干区间（箱），并计算每个箱内 $x$ 与因变量 $y$ 的均值，继而绘制这些均值点。该方法能清晰展现 $x$ 对 $y$ 的条件均值函数 $E[y|x]$ 的形态，主要服务于三个目的：(1) 非参数关系探索；(2) 对参数模型（如线性设定）进行初步诊断；(3) 在控制其他变量后，展示核心变量的“净”相关关系。

然而，Cattaneo等人（2024， *AER*）在‘On Binscatter'一文中指出，该工具在广泛应用的同时缺乏坚实的统计理论基础，且业界通行的“两步残差化”实践存在根本缺陷。该做法受Frisch-Waugh-Lovell定理启发，先分别将 $y$ 和 $x$ 对控制变量 $w$ 做线性回归得到残差，再对残差绘图。作者强调，此定理仅在模型完全线性的前提下成立。一旦真实关系包含非线性成分，此预处理步骤将扭曲 $x$ 与 $y$ 之间的关联图形，导致误导性结论。

为此，Cattaneo等人系统性地重构了Binscatter的方法论框架，其主要贡献包含三个递进的层面：

**第一，基于半参数模型的协变量调整。** 作者摒弃了两步残差法，提出了一个**联合估计**的半参数模型：$y_i = g(x_i) + w_i'\gamma + \epsilon_i$，其中 $g(\cdot)$ 为 $x$ 的非参数函数（通过分箱或样条逼近），$w_i'\gamma$ 为控制变量的线性参数部分。估计得到的函数 $\hat{g}(x)$ 代表了在控制变量处于样本均值水平时，$x$ 对 $y$ 的偏效应，从而在理论上一致地恢复了真实的函数形态与 $x$ 的原始尺度。

**第二，最优分箱数量的数据驱动选择。** 分箱数 $J$ 的选取长期依赖研究者主观判断。作者基于非参数估计中的偏差-方差权衡原理，推导了使积分均方误差最小化的最优分箱数 $J_{IMSE}$ 公式，并提供了数据驱动的选择算法，使这一关键参数的选择标准化、最优化。

**第三，构建一致置信带以支持统计推断。** 传统Binscatter仅提供点估计，无法量化不确定性或检验整体函数形态。作者利用非参数推断理论，构建了覆盖整个函数 $g(x)$ 的**一致置信带**，并采用稳健偏差校正技术处理估计偏差。这使得研究者能够正式检验诸如“函数是否为线性”的假设（若能找到一条直线完全落于置信带内，则不能拒绝线性原假设），将分析从视觉观察提升至严谨的统计推断。

本文后续结构如下：第二部分详述上述理论框架；第三部分介绍其对应的Stata命令 `binsreg`；第四部分通过代码实例进行演示。

## 2. 理论框架与方法比较

### 2.1 Binscatter 的基准构造与传统协变量处理

#### 2.1.1 Binscatter：作为非参数级数估计量的形式化

分仓散点图（binscatter）可被形式化为一种**非参数级数（series）估计量**。Cattaneo 等将其刻画为：在**数据驱动的分仓划分**之上，实施**零阶（piecewise-constant）分段多项式回归**所得到的回归函数估计。下文给出该估计量的标准构建过程。

设观测样本 $\{(y_i,x_i)\}_{i=1}^n$ 为独立同分布（i.i.d.）。研究目标是估计未知的条件均值函数
$$
v_0(x)=\mathbb{E}[\,y_i\mid x_i=x\,].
$$

##### 2.1.1.1 分仓划分

首先将解释变量 $x$ 的支撑集划分为 $J$ 个两两不交的区间（bins）。实践中通常采用**经验分位数划分**以使各仓样本量大致均衡。令 $x_{(i)}$ 表示样本的第 $i$ 个次序统计量，并定义分仓方案 $\Delta=\{\mathcal{B}_1,\mathcal{B}_2,\ldots,\mathcal{B}_J\}$。第 $j$ 个分仓可形式化为
$$
\hat{\mathcal{B}}_{j}=
\begin{cases}
\left[ x_{(1)},\, x_{(\lfloor n/J\rfloor)}\right] & \text{if } j=1,\\[4pt]
\left( x_{(\lfloor n(j-1)/J\rfloor)},\, x_{(\lfloor nj/J\rfloor)}\right] & \text{if } j=2,\ldots,J-1,\\[4pt]
\left( x_{(\lfloor n(J-1)/J\rfloor)},\, x_{(n)}\right] & \text{if } j=J.
\end{cases}
$$

##### 2.1.1.2 基函数表示

为将分仓结构转化为可估计的回归模型，在任意评估点 $x$ 处构造一组基函数（此处为分仓指示函数）。定义基向量
$$
b(x)=\Bigl(
\mathbb{1}\{x\in\hat{\mathcal{B}}_1\},\,
\mathbb{1}\{x\in\hat{\mathcal{B}}_2\},\,
\ldots,\,
\mathbb{1}\{x\in\hat{\mathcal{B}}_J\}
\Bigr)^{\top}.
$$
其中$\mathbb{1}\{\cdot\}$ 为示性函数：当括号内事件成立时取值为 1，否则为 0。由分仓区间互不重叠可知，对任意观测 $x_i$，向量 $b(x_i)$ 中**恰有一个分量为 1**，其余均为 0，从而实现对样本空间的正交划分。

##### 2.1.1.3. 零阶分段回归与 OLS 估计

Binscatter 对 $v_0(x)$ 的估计是一个**分段常数函数**（零阶分段多项式）。定义估计量
$$
\hat v(x)=b(x)'\hat\xi,
$$
其中系数向量 \(\hat\xi\) 由普通最小二乘法给出：
$$
\hat\xi=\arg\min_{\xi\in\mathbb{R}^J}\sum_{i=1}^n\bigl(y_i-b(x_i)'\xi\bigr)^2.
$$
由于 $b(x_i)$ 在每个观测上仅激活一个分仓指示，最小二乘解具有直接的样本均值解释：$\hat\xi_j$ 等于落入第 $j$ 个分仓的样本 $y_i$ 的均值。因此，图中呈现的 $J$ 个“散点”本质上对应于 $\hat v(x)$ 在各分仓区间内的常数取值。

在适当的正则条件下，随着样本量 $n$ 增大且分仓数 $J$ 以合适速率增长（以平衡逼近偏误与抽样方差），该级数估计量可一致收敛于真实的条件均值函数 $v_0(x)$。

#### 2.1.2 残差化协变量控制的局限性

在经典的参数化线性回归框架下，Frisch–Waugh–Lovell（FWL）定理表明：在模型
$$
y_i = x_i\beta + w_i'\gamma + \varepsilon_i
$$
中，$\beta$ 的 OLS 估计等价于先分别将 $y_i$ 与 $x_i$ 对 $w_i$ 做线性回归取残差，再用残差回归。这一等价性依赖于目标参数本身是**线性系数**，以及模型在 $x$ 上的**线性结构**。

受此启发，实务中常见的“残差化 Binscatter”操作是：

1) 定义线性投影算子 $L(\cdot\mid w)$，表示变量 $a$ 在 $w$ 张成空间上的最佳线性逼近：
   $$
   L(a\mid w)=\Pi_{(1,w)}a.
   $$

2) 构造残差：
   $$
   \tilde y_i = y_i - L(y\mid w_i),\qquad
   \tilde x_i = x_i - L(x\mid w_i).
   $$

3) 将 $(\tilde y_i,\tilde x_i)$ 代入无协变量的 binscatter 估计量进行分仓绘图。

Cattaneo 等人指出，上述做法在大样本下收敛于如下“残差—残差”的条件期望对象：
$$
\plim\,\hat v_{\text{resid}}(\cdot)
=\mathbb{E}\!\left[\,y_i-L(y\mid w_i)\ \big|\ x_i-L(x\mid w_i)\,\right].
$$
该极限一般**并不等于**我们在半线性/部分线性设定下希望识别的条件均值曲线 $\mu_0(x)$。

为说明偏差来源，考虑真实的部分线性模型：
$$
y_i=\mu_0(x_i)+w_i'\gamma_0+\varepsilon_i.
$$
利用线性投影的线性性质，有
$$
L(y\mid w_i)=L(\mu_0(x_i)\mid w_i)+w_i'\gamma_0,
$$
从而 $y$ 的残差可写为
$$
\tilde y_i
= y_i-L(y\mid w_i)
= \mu_0(x_i)-L(\mu_0(x_i)\mid w_i)+\varepsilon_i.
$$
可见，残差化确实“剔除”了线性控制项 $w_i'\gamma_0$，但同时也引入了
$\mu_0(x_i)$ 相对于 $w_i$ 的线性投影项 $L(\mu_0(x_i)\mid w_i)$。因此，残差化 binscatter 实际可视化的是：

- 横轴：$\tilde x_i = x_i - L(x\mid w_i)$；
- 纵轴：$\tilde y_i \approx \mu_0(x_i)-L(\mu_0(x_i)\mid w_i)$（忽略噪声项）。

这是否等价于 $\mu_0(x)$，关键取决于 $\mu_0(\cdot)$ 的函数形式：

##### 2.2.2.1  线性情形：残差化与目标函数一致

若 $\mu_0(x)=x\beta$，则
  
  $$
  L(\mu_0(x_i)\mid w_i)=L(x_i\beta\mid w_i)=\beta\,L(x_i\mid w_i),
  $$

  因而

  $$
  \tilde y_i = \beta\bigl(x_i-L(x\mid w_i)\bigr)+\varepsilon_i=\beta\,\tilde x_i+\varepsilon_i,
  $$
  残差化方法可正确恢复线性斜率 $\beta$。这正是 FWL 在 OLS 中成立的根本原因。

##### 2.2.2.2 非线性情形：残差化导致目标曲线扭曲

若 $\mu_0(\cdot)$ 非线性，则一般不存在
  $$
  L(\mu_0(x_i)\mid w_i)=g\!\left(L(x_i\mid w_i)\right)
  $$
  之类的“可交换”关系；更直观地说，线性投影算子无法穿透非线性变换（例如通常有 $L(x_i^2\mid w_i)\neq L(x_i\mid w_i)^2$）。因此，$\mu_0(x)$ 的形状会被项 $L(\mu_0(x_i)\mid w_i)$ 扭曲，导致残差化 Binscatter 在一般情形下并不描绘目标函数 $\mu_0(x)$，从而可能产生误导性的图形结论。

#### 2.2.3 改进方案：部分线性框架下的一步估计

为克服“残差化（residualization）”在模型设定偏误下可能导致的不一致性，Cattaneo 等提出在**半参数部分线性模型**框架下直接构造协变量调整的 Binscatter 估计量。其核心假定为：
$$
\mathbb{E}[y_i \mid x_i, w_i] \;=\; \mu_0(x_i) + w_i' \gamma_0,
$$
其中，$\mu_0(x)$ 是关于核心解释变量 $x$ 的未知（可高度非线性）函数，$w_i'\gamma_0$ 则刻画控制变量 $w_i$ 的线性影响。该设定**不要求**对 $x$ 的效应施加线性限制，从而避免残差化方法隐含的线性结构依赖。

在估计上，方法采用“一步法”联合估计非参数部分与参数部分：令 $b(x_i)$ 为对 $x$ 的分仓/基函数向量，则
$$
(\hat\beta,\hat\gamma)
\;=\;
\arg\min_{\beta,\gamma}
\sum_{i=1}^n
\bigl(y_i - b(x_i)'\beta - w_i'\gamma\bigr)^2.
$$
其中，$\beta$ 捕捉 $x$ 的非线性形状（由分仓基函数逼近），$\gamma$ 捕捉协变量的线性效应；在上述半参数 DGP 下，联合最小二乘可保证 $\hat\beta,\hat\gamma$ 的一致性。

为可视化在“控制其他协变量”的条件下 $x$ 与 $y$ 的偏均值关系，通常将 $w$ 固定在某一代表性水平（例如样本均值 $\bar w$），从而得到用于绘图的估计函数：
$$
\hat\Upsilon(x)
\;=\;
b(x)'\hat\beta + \bar w'\hat\gamma.
$$
该函数刻画了在协变量取均值时，$x$ 对 $y$ 的条件期望（偏均值）关系。

**方法优势：**

1) 允许 $\mu_0(x)$ 任意非线性，并能在 $x$ 与 $w$ 存在复杂相关性时稳健恢复 $\mu_0(x)$ 的形状；
2) 直接在原始尺度上拟合与作图，保留 $x$ 的取值范围与经济含义，避免残差化导致的横轴压缩与解释困难。

### 2.2 最优分箱数 $J$ 的选择原则与实现

根据前文分析，我们已经知道$J$ 过小则分箱过宽、估计更平滑但逼近偏差上升，难以刻画 $\mu_0(x)$ 的局部非线性；$J$ 过大则分箱过窄、逼近偏差下降但单箱样本量不足，抽样方差上升并导致图形噪声显著。

为在全局意义下权衡上述取舍，作者以积分均方误差（IMSE）作为目标函数，选择使估计函数 $\Upsilon(x)$ 与真值 $\Upsilon_0(x)$ 的总体误差最小的 $J$：
$$
\text{IMSE}(J)=\mathbb{E}\!\left[\int\bigl(\Upsilon(x)-\Upsilon_0(x)\bigr)^2 f_X(x)\,dx\right].
$$
在对 IMSE 的主导项展开并最小化后，可得到最优分箱数的闭式近似解：
$$
J_{\text{IMSE}}\approx \mathcal{B}_n\,\mathcal{V}_n^{-1/3}\,n^{1/3},
$$
其中 $\mathcal{B}_n$ 为与目标函数曲率（如二阶导数）相关的偏差常数，$\mathcal{V}_n$为与噪声强度相关的方差常数，$n$ 为样本量。该结果表明 $J$ 不应固定：样本量增大时最优 $J$ 随 $n^{1/3}$ 增长；在噪声更强时倾向于取更小的 $J$ 以抑制方差，在函数更“弯曲”或更非线性时则应取更大的 $J$ 以降低逼近偏差。

### 2.3 可视化不确定性：RBC 一致置信带

在给定 $J_{\text{IMSE}}$ 并得到点估计 $\hat\Upsilon(x)$ 后，需要判断图形形状是否具有统计显著性。为此，Cattaneo 等提出对目标函数 $\Upsilon_0(x)$ 构造**一致置信带**（UCB），以在全域 $x\in\mathcal{X}$ 上同时满足覆盖率：
$$
\lim_{n\to\infty}\Pr\bigl(\Upsilon_0(x)\in I(x),\ \forall x\in\mathcal{X}\bigr)=1-\alpha.
$$
一致置信带不同于逐点误差棒：它控制的是“整条曲线”的同时覆盖，从而避免逐点区间的多重比较失真。

难点在于：$J_{\text{IMSE}}$ 下偏差与方差同阶，直接据此构造置信带会因偏差而欠覆盖。作者采用**稳健偏差校正**（RBC）：先估计并扣除主导偏差，再将偏差估计带来的额外不确定性纳入方差。最终置信带为
$$
I_{\text{RBC}}(x)=\hat\Upsilon_{\text{BC}}(x)\ \pm\ c_{\text{RBC}}\sqrt{\hat\Omega_{\text{RBC}}(x)},
$$
其中 $\hat\Upsilon_{\text{BC}}(x)$ 为去偏中心线，$\hat\Omega_{\text{RBC}}(x)$ 为稳健方差，$c_{\text{RBC}}$ 为经极值型校准（常用模拟）得到的全域临界值。

该置信带可用于整体推断：若某条直线整体落入带内，则无法拒绝线性；若置信带呈系统性弯曲并排除任何直线，则支持显著非线性；若上下界随 \(x\) 同向变化，则为单调性提供证据。

## 3. 命令介绍

为便于不同研究者在各自的计算环境中复现与扩展，作者团队提供了 `binsreg` 在 Python、R 与 Stata 等平台上的实现版本，具体信息可参见官方代码仓库：-[binsreg](https://github.com/nppackages/binsreg)。

本文的实证复现与结果展示以 **Stata** 实现为准。

### 3.1 命令安装

**net install**

```stata
net install binsreg, from("https://raw.githubusercontent.com/nppackages/binsreg/master/stata") 
```

### 3.2 命令语法

在 Stata 中可通过 `net install` 从作者维护的在线安装源获取并安装 `binsreg`：

```stata
binsreg depvar indvar [othercovs] [if] [in] [weight] [ , deriv(v)
            at(position)
            absorb(absvars) reghdfeopt(reghdfe_option)
            dots(dotsopt) dotsgrid(dotsgridoption) dotsplotopt(dotsoption)
            line(lineopt) linegrid(#) lineplotopt(lineoption)
            ci(ciopt) cigrid(cigridoption) ciplotopt(rcapoption)
            cb(cbopt) cbgrid(#) cbplotopt(rareaoption)
            polyreg(p) polyreggrid(#) polyregcigrid(#)
            polyregplotopt(lineoption)
            by(varname) bycolors(colorstylelist) bysymbols(symbolstylelist)
            bylpatterns(linepatternstylelist)
            nbins(nbinsopt) binspos(position) binsmethod(method) nbinsrot(#)
            samebinsby randcut(#)
            pselect(numlist) sselect(numlist)
            nsims(#) simsgrid(#) simsseed(seed)
            dfcheck(n1 n2) masspoints(masspointsoption)
            vce(vcetype) asyvar(on/off)
            level(level) usegtools(on/off) noplot savedata(filename) replace
            plotxrange(min max) plotyrange(min max) twoway_options ]
```

部分重点选项：

- `absorb(absvars)`：吸收（高维）固定效应，用于控制不可观测异质性（如地区、年份固定效应）。
- `vce(vcetype)`：指定方差估计方式（常用 `vce(robust)` 或 `vce(cluster clustvar)`；也可多维聚类）。
- `at(position)`：设定协变量取值用于作图/预测（最常用 `at(mean)`，即把控制变量固定在样本均值处画偏均值曲线）。
- `nbins(nbinsopt)`：指定分箱数 \(J\)（可手动给定或用于覆盖默认选择）。
- `binsmethod(method)`：指定分箱数选择方法（用于数据驱动选择 \(J\)，与偏差—方差权衡相关）。
- `dots(dotsopt)`：绘制分箱散点（binned dots）及其构造方式（图中的“点”）。
- `line(lineopt)`：绘制分箱拟合线/分段多项式拟合（图中的“主曲线”）。
- `ci(ciopt)`：绘制逐点置信区间（pointwise CI，误差棒/区间）。
- `cb(cbopt)`：绘制一致置信带（uniform confidence band，用于对整条曲线进行全域推断）。
- `nsims(#)`：用于置信带（尤其 `cb()`）的模拟次数，次数越多越稳定但更耗时。
- `simsseed(seed)`：设置模拟随机种子，保证置信带结果可复现。
- `by(varname)`：分组作图/估计（对不同组分别画曲线用于比较）。
- `samebinsby`：在 `by()` 分组下强制使用相同分箱切点，提升组间可比性。
- `polyreg(p)`：叠加 \(p\) 阶全局多项式回归作为对照线（便于对比线性/非线性形状）。

## 4. Stata 实操：二十世纪税收与创新的关系（Akcigit et al., 2022. QJE）

本节旨在系统复现 Akcigit, Grigsby, Nicholas, and Stantcheva（2022）发表于 *The Quarterly Journal of Economics* 的论文《Taxation and Innovation in the Twentieth Century》（下文简称 AGNS 2022）。该研究的核心旨在识别并量化税收激励对创新产出的因果效应，其中心研究问题为：边际税率的变化是否会对企业及个人的创新活动产生显著的抑制或激励作用。

为验证这一假说，AGNS（2022）构建了以下核心变量：

- **被解释变量** \( y \) ：以各州年度专利授予数量的自然对数度量创新产出；
- **核心解释变量** \( x \) ：前90%收入阶层的边际净税率，定义为 \( 1 - \tau \) 的自然对数，其中 \( \tau \) 为边际税率。

在识别策略上，为缓解遗漏变量偏误，基准模型纳入了包括滞后的企业税率、人口密度、人均实际GDP、研发税收抵免强度等一系列时变州级经济特征作为控制变量。此外，模型通过引入州固定效应与年份固定效应，以控制不随时间变化的州际异质性与全国性的时间趋势，使控制变量的总维度扩展至约 \( d = 113 \)。

### 4.1 方法的比较：残差化与协变量调整

在高维固定效应模型设定下，若试图对条件期望函数 \( \mathbb{E}[y \mid x, w] \) 进行完全非参数估计，将面临严重的“维数诅咒”问题，导致估计精度急剧下降。因此，实践中通常采用半参数模型中的部分线性框架进行可视化与推断，即在非参数地拟合 \( x \) 与 \( y \) 关系的同时，以参数形式控制高维协变量 \( w \) 的影响。

传统的“两阶段残差化”方法（Two-Stage Residualization）是达成此目的的常用策略。然而，如下图所示，该方法在处理高维固定效应时存在显著缺陷。

```stata
*----- 图像设置 -----*
clear all
set scheme s2color
*----- 基础设置 -----*
global main "/Users/lmm/Documents/Binscatter论文复现/"
global data "$main/02_data_raw"
global temp "$main/03_temp"
global output "$main/04_output"
*----- 设置工作路径&读取数据 -----*
cd "$main"
use "$data/CCFF_2024_AER--AGNS.dta", clear

*-------------------- A：直接使用残差化方法绘图 --------------------*
*-------------------- A1. 残差化因变量 lnpat --------------------*
areg lnpat top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 i.statenum ///
    [aw=pop1940], a(year)
predict resid_lnpat, res
egen mean_lnpat = mean(lnpat)
replace resid_lnpat = resid_lnpat + mean_lnpat  // "均值校准"：将残差平移回原始量级附近
*-------------------- A2. 残差化核心解释变量 mtr90_lag3 --------------------*
areg mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 i.statenum ///
    [aw=pop1940], a(year)
predict resid_mtr90_lag3, res
egen mean_mtr90_lag3 = mean(mtr90_lag3)
replace resid_mtr90_lag3 = resid_mtr90_lag3 + mean_mtr90_lag3   // 同样做均值校准，使横轴变量回到原量级附近
*-------------------- A3. 残差化 binscatter：叠加线性对照线 --------------------*
binsreg resid_lnpat resid_mtr90_lag3 [aw=pop1940], ///
    nbins(50) polyreg(1) ///
    dotsplotopt(mcolor(dkorange)) ///
    polyregplotopt(lcolor(black)) ///
    savedata($temp/tmpAGNS1binscatter) replace ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(5.7(0.1)6.2, nogrid) xlabel(-0.525(0.05)-0.325)
graph export "$output/AGNS_binscatter.png", as(png) replace width(800) height(600)

*-------------------- B. 修正横轴范围后的残差化方法 --------------------*
*-------------------- B1. 生成协变量调整 binsreg 的 savedata --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], ///
    nbins(50) absorb(statenum year) at(mean) polyreg(1) ///
    savedata($temp/tmpAGNS1binsreg) replace
preserve
*-------------------- B2. 合并两套 savedata，并标记来源 --------------------*
use "$temp/tmpAGNS1binscatter", clear
gen binsreg = 0
append using "$temp/tmpAGNS1binsreg"
replace binsreg = 1 if missing(binsreg)  // 1：协变量调整 binsreg 数据（append 进来时 binsreg 变量缺失）
*-------------------- B3. 正确的表示原来的比例尺（虚线） --------------------*
local new = _N + 2
set obs `new'    // 扩充两行观测作为线段的端点
replace poly_x = -0.8 if _n == `new'-1
replace poly_x =  0   if _n == `new'   // 使横轴范围与协变量调整图的横轴范围一致
replace binsreg = 2 if _n >= `new'-1   
gen poly_fit2 = .
reg poly_fit poly_x if binsreg==0 // 用残差化方法的 poly_fit 与 poly_x 拟合线性关系
replace poly_fit2 = _b[_cons] + _b[poly_x]*poly_x if binsreg==2
*-------------------- B4. 最终绘图 --------------------*
tw ///
    (scatter dots_fit dots_x if binsreg==0, color(dkorange)) ///
    (line poly_fit poly_x if binsreg==0, color(black) lwidth(medthick)) ///
    (line poly_fit2 poly_x if binsreg==2, color(black) lpattern(dot) lwidth(medthick)), ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(5(.5)8) legend(off)
graph export "$output/AGNS_covariateAdjustments_binscatter.png", replace as(png) width(800) height(600)

restore
```

<div style="display:flex; gap:16px; align-items:flex-start; margin-bottom: 20px;">
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.1.1 (A)：残差化图示（失真尺度）</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_binscatter.png" alt="图4.1.1 (A)：残差化图示（失真尺度）" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：横轴为残差化处理后的解释变量，其经济含义的原始尺度已丢失。</div>
  </div>
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.1.1 (B)：残差化图示（校正尺度）</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_covariateAdjustments_binscatter.png" alt="图4.1.1 (B)：残差化图示（校正尺度）" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：横轴虽已校正回原始单位，但有效变异范围因残差化过程被严重压缩。</div>
  </div>
</div>

图4.1.1 表明，残差化方法会同时扭曲估计函数关系的**几何形态**与解释变量的**有效定义域**。具体而言，当模型包含高维固定效应时，核心解释变量 \( x \) 的大部分变异被吸收，导致其残差 \( \tilde{x} \) 的取值范围相对于原始经济区间 \([-0.9, 0]\) 被极端压缩。其结果是：第一，可视化图形的横轴尺度失去直接的经济解释力；第二，即便将坐标轴重新标定为原始单位，由于 \( \tilde{x} \) 的变异主要源于噪声，图示也难以提供关于条件均值函数形态的有效信息，点云呈现近似随机散布。

作为对比，采用 **协变量调整回归法（Covariate-Adjusted Regression）** 与分箱估计量，可以直接在原始解释变量尺度上进行估计与可视化，从而避免了上述问题。

```stata
*-------------------- C 协变量调整 binsreg（不含线性对照线） --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], ///
    nbins(50) ///
    absorb(statenum year) replace at(mean) ///
    dotsplotopt(mcolor(black)) ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量") xtitle("前 90% 收入者的边际净税率") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0)
graph export "$output/AGNS_covariateAdjustments_binsreg_woutRegLine.png", replace as(png) width(800) height(600)


```

<div style="margin-bottom: 20px;">
  <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.1.2：协变量调整分箱估计（Covariate-Adjusted Binscatter）</div>
  <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_covariateAdjustments_binsreg_woutRegLine.png" alt="图4.1.2：协变量调整分箱估计" style="width:100%; max-width:600px; height:auto; display:block; margin: 0 auto; border: 1px solid #ddd;">
  <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：横轴完整保留了边际净税率的原始经济区间 \([-0.9, 0]\)，分箱局部均值点清晰地呈现出线性趋势。</div>
</div>

如图4.1.2所示，在恢复至具有明确经济含义的横轴尺度后，分箱局部均值点呈现出高度线性的排列模式。这一结果与AGNS（2022）基准回归中采用的线性设定高度一致。两者的对比表明，传统残差化图示中观察到的“非线性”与“数据散乱”很可能是该数据处理步骤所引入的**统计假象（Statistical Artifact）**；而协变量调整后的可视化则为此项研究采用线性参数模型提供了更为可靠和直观的图形证据。

### 4.2 最优分箱数$(J）$的统计意义与选择

如理论部分所述，在非参数或半参数的分箱估计中，分箱数量 $J$ 是一个关键的平滑参数。其选择并非任意，而是在偏差与方差之间进行权衡的根本性决策，直接影响估计量的有限样本性质及其可视化呈现的可靠性。不同 \(J\) 值下估计结果的差异可通过以下 Stata 代码生成的图示进行对比分析。

```stata
use $data/CCFF_2024_AER--AGNS, clear

*-------------------- 场景A：分箱数过少(J=5，欠平滑） --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], absorb(statenum year) nbins(5) at(mean) polyreg(1) ///
    dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量（对数）") xtitle("前90%收入者的边际净税率（对数）") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9)
graph export "$output/AGNS_nbins5.png", replace as(png) width(800) height(600)

*-------------------- 场景B：基于积分均方误差（IMSE）的最优分箱 --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], absorb(statenum year) vce(cluster fiveyrblockbystate year) ///
    randcut(1) at(mean) polyreg(1) ///
    dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量（对数）") xtitle("前90%收入者的边际净税率（对数）") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9)
graph export "$output/AGNS_nbinsOptimal.png", replace as(png) width(800) height(600)

*-------------------- 场景C：分箱数过多（J=50，过平滑） --------------------*
binsreg lnpat mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
    [aw=pop1940], absorb(statenum year) nbins(50) at(mean) polyreg(1) ///
    dotsplotopt(mcolor(black)) polyregplotopt(lcolor(forest_green)) ///
    graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
    ytitle("专利数量（对数）") xtitle("前90%收入者的边际净税率（对数）") ///
    ylabel(, nogrid) xlabel(-0.9(0.3)0) ylabel(6(1)9)
graph export "$output/AGNS_nbins50.png", replace as(png) width(800) height(600)
```

上述代码生成了三种不同平滑程度下的条件均值函数估计图，其对比结果如下：

<div style="display:flex; gap:16px; align-items:flex-start; margin-bottom: 20px;"> 
  <div style="flex:1;"> 
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.2.1 (A)：分箱数过少 (J=5))</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_nbins5.png" alt="图4.2.1 (A)：分箱数过少" style="width:100%; height:auto; border: 1px solid #ddd;"> 
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">欠平滑状态：箱体过宽导致局部常数近似失效，估计偏差主导。
    </div> 
  </div>
  <div style="flex:1;"> 
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.2.1 (B)：IMSE最优分箱</div> 
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_nbinsOptimal.png" alt="图4.2.1 (B)：IMSE最优分箱" style="width:100%; height:auto; border: 1px solid #ddd;"> 
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">偏差-方差权衡：最优的 J，平衡平滑性与保真度。
    </div> 
  </div> 
  <div style="flex:1;"> 
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.2.1 (C)：分箱数过多 (J=50)</div> 
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_nbins50.png" alt="图4.2.1 (C)：分箱数过多" style="width:100%; height:auto; border: 1px solid #ddd;"> 
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">过平滑状态：箱体过窄导致样本稀疏，估计方差急剧增大。
    </div> 
  </div> 
</div>

通过对比可以得出以下统计见解：

- **分箱过少 (欠平滑)**：如图4.2.1(A)所示，当 $J$ 取值过小时，每个箱体的区间过宽，“箱内条件均值为常数”（或低阶多项式）的局部近似过于粗糙。这导致估计量的**近似偏误**显著增大，可能过度平滑或“抹平”真实条件均值函数中潜在的细微变化或非线性结构，使图形呈现欠拟合状态。
- **分箱过多 (过平滑)**：如图4.2.1(C)所示，当 $J $取值过大时，每个箱体包含的观测值数量锐减。这虽然减少了模型设定带来的偏误，但导致估计量的**抽样方差**急剧增大，表现为图形出现剧烈的“抖动”或不稳定性。这种不稳定性并非源于真实的数据生成过程，而是估计精度不足所致，可能误导研究者识别出不存在的局部波动。
- **最优分箱选择**：如图4.2.1(B)所示，采用基于积分均方误差（Integrated Mean Squared Error, IMSE）准则数据驱动选择的$J$，能够在偏误与方差之间实现最优权衡。该准则旨在最小化估计函数与真实函数之间差异的整体期望（积分意义下的均方误差），从而使得到的非参数估计量既能有足够的灵活性去捕捉潜在的函数形态，又能保持必要的稳定性以进行可靠的统计推断。

因此，选择适当的分箱数 $J$ 并非一个简单的可视化问题，而是非参数估计理论中的核心步骤。一个经严谨准则（如IMSE）选择的$J$，是确保可视化结果能够作为有效统计证据、准确揭示变量间潜在数据结构的前提。

### 4.3 基于置信带的统计推断

需要再次强调的是，分箱散点图（binscatter）中所呈现的每一个“点”，本质上并非原始观测值的直接绘制，而是对潜在条件期望函数 \( \mathbb{E}[y \mid x, w] \)（或其经过协变量调整后的形式）的一种非参数局部估计量。因此，对图形所展示关系之不确定性的量化，应围绕**估计函数本身的抽样分布**进行，而非局限于单一参数的显著性检验。

为实现对函数形态的有效推断，关键在于选择适当的分箱数量 \( J \)。其选择标准需满足双重目标：一方面，\( J \) 必须随样本量 \( n \) 的增加而增长，以渐近消除因分箱离散化而产生的“模型设定偏误”（specification bias）或“近似误差”；另一方面，\( J \) 的增长速率又须受到约束，以确保每个箱体内拥有足够的观测值来控制估计量的方差。简言之，最优的 \( J \) 应在偏差与方差之间寻求权衡，从而使分箱估计量能够一致地逼近真实的条件均值函数。下图展示了不同 \( J \) 选择策略下估计结果的差异：

```stata
*-------------------- 固定分箱数 J=5 的 binscatter + 置信带--------------------*
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 /// 
[aw=pop1940],  ///
randcut(1) absorb(statenum year) nbins(5) ///
at(mean) cb(0 0) dotsplotopt(mcolor(black)) ///
polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) ///
plotregion(lcolor(black)) ytitle("专利数量") ///
xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/AGNS_fixedJconfidenceBandwDots.png", replace as(png) width(800) height(600)

*-------------------- 采用数据驱动的J的 binscatter + 置信带 --------------------*
use "$data/CCFF_2024_AER--AGNS.dta", clear
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
 [aw=pop1940], ///
vce(cluster fiveyrblockbystate year) absorb(statenum year) ///
randcut(1) cb(T) replace at(mean) dotsplotopt(mcolor(black)) ///
polyregplotopt(lcolor(forest_green)) graphregion(color(white) margin(large)) ///
plotregion(lcolor(black)) ytitle("专利数量") ///
xtitle("前 90% 收入者的边际净税率") ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/AGNS_confidenceBandwDots.png", as(png) replace 
```

<div style="display:flex; gap:16px; align-items:flex-start; margin-bottom: 20px;">
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.3.1：固定小J下的置信带</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_fixedJconfidenceBandwDots.png" alt="图4.3.1：固定小J下的置信带" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：当 J  固定且取值过小时，估计更类似于“分组比较”，因过宽的箱体导致显著的平滑偏误，可能掩盖真实的函数形态。</div>
  </div>
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:13px; text-align: center;">图4.3.2：数据驱动选择J下的置信带</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_confidenceBandwDots.png" alt="图4.3.2：数据驱动选择J下的置信带" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：当  J  根据数据（如依均方误差（IMSE）准则）选择时，估计能更灵活地捕捉数据的潜在结构，为统计推断提供有效基础。</div>
  </div>
</div>

如图4.3.1所示，当 \( J \) 取值过小且固定不变时，估计过程实质上退化为粗糙的“分组比较”。过宽的箱体导致估计量的偏差增大，以至于从图形上难以辨别变量间的系统性关系。例如，可以构造一条水平参考线，使其完全位于各个箱体估计点对应的置信区间之内。这一现象在统计上意味着，我们无法拒绝“不同税率区间内专利数量的条件均值无显著差异”的原假设。

反之，如图4.3.2所示，当 \( J \) 根据数据驱动准则（如前文推导的基于积分均方误差（IMSE）的最优分箱数 \( J_{\text{IMSE}} \)）选取时，估计量能够更精细地反映数据的局部特征。此时，围绕估计函数构建的置信带为图形化统计推断提供了严谨的框架。具体而言，我们可以利用置信带检验关于函数形态的特定假设。例如，若要检验“条件均值函数是否为线性”这一原假设，可观察是否存在一条完整的线性函数完全位于置信带内部。下图演示了这一推断过程：

```stata

*-------------------- 置信带（confidence band）及辅助直线/水平线展示--------------------*
use $data/CCFF_2024_AER--AGNS, clear

*-------------------------- 用 binsreg 生成“同时置信带” ------------------------------*
binsreg lnpat  mtr90_lag3 top_corp_lag3 lreal_gdp_pc lpopulation_density rd_credit_lag3 ///
      [aw=pop1940], ///
      vce(cluster fiveyrblockbystate year) ///
      absorb(statenum year) ///
      randcut(1) ///
      cb(T) ///
      at(mean) ///
      nodraw ///
      savedata($temp/tmp3) replace ///
      nsims(5000) simsgrid(200)

*-------------------- 载入 binsreg 保存的“作图数据” --------------------*
use $temp/tmp3, clear

*-------------------- 构造一条“参考水平线” --------------------*
*    做法：在左侧区间（CB_x < -0.8）里，取置信带上界 CB_r 的最大值作为常数 y
*    解释：用于在图中直观展示某条水平线与置信带/分箱点的相对位置（辅助形状比较）
qui summ CB_r if CB_x < -0.8
gen horline = `r(max)'
*-------------------- 构造一条“线性趋势线” --------------------*
gen CB_avg = 0.4*CB_r + 0.6*CB_l
reg CB_avg CB_x
predict inline, xb 
*-------------------- 图1：置信带（阴影）+ 分箱点（黑点）+ 参考水平线（黑色虚线） 
tw ///
 (rarea CB_l CB_r CB_x, fcolor(edkblue) fintensity(20) lwidth(none)) ///
 (scatter dots_fit dots_x, mcolor(black)) ///
 (line horline CB_x, lcolor(black) lpattern(dash)), ///
 graphregion(color(white) margin(large)) ///
 plotregion(lcolor(black)) ///
 ylabel(, nogrid nogextend) ///
 legend(off) ///
 ytitle("专利数量") ///
 xtitle("前 90% 收入者的边际净税率")
graph export "$output/AGNS_confidenceBandwHorLine.png", replace as(png) width(800) height(600)
*-------------------- 图2：置信带（阴影）+ 参考直线（红色）） --------------------*
tw ///
 (rarea CB_l CB_r CB_x, fcolor(edkblue) fintensity(20) lwidth(none)) ///
 (line inline CB_x, lcolor(red)), ///
 graphregion(color(white) margin(large)) ///
 plotregion(lcolor(black)) ///
 ylabel(, nogrid nogextend) ///
 legend(off) ///
 ytitle("专利数量") ///
 xtitle("前 90% 收入者的边际净税率")
graph export "$output/AGNS_confidenceBandwLine.png", replace as(png) width(800) height(600)
```

<div style="display:flex; gap:16px; align-items:flex-start;">
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.3.3：对水平关系的检验</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_confidenceBandwHorLine.png" alt="图4.3.3：对水平关系的检验" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：一条水平线被置于置信带中，用于检验“无关系”的原假设。</div>
  </div>
  <div style="flex:1;">
    <div style="font-weight:600; margin-bottom:8px; text-align: center;">图4.3.4：对线性关系的检验</div>
    <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedAGNS_confidenceBandwLine.png" alt="图4.3.4：对线性关系的检验" style="width:100%; height:auto; border: 1px solid #ddd;">
    <div style="text-align: center; font-size: 0.9em; color: #555; margin-top: 5px;">注：一条倾斜的直线被置于置信带中，用于检验“线性关系”的原假设。</div>
  </div>
</div>

如图4.3.4所示，可以构造一条完整的线性函数，使其完全位于由数据驱动选择 \( J \) 所构建的置信带内部。根据置信带的定义，若原假设（函数为线性）成立，则构造的线性函数有较高概率被包含在内。因此，基于该图示证据，**我们不能拒绝“边际净税率与专利数量对数之间存在线性关系”的原假设**。这一图形化检验方法，为AGNS（2022）研究中采用线性回归设定提供了直观且严谨的辅助证据。

## 5. Stata 实操2：顶尖发明家生产率与高科技集群规模的关系 (Moretti et al., 2021, *AER*)

本节应用改进的Binscatter框架，对Moretti等人（2021）关于高科技集群规模对顶尖发明家创新产出影响的研究进行重新评估，以演示一个完整的图形统计推断流程。

Moretti等人的原始研究利用Binscatter初步揭示了集群规模可能存在边际效应递减的模式，但其图形构建依赖于简易的残差化方法。本研究将此案例作为检验Binscatter方法正确实施的第二个实证场景。为保持行文简洁，本节仅展示核心分析代码（对应图5.5与图5.6），完整复制代码详见项目仓库 [On Binscatter](https://github.com/lmm51315-pixel/)。下文列出了生成关键结果图的Stata代码示例：

```stata
use $data/CCFF_2024_AER--M2, clear
sort x
*-------------------- 图5.5：置信带  --------------------*
binsreg y x, absorb(year bea zd) vce(cluster cluster1) ///
randcut(1) cb(T) cbplotopt(fcolor(edkblue) ///
fintensity(20) lwidth(none)) dotsplotopt(mcolor(black)) ///
graphregion(color(white) margin(large)) ///
plotregion(lcolor(black)) ytitle("顶尖发明家创新产出") ///
xtitle("高科技集群规模") ///
ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/M_confidenceBand.png", replace  as(png) width(800) height(600)

*-------------------- 图5.6：置信带（高纬固定效应） --------------------*
binsreg y x, absorb(year bea zd class cluster1 cluster_bea_class cluster_zd_year cluster_class_year ///
inventor cluster_bea_year org_new) vce(cluster cluster1) ///
randcut(1) cb(T) cbplotopt(fcolor(edkblue) fintensity(20) ///
lwidth(none)) dotsplotopt(mcolor(black)) ///
graphregion(color(white) margin(large)) plotregion(lcolor(black)) ///
ytitle("顶尖发明家创新产出") ///
xtitle("高科技集群规模") ylabel(, nogrid) nsims(5000) simsgrid(200)
graph export "$output/M_confidenceBandFullSpec.png", replace as(png) width(800) height(600)
```

<div style="display:grid; grid-template-columns: 1fr 1fr; gap:16px; align-items:start;"> <div> <p><strong>(图5.1) 原始散点图</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_scatter.png" alt="原始散点图" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：每个点代表一个城市-技术领域-年份单元的人均专利产出与集群规模。图形显示总体正相关，但关系呈现显著非线性，尤其是在大规模集群处。</em></p> </div> <div> <p><strong>(图5.2) 论文再现 (Moretti, 2021)</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_binscatter.png" alt="论文再现图" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：控制年份、领域、城市固定效应后的残差化Binscatter及线性拟合。线性拟合斜率为正，但Binscatter估计点（橙色）提示潜在非线性。</em></p> </div> <div> <p><strong>(图5.3) 有误的残差化Binscatter</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_binscatter_pwc.png" alt="有误残差化结果" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：采用不恰当的残差化方法得到的图形扭曲且不规则，难以提供可靠的经济学解释。</em></p> </div> <div> <p><strong>(图5.4) 正确的协变量调整对比</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_covariateAdjustments.png" alt="协变量调整对比" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：黑色点为本文半参数框架下的Binscatter估计，橙色点为(c)中有误方法的结果。正确方法下的关系曲线更为平滑，呈现先上升后趋缓的形态。</em></p> </div> <div> <p><strong>(图5.5) 条件均值函数的置信带 (基础设定)</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_confidenceBand.png" alt="基础设定置信带" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：基于半参数Binscatter绘制的95%均匀置信带。阴影区域显示，在中小规模集群处条件均值函数斜率较陡，随后渐趋平缓，形态显著偏离直线。</em></p> </div> <div> <p><strong>(图5.6) 条件均值函数的置信带 (完整高维固定效应)</strong></p> <img src="https://fig-lianxh.oss-cn-shenzhen.aliyuncs.com/undefinedM_confidenceBandFullSpec.png" alt="完整设定置信带" style="width:100%; height:auto; border:1px solid #ddd;"> <p><em>注：控制原文最严格设定中的所有高维固定效应后，计算得到的Binscatter及置信带。非线性形态依然稳健存在。</em></p> </div> </div>

**综合发现与推断**：图5系列结果表明，在采用正确的协变量调整方法后，高科技集群规模对顶尖发明家人均产出的边际效应并非恒定，而是呈现出清晰的非线性模式——先经历一个边际效益递增阶段，随后转为递减。具体而言，在集群规模较小时，产出提升平缓；当集群达到中等规模时，集聚带来的知识溢出等正外部性显现，产出显著提升；然而，在集群规模极大时，边际效益增长放缓甚至可能出现轻微回落（如置信带右端所示）。通过对比发现，传统的残差化方法可能低估了这种非线性的强度。

基于置信带的统计推断为拒绝线性关系提供了正式依据。如图(5.6)所示，在整个解释变量取值范围内，无法找到一条直线能够完全被包含在均匀置信带之内。尤其是在中等集群规模区域，条件均值函数明显位于任何可能线性拟合的上方；而在最大规模区域，其增长则低于线性外推的预测。这一发现与Moretti等人（2021）关于效应递减的定性结论相一致，但本分析通过严谨的半参数推断框架，进一步量化并统计地证实了该非线性关系的具体形态与显著性。

## 6. 结论
本文通过系统的理论构建与实证检验证明，常用的Binscatter方法能够被系统性升级，转型为一个适用于严谨统计推断的非参数分析工具。

实证应用的结果证实了上述方法论改进的必要性与实质性影响。传统不严谨的操作可能产生具有误导性的可视化结论。例如，在评估边际税率对创新的影响时，错误方法得到的图形几乎无法识别出任何清晰模式；而应用本文框架，则能明确揭示并统计证实两者间显著的线性负相关关系，从而为原有的政策结论提供了更稳健的图形证据支持。再如，对高科技集群规模效应的再分析，经本文方法揭示出明确的边际效应递减非线性动态。这一细微但稳健的非线性特征对区域创新政策设计具有关键启示：针对处于不同发展阶段（小型初创期与大型成熟期）的集群，所应采取的扶持政策在目标与力度上应进行差异化设计。

## 7. 参考文献

本文的理论和实证部分主要参考下面的几篇文献。

> Cattaneo M D, Crump R K, Farrell M H, et al. On binscatter[J]. American Economic Review, 2024, 114(5): 1488-1514. [Link](https://www.aeaweb.org/articles?id=10.1257/aer.20221576)

> Akcigit U, Grigsby J, Nicholas T, et al. Taxation and innovation in the twentieth century[J]. The Quarterly Journal of Economics, 2022, 137(1): 329-385.[Link](https://academic.oup.com/qje/article-abstract/137/1/329/6292271)

> Moretti E. The effect of high-tech clusters on the productivity of top inventors[J]. American Economic Review, 2021, 111(10): 3328-3375. [Link](https://www.aeaweb.org/articles?id=10.1257/aer.20191277)
