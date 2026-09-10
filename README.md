# EPR-1 计算包（R4，2026-08-12）

课题 EPR-1"量子熵锥的完整刻画"（方案 A：n = 5 三锥判定化）的自足计算包。
无 psitip 依赖；全部输入数据与 Normaliz 静态二进制已捆绑；断点可续；进度可看。

## 0. 状态速览（先读这个）

**已裁决 / 已完成（不要重做）：**
- [x] s12 裁决（R3，双机复现）：psitip 清单缺 **Ingleton(39)**，判决
  PSITIP-INCOMPLETE；补救已入包（`CLR_H_fixed` 1,905 行为 s4 默认，
  `templates28` / **pure28** 23,355 行为量子侧现役）。
- [x] **当前头条 18/19**：pure28 下 19 条 HEC₅ 射线 18 条 extreme，非极只剩
  **#19**（精确秩 29）——与 SHC 论文宣称的 18 一致。**R11 已给出规范证书**：
  #19 的最小面是二维，两棱精确求出，**r₁₉ = 1·q₁ + 2·q₂**（q₁ = HEC#4，
  两棱各秩 30；data/cert19_exact.npz，G14 守护）。
- [x] **31 ↔ 31 向量级对映已机器认证（R12）**：BCHS 的"31 QLR 不等式"
  = DFZ 清单换弱单调（原文注钉出）；朴素族的 S₆ 类 = 34（= R2 之谜的
  34），其中 3 个弱单调类有效但冗余，余 31 与我们的 31 个 facet 类
  **集合级相等**——⑥ 之谜全解，无需任何新转录。
- [~] **A1 本体战役首日（R15，2026-09-10）**：s4c_qlr 在 QLR₅ 精确不可冗余
  H 表示（`qlr5_H_facets10860`，10,860 行 = R11 的 facet 实例数）上、以 59+q₂
  为种子做 S₆ 邻接分解，170 秒扩展 12 个最轻顶点图 → **1,924 个认证极轨道**
  （远未饱和：末图 1,067 邻居中 91% 为新轨道；坐标最大 162）。**量级结论：
  QLR₅ 极射线数比 CLR₅ 大数个数量级，全量枚举按当前引擎不可及**——这正是
  BCHS 2021"算力所限"的成因。可陈述成果：**轨道数下界 ≥ 1,924（已知 59）**、
  带证书目录 `qlr5_orbits_partial`（G17 守护）、S7 优先测试集
  `qlr5_new_small200`（新轨道中坐标最小的 200 条，max 坐标 3..13）。
  **解读护栏**：新轨道 1,864 条零条被已知可实现池见证——这只说明池的熵尺度
  （≤12）够不到新射线（坐标达 162），**不是对 S7 的证据**。
- [x] **59 = 40+17+2 台账全链条独立复现（R13）**：S₆ 口径重叠恰 HEC#1、
  图态新极轨道恰 2（G11/G15，后者以 (4.4) 全轨道最小 −1 鉴定）；
  59 资产入包（`qlr59_reps`，G16 守护）。**定理注记：31↔31 facet 相等
  ⇒ pure28 锥 = QLR₅**；q₂ 为 QLR₅ 真极射线且避开全部可实现池
  = **S7 猜想的具体试金石**（A2 首要对象）。
- [x] **#18 悬案裁决（R12）**：SHC Table 3 机器解析 + S₆ 轨道级双射：
  **论文点名的第 18 条 = 我们的 #19**（G15 守护）——18/19 头条与
  R11 精确证书 r₁₉ = HEC#4 + 2·q₂ 至此与文献逐环咬合；仓库序与论文序
  在 18/19 位互换是旧担忧的根源。59 = 40+17+2 台账公式已由原文全解。
- [x] **s4/s5 已达成（R7）：162 个 CLR₅ 极射线轨道全部找到**，轨道大小合计
  = **7,943 全射线**，与文献双口径精确吻合；每轨道代表逐一认证极性
  （紧秩 30）。引擎史：Normaliz primal 爆炸（R4）→ lrs 全量长尾不可行
  （R6 服务器 12.7 h / 9,700 万基 / 12% 射线，估计器低估 68 倍，勘误 C3）
  → **s4c 带对称邻接分解**十几分钟收官。数据已捆绑
  （`clr5_orbit_reps162` / `clr5_rays7943`，G12 关卡守护）。
- [x] **s6 已达成（R7）：162 代表全部在量子锥内，恰 40 个保持 QLR-extreme**
  ——台账头数精确复现，且经格拉姆精确秩认证（非极分布 25×56、26×25、
  27×18、28×23）。
- [x] **完备性互证已在逻辑层闭环（R9）**：DFZ §4 计数原句已钉（7,943 条 /
  162 轨道）+ 极射线集合唯一性 + 我们的 162 份精确证书 ⇒ 两集合必然相等；
  论文印出的三个射线例逐坐标命中 3/3。rays5 文件降级为可选加固
  （获取途径与邮件模板见 §7① 与迭代9报告）。122/162 顶点图闭合与
  7,943 总数吻合为计数提供独立佐证；mplrs 从零全量仅为发表奢侈品。
- [x] 27 类互不冗余证书、pure(27) 16/19 为历史口径，保留作回归关卡。

**待办按序（详见 §7）：** 服务器一条命令 `sh scripts/run_b_batch.sh`
（s7 → s8 → s10 → s11，约 25 分钟）→ 回传 s11 tar；离线：宋博士 11 条
单源抽查、（可选）rays5 加固；下轮分析：#19 ↔ SHC Table 3 对号 +
BCHS 31 对账。

## 1. 环境需求

- python3 ≥ 3.8，`pip install -r requirements.txt`（numpy、scipy）
- s7 需要 gcc（OpenMP 可选但推荐；无则单线程并有警告）
- lrs 仅在复算/取种子时需要（可选）：`apt-get install -y lrslib`
  （或 conda-forge lrslib）。交叉验证用 Normaliz：**包内捆绑静态二进制 `bin/normaliz`（v3.11.1，GPL-3）**，
  脚本按 `$NORMALIZ` → PATH → 捆绑件顺序解析。想固定用捆绑件（避免 PATH 上
  其它版本抢先，R4 服务器实际用的是 conda 的 3.11.0）：
  `NORMALIZ=$PWD/bin/normaliz bash scripts/s4_run_normaliz.sh ...`
- 仅限 Linux x86_64（并行层依赖 fork；捆绑二进制为 Linux 静态链接）

## 2. 快速开始（25 vCPU / 90 GB 参考耗时）

```
# ---- A. 验证已达成结果(可选, 共 ~3 分钟) ----
python3 selftest.py --full            # 16 关卡, ~35 s, 必须全绿再往下
python3 scripts/s1_build_qlr.py --variant pure28 --workers 25     # 秒级
python3 scripts/s2_judge.py QLR_H_pure28.npy                      # 期望 violated 0
python3 scripts/s3_rank19.py QLR_H_pure28.npy                     # 期望 18/19, 仅 #19 秩 29
python3 scripts/s12_dfz_reconcile.py                              # 期望 PSITIP-INCOMPLETE/ing39

# ---- B. 当前服务器工作: 一条命令 ----
sh scripts/run_b_batch.sh
#   = s7(两层 0 违反) -> s8(#19 证书; 无条件池优先, 失败自动升级并标注条件层)
#     -> s10(facet classes N + redundant M = 84) -> s11 打包
#   结束后把 s11 生成的 tar.gz 上传回对话

# ---- C. (可选加固) rays5 逐向量比对 ----
python3 scripts/s5b_diff_rays.py rays5     # 期望末行 MATCH; 文件获取见 §7①
```

预期数字：s2=0 违反；s3=18/19 且仅 #19 秩 29（sha16 `827079728c051bfe`）；
s7 两层均 0 违反；s10 输出 facet 类数供对账 BCHS 31；s5b 末行 MATCH。
**枚举已完结**：162 轨道 / 7,943 射线随包分发（G12 守护），无需重跑；
复算路径（可选）：s4_clr_prepare → s4b 短跑取种子 → s4c 批式续跑至闭合。

## 3. 十二阶段一览

| 阶段 | 命令要点 | 作用 / 期望 |
| --- | --- | --- |
| selftest | `selftest.py [--full] [--workers N]` | 14 关卡回归；服务器上第一件事 |
| s1 | `--variant pure28`（当前）/ pure、v3、pure2（历史） | 构建 H 表示；秒级 |
| s2 | `s2_judge.py <H.npy>` | 760 图态裁判；pure* 必须 0 违反 |
| s3 | `s3_rank19.py <H.npy> [--extra-rows X]` | 19 射线紧秩判定；pure28 → 18/19 |
| s4 | **`s4c_adjacency.py`**（`--init` 种子 → `--state` 批式续跑） | 极射线枚举主引擎（带对称邻接分解；R7 实战 162 收官）；s4b(lrs)/s4(normaliz) 为种子与交叉验证 |
| s5 | `s5_orbits.py clr5.out [--raw] [--expect 162]` | 极射线 → S₅ 轨道；对账 162 |
| s6 | `s6_sweep_clr_orbits.py reps.npy --H pure28` | CLR 极射线过量子锥；对账台账 40 |
| s7 | `--H pure28 [--f3-limit N] [--gf2-samples N]` | 扩可实现池（F₃ 全枚举 + 12 顶点采样） |
| s8 | `--H pure28 --pool <s7 输出>` | **#19** 的面限制分解证书（默认 targets1） |
| s9 | `s9_sixvar.py <人工转录 CSV>` | 6 变量清单代入（备选；铁律：只抄文献） |
| s10 | `--H pure28 [--limit N]` | facet 类 → 对账 BCHS 31 |
| s11 | `s11_collect.py <globs> --out tar.gz` | 结果打包 + 规范 sha 清单 |
| s12 | `s12_dfz_reconcile.py`（默认捆绑 CSV） | 复验 R3 判决（PSITIP-INCOMPLETE / ing39） |

## 4. 进度查看

- selftest：逐关卡 PASS/FAIL。
- s1 等并行构建：每 ≤5 秒一行 `k/n chunks + ETA`（前 2 秒 ETA 显示 `--`）。
- s4c（复算时）：每扩展一个代表打一行 `tight=… 邻居数/新轨道/总数/已闭`，
  逐代表落盘（`--state`），随时中断续跑；巨型图超 `--rep-timeout` 会明说跳过。
- lrs/Normaliz（种子或交叉验证时）：`wc -l *.lrs.out` / `tail -f clr5.log`；
  终止用 `pkill -x lrs` / `pkill -x normaliz`（**不要** `pkill -f`，R2 曾误杀）。
- s7：两个 C 工具向 stderr 报计数（F₃ 每 2²⁰ 配置一行带百分比；GF(2) 每 1M 样本一行）；
  字节核对：F₃ 全枚举总输出 = 14,348,907 × 31 ≈ 445 MB。
- s10：每 10 类一行。断点续跑分片启动时报告已完成/待跑数。
- 内存守护：`EPR1_MAX_RSS_GB` 只管本包 worker，**不管 Normaliz**。

## 5. 常见坑（真实发生过的在前）

1. **建错变体**：`--variant pure` 是历史 27 模板口径（16/19）。当前一律 `pure28`。
   已发生一次（R3 服务器日志）；s3 若输出 16/19 且 #17/#18 秩 29，先查这个。
2. **s4 无 normaliz**：脚本会明确报错并给两条出路（捆绑件 / conda）。捆绑件若报
   权限，`chmod +x bin/normaliz`。
3. **枚举引擎史（三次裁决，勿走回头路）**：Normaliz primal 中间爆炸
   （91.5M 超平面 @ gen 62/1905，R4）；lrs 全量长尾不可行（12.7 h /
   9,700 万基 / 12% 射线，估计器低估 68 倍——勘误 C3，R6）；
   **s4c 带对称邻接分解十几分钟收官（R7）**。结果数据已随包分发，
   常规使用不需要任何枚举。
4. **s5 MISMATCH**：先确认吃的是 `CLR_H_fixed` 产物（1,905 行；`clr5.in` 头部
   `inequalities 1905`）。旧 `CLR_H`（1,875）会差一截。
5. **rays5 原始数据**：Zeger 站 `code.ucsd.edu/zeger/linrank/` 禁自动抓取，
   需人工浏览器下载；拿到后 `s5_orbits.py rays5 --raw` 可与我们的枚举逐条 diff。
6. **ETA 首行观感**：并行首个分片完成后的外推值不可信，几秒后自然收敛。
7. **s7 无 gcc / 无 OpenMP**：分别是明确报错与降级警告，见 §1。
8. **s9 铁律**：6 变量不等式只允许从论文转录（模板见 `data/sixvar_template.csv`），
   逐变元平衡校验不过会拒收——那几乎总是抄写错误。
9. **估计器会骗人**：lrs 深度-2 树规模估计在极不平衡树上系统性低估
   （R6 实测低估 68 倍）。长跑决策看 checkpoint 的基数计与射线增速，
   不看估计值。
10. **巨型顶点图**：s4c 对 tight ≳ 1,000 的超对称射线的图枚举与原锥同难，
   会按 `--rep-timeout` 跳过并如实报告"closure incomplete"。
11. **小时级任务只在服务器上跑**：交互式容器/笔记本会话会在轮次间收割后台进程
   （R5 实测两次，皆截断在半行）。lrs 无状态，被杀直接重启即可；已出流的射线
   仍是真数据（s5 解析器自动丢弃截断的末行）。

## 6. 数据清单（data/）

| 文件 | 行数 | 含义 |
| --- | :-: | --- |
| `M_LR.npy` | 1,790 | psitip 提取的 LR 实例族（27 类；历史保留） |
| `CLR_H.npy` / `CLR_H_fixed.npy` | 1,875 / **1,905** | 经典锥 H 表示（修复版 = +ing39 轨道 30 条） |
| `REF28.npy` + `dfz_ref28.csv` | 28 | 文献参考族（24 DFZ + 4 Ingleton 形式；CSV 含来源标签与单源待查清单） |
| `templates27.npy` / `templates28.npy` | 27 / 28 | 代入模板（28 = 当前） |
| `ING39_orbit.npy` | 30 | Ingleton(39) 的 S₅ 实例 |
| `graphstate_vecs.npy` | 760 | 6 顶点 GF(2) 图态熵向量（裁判） |
| `hec5_rays_maskorder.npy` | 19 | HEC₅ 极射线轨道代表（掩码序） |
| `targets1.npy` / `targets3.npy` | 1 / 3 | s8 目标：#19（当前）/ #17,18,19（历史） |
| `sep_certs_27.npz` | 27 | R3 互不冗余精确整数证书（历史记录） |
| `clr5_orbit_reps162.npy` | **162** | CLR₅ 极射线轨道代表（R7 邻接分解产出，逐一秩 30 认证） |
| `clr5_rays7943.npy` | **7,943** | 全射线集（Σ轨道；与文献双口径吻合；G12 守护） |
| `manifest_shas.json` | — | 全部行族的规范 sha256（`sha_rows` 口径） |

构建产物（`QLR_H_*.npy`）不捆绑，由 s1 现场生成并对 sha。

## 7. 下一步：做什么、怎么做

**⓪ A1 战役续跑（服务器，小时级，可断点）**
```
python3 scripts/s4c_qlr.py --init data/qlr_seeds60.npy --state qlr_adj.npy   # 或从容器状态续
python3 scripts/s4c_qlr.py --state qlr_adj.npy --budget 14400 --rep-timeout 1800 --max-tight 400
```
目标不是闭合（不可及），而是：更紧的下界、更多小坐标测试用例、增长曲线
（每批记录"总轨道数 vs 已扩展数"）。产物：`qlr_adj.npy` 状态文件回传。
A1 按计划 §5.1 止损条款改述为 **A1′：认证部分目录 + 下界 + 结构统计**。

**⓪′ GitHub 仓库（建议：现在开私有库，成文时转公开）** 见 §9。


**① 162 完备性互证——已在逻辑层完成（R9），文件降级为可选**
命题：极射线集合是锥的内蕴不变量。我们的 162 个轨道代表两两不等价且各带
精确紧秩 30 证书（⊆ ExtRays(C)）；DFZ §4 原文明确 C 的极射线恰 7,943 条 /
162 轨道；故两个集合**必然相等**——逐向量 diff 在逻辑上冗余。抽验加固：
论文自己印出的三个具体射线（F³ 五子空间例、U₂,₄ 例、U₂,₅ 例）均逐坐标
命中我们的 7,943 集。诚实注记：计数前提是 DFZ 的 cddlib 计算（2009，
每次 31 维计算耗时 2–3 天）；我们 122/162 的邻接闭合与总数吻合独立佐证。
**可选加固**（若仍想拿文件）：曼城本机/eduroam 浏览器开
`http://code.ucsd.edu/zeger/linrank/rays5`；或宋博士从武大网络 `curl -O`；
或邮件向 zeger@ucsd.edu 索取（模板见迭代9报告）。拿到后
`python3 scripts/s5b_diff_rays.py rays5` 一键比对。

**② 服务器会话（一条命令，约 25 分钟）**
`sh scripts/run_b_batch.sh`。逐段判据：
- s7：两层各打印 `pure violations 0`；产物 `s7_out/realizable_gf2_12v.npy`
  （无条件层）与 `s7_out/realizable_f3_conditional.npy`（**条件层**，其有效性
  前提见坑清单，勿混淆）。
- s8：先用无条件池；`certificate saved` = #19 非极性证书落 `s8_certs.npz`
  （λ 权重即证书）。若 `no certificate`，脚本自动合并条件层重试并把结果
  **另存** `s8_certs_conditional.npz` + 标注前提——两种证书严谨等级不同，
  报告时须分开陈述。两者皆无 → 记开放，升级 s9（6 变量清单，人工转录）。
- s10：末行 `facet classes N  redundant M`（N + M = 84）。
- s11：自动打包（含三份日志）——**把这个 tar.gz 上传回对话**。

**③ 回传后（Claude 下轮分析）**
用 s10 的 `pure_facet_classes.npz` 做 BCHS 31 对账（⑥ 32/34 之谜）；
做 #19 ↔ SHC 1903.09148 Table 3 逐坐标对号（hec lex-by-size ↔ 掩码序的
映射公式已备案）。

**④ 外部时钟**
SHC 截稿 2026-08-25：届时弹药 = 162+7,943 双口径 + 完备性逻辑闭环 +
论文三例命中 + 40 台账精确认证 + ing39 定位 + 18/19 头条 + #19 证书。

**⑤ 宋博士：11 条单源抽查（离线，随时可做）**
对照 DFZ 论文式 (10)–(23) 核对 `data/dfz_ref28.csv` 中 source=psitip 的
11 条：dfz10, 11, 13, 14, 15, 16, 19, 20, 21, 22, 23（CSV 注释行印有每条的
原始 I-记法，逐项对系数即可）。结论不影响既有判决，但完成转录闭环。

**⑥ 可选奢侈品：从零完备性证明**
mplrs 并行全量（`apt`/源码装 mplrs + openmpi，25 进程，天级，可断点）；
或对 ~40 个巨型顶点图（tight ≳ 1,000）单独攻坚。仅发表需要时再做。
## 8. 版本

见 `CHANGELOG.md` 与 `VERSION`。R4 要点：捆绑 Normaliz、13 关卡自检、
targets1 默认、越界与 Bareiss 精确性断言、全脚本 `--help` 可用、本 README 重写。

## 9. GitHub 仓库工作流（R15 评估结论：开，先私有）

**为什么开**：轮次打包已暴露两类 git 会当场抓住的问题（C4 静默补丁失效 =
一次 `git diff`；C5 日期惯性 = 提交时间戳）；服务器同步从"上传 zip"变成
`git pull`；CI 每次推送自动跑 20 关卡；成文时挂 Zenodo DOI；给 SHC 仓库开
issue/PR 本来就需要 GitHub 身份。**为什么私有**：q₂ 与新射线目录在人工
新颖性留痕与成文前不宜暴露给活跃的 HR 生态；私有库免费且可加协作者。

**怎么做（你本机，约 10 分钟；我不接触任何令牌/凭据）**：
1. GitHub 新建私有仓库 `epr1-qlr5`（不勾选自动生成 README）；邀请宋博士为协作者。
2. 解包本 tar 后在 `kit/` 内：
   `git init && git add -A && git commit -m "R15: kit with 20 gates, QLR5 partial catalogue"`
   `git branch -M main && git remote add origin git@github.com:<你>/epr1-qlr5.git && git push -u origin main`
3. 已随包提供：`.gitignore`（忽略可重建的 pure28、池、日志、状态文件）、
   `.github/workflows/selftest.yml`（Actions 装 lrslib + 依赖后跑 `selftest.py --full`）、
   `LICENSE`（代码 MIT；数据 CC-BY-4.0；捆绑 Normaliz 为 GPL-3 二进制，来源见 bin/）、
   `CITATION.cff`（占位，成文时填）。
4. 服务器：`git clone` → `pip install -r requirements.txt` → 跑批 → 结果文件用
   `git lfs`（>50 MB 的池）或作为 Release 附件回传；小结果直接 commit 到 `results/` 分支。
5. 与我的往返：私有库我无法拉取（不持凭据）——继续上传 zip 即可；或临时
   公开一个只含结果的镜像库供 codeload 抓取。转公开时：加 Zenodo 集成拿 DOI。
