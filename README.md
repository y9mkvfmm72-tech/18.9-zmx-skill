# 18.9-zmx-skill

Zemax OpticStudio（老版 18.9 世代，`.zmx` 文件）镜头设计技能：
**自写光线追迹验算 + 纯 Python 优化处方 + .zmx 文件生成**。
生成的 `.zmx` 写 `VERS 181119`，同时兼容 18.9 与 2026R1。

## 这个仓库是什么

一份"踩坑沉淀"：本机（Windows + Zemax 18.9）多年实操中定位过的
**静默错误**（不报错、只给错答案）与可复跑的验证关卡。
所有结论都有本机实测数字，不是文档摘抄。

主文档：[`SKILL.md`](SKILL.md)

## 目录

```
SKILL.md          技能正文（版本选择、.zmx 格式、16 条静默 bug、优化流程）
scripts/          可复跑验证关卡（纯 Python 标准库，无第三方依赖）
  check_design_gates.py    六道关卡 A–F（追迹器+近轴+口径）
  matrix_firstorder.py     2×2 矩阵一阶量（第三条独立路径）
  footprints.py            通光口径/暗角检查
  tracer.py                真实光线追迹
  paraxial.py              近轴 y-u 引擎
  naive_ref.py             双胶合朴素独立实现（关卡 A 的对照）
  glass_index.py           SCHOTT.AGF 解析（CD 行交错列，带自检）
  check_zmx_structure.py   .zmx 结构比对（对 Zemax 自己写的文件）
  check_tracer.py          旧关卡：双胶合 + 近似法向对照
  README.md                关卡说明与已知正确基准
```

## 快速开始

```bash
cd scripts
python check_design_gates.py        # 报任何设计数字之前，先让这套全绿
```

## 已知正确基准（三条独立路径互证）

| 量 | 值 |
|---|---|
| Zemax 样板双胶合 EFL（0.633 µm） | **100.0439 mm** |
| 同上，近轴 l′ | **97.9387 mm** |
| 平面斯涅尔 30°、n 1→1.5 | **19.4712°** |

> 早期版本把双胶合 EFL 写成 99.04：那来自一个传递符号写反的近轴引擎。
> 教训本身也写进了 SKILL.md §3.16——**参照数本身也可能是错的，要三角互验。**
