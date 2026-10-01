# Bambu Lab TPU90A インソールの G-code 構造を反映した超弾性 FEM 解析方針

## 1. 目的

Bambu Lab でスライス済みの `.gcode.3mf` から、実際の押出経路を反映したインソール形状を復元し、内部 Sparse Infill（Gyroid）を明示的に含めた状態で超弾性 FEM 解析を行う。

主目的は以下。

- Gyroid 内部の荷重伝達を確認する
- wall / infill / skin 接合部の局所応力・局所ひずみを確認する
- インソールの沈み込み量と反力の関係を得る
- 足裏と床の接触を含む圧縮解析を行う
- TPU90A の大変形特性を Mooney–Rivlin 等の超弾性モデルで反映する

解析環境の第一候補は **PrePoMax v2.6.2 + CalculiX** とする。

---

## 2. 入力 G-code の確認結果

対象ファイル：`Insole-L2.gcode.3mf`

確認できた主なスライス条件は以下。

| 項目 | 値 |
|---|---:|
| 総レイヤ数 | 40 |
| Layer height | 0.20 mm |
| G-code 内 nozzle diameter | 0.40 mm |
| Sparse infill | Gyroid |
| Sparse infill density | 15% |
| Wall loops | 2 |
| Outer wall line width | 約 0.42 mm |
| Inner wall line width | 約 0.45 mm |
| Sparse infill line width | 約 0.45 mm |
| Gyroid の存在レイヤ | Layer 4〜35 |
| 押出セグメント総数 | 約 90,233 |
| Sparse infill セグメント数 | 約 52,514 |

### 注意

今回のファイルでは、**0.20 mm はノズル径ではなく layer height** である。

ノズル径は G-code 内で **0.40 mm** と記録されている。

また、Gyroid は大量の `G2/G3` 円弧コマンドで記述されているため、変換プログラムで `G1` だけを処理すると内部構造が欠落する。

したがって変換器は最低限、以下を扱う必要がある。

- `G1`
- `G2`
- `G3`
- `; FEATURE:`
- `; LINE_WIDTH:`
- layer height
- extrusion amount `E`

---

## 3. 推奨解析ソフト

### 第一候補：PrePoMax + CalculiX

今回の問題は以下の要素を含む。

- 超弾性
- 大変形
- Gyroid の明示形状
- 接触
- 圧縮
- 将来的には自己接触の可能性

このため、Fusion よりも CalculiX 系の方が柔軟である。

PrePoMax は CalculiX の GUI として使用できる。

### Fusion

Fusion でも 2 定数 Mooney–Rivlin を使用できるが、今回のような細かい Gyroid を全体で扱うとメッシュ規模が大きくなりやすい。

また、材料係数の curve fitting は外部で実施する必要がある。

### SolidWorks Simulation

SolidWorks は超弾性材料の試験データから Mooney–Rivlin / Ogden 等の係数をフィットする用途に有用。

材料同定だけ SolidWorks を使い、その係数を CalculiX に移す方法も考えられる。

---

## 4. TPU90A 材料データについて

Bambu Lab TPU90A について詳細な stress–strain データが公開されていないため、初期検討では Formlabs TPU 90A の代表的な stress–strain 曲線を surrogate material として使う。

ただし、以下の違いがある。

- Formlabs TPU90A：SLS 系
- Bambu Lab TPU90A：FDM フィラメント

したがって絶対値の一致を保証するものではなく、初期の傾向解析・モデル成立確認用と考える。

---

## 5. Mooney–Rivlin の暫定係数

Formlabs TPU90A の代表的 stress–strain curve を参考に、非圧縮の 2 定数 Mooney–Rivlin に合わせた暫定値として、次を利用可能。

```text
C10 ≈ -0.44 MPa
C01 ≈ 4.99 MPa
```

画像全体を参考にした場合でも、おおむね以下の範囲となる。

```text
C10 ≈ -0.42 MPa
C01 ≈ 5.3 MPa
```

### CalculiX 入力例

```text
*MATERIAL, NAME=TPU90A
*HYPERELASTIC, MOONEY-RIVLIN
-0.44, 4.99, 0.
```

### 注意点

`C10 < 0` となること自体は、この stress–strain 曲線を 2 定数 Mooney–Rivlin だけで表そうとした結果である。

これは材料が異常という意味ではなく、**2 定数 Mooney–Rivlin の自由度が不足している可能性**を示す。

より精度を上げる場合は以下も検討する。

- Yeoh
- Ogden
- Arruda–Boyce

特に 100% を超える大ひずみ領域まで合わせる場合は、Ogden や Yeoh の方が適する可能性が高い。

---

## 6. 最終的な材料同定で望ましい方法

今回 G-code から Gyroid を明示的に再現するため、材料モデルに入れるのは **Gyroid 構造全体の見かけ特性ではなく、押出された TPU 材料自体の特性**であるべき。

したがって理想的には、Bambu TPU90A で以下を造形する。

- 100% infill
- Gyroid ではなく solid に近い引張試験片
- 実際のインソールと同じ nozzle / layer height / temperature / speed 条件

その試験片から一軸引張 stress–strain curve を取得し、超弾性係数を fit する。

### 重要

Gyroid 15% 試験片の見かけ stress–strain を材料定数として使い、その上で解析モデルにも Gyroid を入れると、Gyroid による柔らかさを二重に計上してしまう。

---

## 7. 推奨する G-code → FEM 変換経路

当初は STL 化も可能だが、PrePoMax / CalculiX を主解析環境とするなら、最終的には以下を推奨する。

```text
Bambu .gcode.3mf
        ↓
G-code parser
        ↓
G1 / G2 / G3 extrusion path
        ↓
voxel occupancy
        ↓
C3D8 / C3D8R hexahedral mesh
        ↓
CalculiX .inp
        ↓
PrePoMax
```

### この方法の利点

STL 経由では、

```text
G-code
→ 表面三角形化
→ STL
→ 四面体メッシュ
→ FEM
```

となり、細いビード部分でメッシュ失敗や幾何誤差が発生しやすい。

一方 voxel から直接 C3D8 系要素にすれば、押出形状を直接 solid element にできる。

---

## 8. 推奨 voxel サイズ

今回の代表寸法は、

- layer height = 0.20 mm
- infill line width ≈ 0.45 mm

である。

### 初期解析

```text
voxel = 0.20 × 0.20 × 0.20 mm
```

この場合、0.45 mm の押出幅に対して XY 方向約 2〜3 要素となる。

初期の傾向解析には使用可能。

### 精密解析

```text
voxel = 0.10 × 0.10 × 0.10 mm
```

にして mesh convergence を確認する。

### 注意

インソール全体を 0.10 mm voxel で解析すると要素数が非常に大きくなる。

例えば 270 × 100 × 8 mm の bounding box なら、単純換算で、

```text
2700 × 1000 × 80
= 216,000,000 voxel
```

となる。

空隙を除いても巨大になるため、まず ROI 解析が適切。

---

## 9. 推奨 ROI

最初は **heel region（踵部）** を切り出して解析する。

例：

```text
約 60 × 60 × 8 mm
```

その後、

1. heel ROI
2. forefoot ROI
3. 必要なら全インソール

の順に拡張する。

---

## 10. FEATURE ごとの element set

G-code の `; FEATURE:` 情報を保持し、CalculiX `.inp` に以下のような element set を作ると便利。

```text
ELSET_OUTER_WALL
ELSET_INNER_WALL
ELSET_SPARSE_INFILL
ELSET_SOLID_INFILL
ELSET_TOP_SURFACE
ELSET_BOTTOM_SURFACE
```

これにより PrePoMax 上で、

- Gyroid だけ表示
- Wall だけ表示
- Sparse infill の最大ひずみだけ確認
- Wall / infill 界面の応力集中確認

などができる。

---

## 11. FDM ビードの取り扱い

第一段階では、隣接ビード同士は完全に融着したものとして扱う。

```text
TPU bead ─┐
          ├─ continuous / bonded
TPU bead ─┘
```

つまり 1 つの連続体として voxel 化する。

実際の FDM では以下の影響がある。

- layer 間接着
- bead 間接着
- microscopic void
- 押出方向による異方性
- 残留応力
- 粘弾性
- ヒステリシス

ただし、最初のモデルではこれらを無視し、

```text
explicit G-code geometry
+ isotropic hyperelastic TPU
+ large deformation
```

で解析を成立させることを優先する。

---

## 12. インソールの荷重条件

単純な一点荷重よりも、以下の接触解析が適切。

```text
        足裏 / loading surface
                 ↓
          prescribed displacement
                 ↓
       ┌────────────────┐
       │    TPU skin    │
       │  Gyroid infill │
       │    TPU skin    │
       └────────────────┘
────────────────────────────
             floor
```

### 推奨：変位制御

初期段階では force control より displacement control を推奨。

例：

```text
足裏剛体：Uz = -2 mm
          Uz = -5 mm
          Uz = -10 mm
```

床側は固定。

```text
Ux = 0
Uy = 0
Uz = 0
```

### 理由

超弾性 + 接触 + Gyroid の局所座屈では、荷重制御より変位制御の方が解析が安定しやすい。

---

## 13. 足裏側のモデル

段階的に複雑化する。

### Stage 1

平板 rigid surface

### Stage 2

踵を模した球面・楕円面 rigid surface

例：半径 40〜60 mm 程度から試す。

### Stage 3

実際の足裏 STL

最終的には実足裏 geometry を使い、接触圧力分布を評価する。

---

## 14. 出力として見るべき量

TPU Gyroid では、最大 von Mises stress だけでは不十分。

推奨する主要出力は以下。

- displacement `U`
- reaction force `RF`
- contact pressure
- principal strain
- logarithmic strain / equivalent strain
- stretch ratio
- wall / infill junction stress
- Gyroid local strain
- bottom reaction distribution

特に重要なのは **局所ひずみ**。

Gyroid は、

- bending
- membrane stretching
- local buckling
- self-contact

によって変形するため、破壊や過大変形評価では応力だけでなく strain / stretch を追う。

---

## 15. 荷重–沈み込み曲線

足裏側を強制変位させ、その反力を合計すると、

```text
沈み込み量 [mm]
vs
反力 [N]
```

の曲線が得られる。

これはインソールの剛性評価として非常に有用。

例：

```text
-2 mm  → RFz = ... N
-5 mm  → RFz = ... N
-10 mm → RFz = ... N
```

実験の圧縮試験と比較すれば、材料係数の校正にも使える。

---

## 16. 推奨する初期 CalculiX 解析設定

| 項目 | 初期設定 |
|---|---|
| Solver | CalculiX |
| GUI | PrePoMax v2.6.2 |
| 要素 | C3D8 / C3D8R voxel mesh |
| voxel size | 0.20 mm |
| 領域 | heel ROI |
| 材料モデル | Mooney–Rivlin |
| C10 | -0.44 MPa |
| C01 | 4.99 MPa |
| D1 | 0 暫定 |
| 幾何非線形 | ON |
| 床 | rigid |
| 足裏 | rigid |
| TPU–床 | contact |
| TPU–足裏 | contact |
| 負荷 | 足裏側 prescribed displacement |
| 初期押込量 | -2 mm |
| 次のケース | -5 mm |
| 主出力 | U, S, strain, RF, contact pressure |

---

## 17. Mooney–Rivlin から先の材料モデル

2 定数 Mooney–Rivlin で Formlabs TPU90A の曲線を完全に再現することは難しい。

そのため、解析が成立した後に以下を比較する価値がある。

### Mooney–Rivlin

- パラメータが少ない
- 初期検討向き
- 中〜大変形域で fitting 能力が不足する場合あり

### Yeoh

- ゴム・TPU の大変形に比較的使いやすい
- 一軸データだけでも fit しやすい

### Ogden

- 大変形曲線への適合性が高い
- パラメータ数が増える
- TPU のような強い非線形性には有望

---

## 18. 解析精度を上げるための将来課題

### 材料側

- Bambu TPU90A の実引張試験
- 100% infill 試験片を使用
- 一軸引張曲線を取得
- 可能なら planar / biaxial データ取得
- Mooney–Rivlin / Yeoh / Ogden を比較

### FDM 側

- 層間界面強度
- 押出方向異方性
- bead cross-section の実測
- 実 line width
- 実 layer height
- 造形温度・速度の影響

### 接触側

- TPU と床の摩擦係数
- TPU と足裏の摩擦係数
- self-contact

### 動的特性

実際の歩行を評価する場合は、将来的に粘弾性や動的解析も必要。

TPU はヒステリシスを持つため、静的 hyperelastic だけでは歩行時のエネルギー損失までは再現できない。

---

## 19. 推奨作業フロー

### Step 1

`.gcode.3mf` を解析し、G1/G2/G3 の extrusion path を抽出。

### Step 2

各 extrusion segment を line width と layer height に基づいて voxel 化。

### Step 3

heel ROI を 0.20 mm voxel で C3D8 / C3D8R 化。

### Step 4

FEATURE ごとに ELSET を生成。

### Step 5

CalculiX `.inp` を生成。

### Step 6

PrePoMax v2.6.2 で読み込み。

### Step 7

TPU90A Mooney–Rivlin を設定。

### Step 8

床と足裏側 rigid surface を設定。

### Step 9

接触を設定。

### Step 10

足裏を -2 mm 強制変位して解析。

### Step 11

収束確認後 -5 mm、-10 mm へ拡張。

### Step 12

RF–displacement、contact pressure、Gyroid local strain を評価。

### Step 13

0.10 mm voxel と比較して mesh convergence を確認。

### Step 14

必要に応じて Yeoh / Ogden に移行。

---

## 20. 次に実装すべき Python プログラム

次の実装目標は、STL ではなく直接 CalculiX `.inp` を生成するスクリプト。

想定インターフェース：

```bash
python bambu_gcode3mf_to_calculix.py \
    Insole-L2.gcode.3mf \
    heel_0p2mm.inp \
    --voxel 0.20 \
    --crop xmin xmax ymin ymax zmin zmax
```

生成内容：

```text
*NODE
*ELEMENT, TYPE=C3D8R
*ELSET, ELSET=OUTER_WALL
*ELSET, ELSET=INNER_WALL
*ELSET, ELSET=SPARSE_INFILL
*ELSET, ELSET=SOLID_INFILL
*MATERIAL, NAME=TPU90A
*HYPERELASTIC, MOONEY-RIVLIN
...
```

必要なら境界 node set も自動生成する。

例：

```text
*NSET, NSET=BOTTOM
*NSET, NSET=TOP
```

---

## 21. 現時点の結論

今回の目的に対して最も実用的な構成は、

```text
Bambu Lab .gcode.3mf
        ↓
G-code geometry reconstruction
        ↓
voxel hexahedral FE mesh
        ↓
CalculiX .inp
        ↓
PrePoMax v2.6.2
        ↓
TPU90A hyperelastic material
        ↓
foot–insole–floor contact analysis
```

である。

最初は Formlabs TPU90A を参考にした暫定 Mooney–Rivlin 係数を用い、Gyroid 構造の応力・ひずみ分布とインソール全体の荷重–沈み込み特性を確認する。

その後、Bambu TPU90A の実造形引張試験データを取得して材料係数を再同定すれば、より定量的な解析へ移行できる。

