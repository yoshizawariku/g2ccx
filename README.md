# g2ccx — Bambu G-codeからTPUの圧縮解析

Bambu Studioのスライス済み `.gcode.3mf` / `.gcode` を、押出経路を反映したvoxel六面体メッシュへ変換し、PrePoMax同梱のCalculiXで解きます。STLの四面体再メッシュは不要です。Python処理はuv環境で動作し、PrePoMaxのGUIを起動せず解析・集計できます。

## セットアップ

リポジトリ直下で実行します。

```powershell
uv sync --locked
uv run python PythonCode/bambu_fea.py inspect sample/Insole-L2.gcode.3mf --preview output/202609301817_gcode_xy_preview/sample_xy.png
```

Python 3.12の `.venv` と依存ライブラリをuvで管理します。CalculiX（PrePoMax同梱の `ccx_dynamic.exe` など）は同梱していないため、各自の環境で場所を指定します。次の優先順で探します。

1. JSON設定の `"solver": "C:/path/to/ccx_dynamic.exe"`
2. 環境変数 `G2CCX_SOLVER`
3. リポジトリ直下の `.g2ccx.local.json`（Git管理外）: `{"solver": "C:/path/to/ccx_dynamic.exe"}`
4. `PATH` 上の `ccx_dynamic` / `ccx` / `ccx_static`

どれも見つからない場合は、設定方法を示すエラーで停止します（解析以外のコマンドには不要です）。

通常のプロジェクト `.3mf` にはG-codeが含まれません。Bambu Studioからスライス後のプレートをG-codeとして書き出してください。複数プレートのG-codeが含まれる場合は、1プレートを別ファイルとして書き出すか、使用する `.gcode` を取り出します。

## 初回の踵ROI選択

```powershell
uv run python PythonCode/bambu_fea.py select-roi --config configs/heel.json
```

XY投影図で踵中心をクリックし、赤枠を確認して **Save ROI** を押します。初期サイズは60×60 mm、厚さ方向は全体です。閉じるだけなら保存しません。パン・ズームの操作中はROIを変更しません。

選択位置をJSONの `crop_mm` に保存し、同じ場所へ `.roi.png` も保存します。以後の解析では再選択不要です。JSONの `roi_size_mm` で矩形サイズを変更できます。切断面は造形座標のX/Y軸に平行です。

```powershell
uv run python PythonCode/bambu_fea.py run --config configs/heel.json
```

`run` は生成→解析→集計を行います。大規模モデルは後述の自由度上限で解析前に停止します。出力ディレクトリの上書きはしません。

## 小規模の実行例

先に実サンプルの6×6 mm領域を0.5 mm圧縮するケースで、環境と一連の処理を確認できます。

```powershell
uv run python PythonCode/bambu_fea.py run --config configs/sample_smoke.json
```

各段階を個別実行する場合:

```powershell
uv run python PythonCode/bambu_fea.py build --config configs/heel.json
uv run python PythonCode/bambu_fea.py solve --config configs/heel.json
uv run python PythonCode/bambu_fea.py postprocess --config configs/heel.json
```

`solve` は集計まで行い、`postprocess` は保存済みDATを再集計します。生成後に解析条件を変えた場合は、新しいディレクトリへ `build` し直してください。INPのハッシュも検査し、別の入力に結果を誤って結び付けません。

実行時設定の `solver`、`threads`、`timeout_seconds`、`max_solver_dofs` だけは再生成せず変更できます。実際の値を `run.json` に残します。

| 設定ファイル | 用途 |
|---|---|
| `configs/heel.json` | 使用者が踵を選ぶ通常の入口。0.20 mm voxel、2 mm圧縮 |
| `configs/sample_smoke.json` | 実サンプル6 mm角、0.5 mm圧縮の動作確認 |
| `configs/sample_small.json` | 同じ6 mm角、2 mm圧縮 |
| `configs/sample_small_fine_steps.json` | 2 mm圧縮を小さい増分で実行 |
| `configs/sample_small_c3d8.json` | 同じ2 mm圧縮をC3D8完全積分要素で比較 |
| `configs/sample_small_0p1.json` | 同じ6 mm角を0.10 mm voxelで生成・比較 |
| `configs/sample_heel_validation.json` | プレビュー右上を踵と仮定した60 mm角の生成検証。位置は使用者の選択で確定すること |

JSON内の相対パスは**JSONファイルの場所**を基準に解決します。ROIの保存では `crop_mm` だけを書き換え、他の記述は変更しません。

**出力フォルダ名**: 設定の `output` に `{timestamp}` を含めると、`build` / `run` / `select-roi` 実行時に現在時刻（`yyyymmddhhmm`）へ展開され、毎回新しいフォルダになります（例: `output/202610011530_sample_roi6mm_smoke_0p5mm`）。`solve` / `postprocess` / `view` など既存結果を使うコマンドは、同じパターンに一致する**最新のフォルダ**を自動で選びます。`output/` 直下はこの形式（`yyyymmddhhmm_何の解析か`）で整理してください。

`heel.json` は実サンプルの小領域で2 mm圧縮が完了した **C3D8、初期増分0.005、最大増分0.02** を明示設定しています。C3D8Rでは同じ大圧縮ケースに数値不安定が生じたためです。以下の表は設定省略時のCLI既定値です。要素種類を変えると剛性も変わるため、C3D8の収束だけを精度保証とはしません。

## 解析設定とモデルの意味

| 項目 | 既定値・扱い |
|---|---|
| 単位 | N、mm、MPa |
| 形状 | 楕円体のビードを経路に沿って掃引した近似、セル中心で占有判定 |
| voxel | `voxel_mm: 0.2`。格子は造形座標原点に固定 |
| 材料 | 等方性Mooney–Rivlin、C10=-0.44、C01=4.99、D1=0 |
| 要素 | `element_type: "C3D8R"`、節点共有による完全融着。`"C3D8"` も指定可 |
| 解析 | 静的・幾何非線形、変位制御 |
| 床・上板 | 剛体、上板Uzを負方向に指定 |
| 接触 | 面対面penalty、摩擦なし。`friction > 0` のとき摩擦カードを追加 |
| 接触剛性 | `contact_stiffness_mpa_per_mm: 1000` |
| 線形ソルバー | `linear_solver: "PARDISO"`。同梱2.22で検証 |
| 増分 | 初期0.01、最大0.05、最小0.000001（全ステップ時間1） |
| 出力頻度 | 5増分ごと＋最終増分。`output_frequency: 1` なら毎増分 |
| 並列 | 最大4スレッドを初期値にする。`threads` で変更 |
| 時間上限 | `timeout_seconds: 7200` |

ROIの切断面は自由境界です。これは**切り出した試験片の圧縮モデル**であり、周囲のインソールの支持を伝えるサブモデリングではありません。底面の中央に近い節点のUx/Uyと、同じY座標で離れた節点のUyを拘束して、面内の剛体移動・回転だけを抑えます。拘束位置と残留反力を記録します。摩擦を追加する場合も、この拘束の影響を確認してください。

材料係数は仕様書の暫定値で、Bambu TPU90Aの実測同定値ではありません。CalculiXでD1=0を指定すると、初期ポアソン比0.475相当の圧縮性が採用されます。完全非圧縮の指定にはなりません。材料・接触・解像度の感度を確認してから定量評価に使用してください。

自己接触、粘弾性、異方性、層間剥離、足裏形状は初期版に含みません。大圧縮で内部構造が接触する段階はこのモデルの適用限界です。入力の2 mmや5 mmという値だけでは収束・物理的妥当性を保証しません。

C3D8Rは薄いビードが厚さ方向1要素となる場合、曲げ・局所変形の精度に限界があります。C3D8は完全積分による過剛性にも注意が必要です。どちらかが収束することだけを根拠に精度を判断せず、要素種類・解像度の比較を行ってください。

## 要素数と接続性の診断

- `max_grid_cells`: 5千万セル。メモリ確保前に検査します。
- `max_elements`: 200万要素。節点・要素配列生成前に検査します。
- `max_solver_dofs`: 50万自由度。ソルバー開始前に検査します。十分な計算資源がある場合に明示的に増やしてください。
- `assembly_estimate_mib` は組立用の粗い見積りです。ソルバーの疎行列分解に必要なメモリは別で、上限保証ではありません。
- 全体の面接続成分が複数ある場合は解析を停止します。孤立部分を勝手に削除・架橋しません。
- `local_nonmanifold_vertices` は、節点近傍のセルが面接続せず辺・頂点で接する箇所数です。全体が単一成分でも発生します。自動修復はせず、細分化と局所結果の確認に使います。

実サンプルの60 mm角・0.20 mmモデルは約90.5万要素、124万節点、372万自由度です。通常の自由度上限を超えるため、そのまま `run` するとモデル生成後に解析を停止します。先にROIを小さくする方法を推奨します。粗いvoxelは細いビードを欠落させる可能性があるため、単に計算を軽くするための大幅な粗大化は避けてください。

## 出力とPrePoMaxでの確認

| ファイル | 内容 |
|---|---|
| `model.inp` | TPU、FEATURE別ELSET、剛体平板、接触、境界条件、解析ステップ |
| `model.frd` | 変位、応力、ひずみ、反力、接触変位・圧力の場 |
| `model.dat` / `.sta` / `.cvg` | 履歴・増分進捗・収束ログ |
| `solver.log` / `run.json` | 実行記録、終了コード、ソルバー版、入力ハッシュ |
| `force_displacement.csv` / `.png` | 保存増分の沈み込み–反力曲線、上下反力の釣合い |
| `feature_fields.csv` | FEATUREごとの積分点応力・ひずみの極値 |
| `summary.json` | 目標変位への到達と完了判定、釣合い検査 |
| `geometry.json` / `input.json` | 接続性、要素数、境界、押出・スライス条件 |
| `config.resolved.json` / `build.json` | 実際に使用した設定と入力ハッシュ |
| `features.npz` | 集計用のFEATURE所属。重複所属あり |
| `geometry.stl` | `export_stl: true` 時の確認用表面 |

PrePoMaxでは `model.inp` をImportしてメッシュ・要素集合を確認し、解析後の `model.frd` をOpenして結果を確認します。アプリ内で解き直す場合、インポート時に未対応カードがないか確認してください。自動パイプラインが実行する元のINPを解析条件の基準とします。専用の `.pmx` プロジェクトは生成しません。

`E` はCalculiXのLagrangeひずみです。CSVでは主ひずみから `sqrt(1+2E)` で主伸長比、`log(1+2E)/2` で主対数ひずみを導出します。応力はvon Mises、ひずみの偏差相当量は `sqrt(2/3 * dev(E):dev(E))` です。これらは破壊判定や塑性ひずみを意味しません。FEATUREの重複領域は複数集合に現れるため、集合間の値を単純加算しないでください。

終了コード0に加え、正常終了ログと目標変位・最終時刻への到達を検査します。CLIは上下反力差が1%を超えても非ゼロで終了します。未収束・タイムアウトは途中結果を残し、`incomplete` として集計します。FRDが存在することだけでは完了を意味しません。

## 結果ビューア（ヒートマップ）

```powershell
uv run python PythonCode/bambu_fea.py view --config configs/sample_smoke.json
uv run python PythonCode/bambu_fea.py view --output output/202609301936_sample_roi6mm_smoke_0p5mm --save output/view.png
```

`model.frd` を読み、XY・XZ・YZの3断面をヒートマップ表示します。X/Y/Zスライダで断面位置、Incrementスライダ（または←→キー）で増分、ラジオボタンでvon Mises応力・主応力・相当ひずみ・変位を切り替えます。下段の荷重–沈み込み曲線に現在の増分を表示します。`--save` はウィンドウを開かず最終増分をPNGに保存します。

- 値は要素8節点の節点値（CalculiXが積分点から外挿）の平均です。FEATURE別の積分点極値は `feature_fields.csv` を参照してください。
- 断面は変形前の格子座標に描画します。カラースケールは既定で最終増分の範囲に固定し、`Scale per increment` で各増分に合わせます。
- 要素は8節点六面体のみ対応します。

### 3D表示

```powershell
uv run python PythonCode/bambu_fea.py view --3d --output output/202609301936_sample_roi6mm_smoke_0p5mm
```

PyVista（VTK）で六面体メッシュを変形形状に重ねて表示します。マウスで回転・ズーム、下部のスライダで増分・変形倍率・Y方向のカット（内部の応力確認）を操作します。キー `1`〜`6` で表示量、`←→` で増分、`d` で変形/未変形、`s` でカラースケールの固定/増分ごとを切り替えます。`--save out.png` でウィンドウなしの画像出力ができます。

インソール全体の中に解析ROIを置いて見るには `--context` を付けます。

```powershell
uv run python PythonCode/bambu_fea.py view --3d --context --output output/<yyyymmddhhmm>_heel_roi60mm_voxel0p2_2mm
```

G-code全体を0.2 mmでvoxel化してXYのみ0.4 mmに間引いた灰色の表面を描き、ROI部分だけを応力コンターに置き換えます。初回は数分かかり、結果は出力ディレクトリの `context_*.npz` に保存されます。半透明の青い板は、ROIに与えた剛体上板（圧縮量ぶん下降）です。**荷重がかかるのはROIの上面全体のみで、インソール全体の応力解析ではありません**（全体を0.2 mmで解くと約400万要素で、上限50万自由度を大きく超えます）。

## 足裏全体の荷重解析（均質化モデル）

インソール全体をビード解像度で解くことはできないため、G-codeの内部構造（ジャイロイド・グリッド等）を**タイル圧縮で見かけの材料に均質化**し、外形全体を1.2 mm角の柱×6層の連続体として解きます。Moticon OpenGoの16ch圧力を上面荷重として与えます。設定は `configs/insole_global.json`。

```powershell
uv run python PythonCode/bambu_fea.py insole prepare   --config configs/insole_global.json   # 柱の高さ・充填率をクラス分けし、代表タイルを選定
uv run python PythonCode/bambu_fea.py insole tiles     --config configs/insole_global.json   # クラスごとの10 mm角タイルを0.2 mm voxelで圧縮（並列・数時間）
uv run python PythonCode/bambu_fea.py insole calibrate --config configs/insole_global.json   # 応力–ひずみ曲線に圧縮性Ogdenを当てはめ materials.json
uv run python PythonCode/bambu_fea.py insole build     --config configs/insole_global.json   # 荷重ケースごとにINPを生成
uv run python PythonCode/bambu_fea.py insole solve     --config configs/insole_global.json   # 全ケースを解く（1ケース数分）
uv run python PythonCode/bambu_fea.py view --3d --output output/<yyyymmddhhmm>_insole_global_gait_loads/cases/standing
```

- **インフィル形状の反映**: 各柱の高さ・全体/底面/上面の充填率でクラス分けし、クラスごとの代表タイルを**実際のG-codeのままビード解像度で圧縮**して材料を同定します。グリッドなど別のG-codeに差し替えて `prepare` からやり直すと、応答の違いが全体解析に反映されます。
- **荷重**: `260410130732_gait.txt` の左足16chの圧力を、`plan/Insoleセンサ座標.txt` の位置（踵端からの距離y、内側が正のx）でインソール座標に置き、ガウス核で補間します。合計が指定の力になるよう規格化します。足の向きは外形の主軸から自動判定します（踵はX,Y大の側、内側はくびれの深い側）。
- **ケース**: `standing`（片脚で0.5×体重、立脚中期の平均パターン）、`heel_strike`（踵荷重最大）、`push_off`（前足部荷重最大）、`walk_01`〜`walk_08`（代表的な1立脚を8等分）。歩行時の力は、記録の立脚ピーク平均が1.1×体重（`peak_factor`）になるよう換算します。体重は `body_mass_kg`（既定53 kg）。
- **境界条件**: 底面は床に接して上下方向拘束（摩擦なし、面内の剛体移動のみ拘束）、上面の圧力は柱面積を4節点に配分した鉛直の節点力で与えます（ひずみが小さいため非追従）。反力と荷重の差は `cases/*/summary.json` の `balance_relative`（全11ケースで0.1%未満）。CalculiXの残差判定は節点力の平均値で規格化されますが、この規模では平均が微小で既定値だと収束判定が成立しないため、`*CONTROLS` で許容を緩めています（許容を変えても結果は一致することを確認済み）。

**材料同定の範囲**: タイル圧縮は発散で途中終了したクラスがあり（クラス0は圧縮が約3%・14 kPaまで）、足裏荷重の範囲（〜0.25 MPa）の応力–ひずみ曲線に圧縮性Ogden（初期ポアソン比0.15固定）を当てはめています。同定範囲を超える圧力が掛かる柱の数は `cases/*/build.json` の `classes_beyond_calibration` に記録され、その部分の結果は信頼できません（クラス0は踵接地・歩行の複数ケースで該当）。

**限界**: 見かけの等方性超弾性なので、ビード単位の局所応力は出ません。タイル圧縮は側面が自由な小片で、周囲の拘束を含みません。柱内の高さ方向の構造（表層と内部の剛性差）は平均化されます。荷重の力への換算とセンサ配置は仮定を含みます。

## STLのみの既存CLI

```powershell
uv run python PythonCode/bambu_gcode3mf_to_fea_stl.py sample/Insole-L2.gcode.3mf output/roi.stl --voxel 0.2 --crop 195 201 190 196 0 8
```

既存CLIも共通voxel処理を使用します。旧実装の点の丸めとstamp方式からセル中心判定へ変更したため、以前のSTLと形状が完全に同一とは限りません。ROI境界は最も近い格子面へスナップし、実際の境界はレポートに記録します。

## テスト

```powershell
uv run pytest -q
uv run python PythonCode/verify_pipeline.py --output output/block_check_0p2 --pitch 0.2
uv run python PythonCode/verify_pipeline.py --output output/block_check_0p1 --pitch 0.1
```

検証ブロックは2 mm角、0.2 mm圧縮です。新しい出力先を使用してください。これは接触・材料入力・反力集計の確認であり、インソールの局所応力がメッシュ収束したことを示す試験ではありません。実行済みの検証状況は `VALIDATION.md` に記録します。

参考: [解析方針](Insole_TPU90A_PrePoMax_FEA.md)、[CalculiX 2.22マニュアル](https://www.dhondt.de/ccx_2.22.pdf)。
