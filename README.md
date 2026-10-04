# Concatemer Decode Service

从滚环扩增（RCA）产生的串联条码读段中，在插入、缺失和替换噪声下恢复共同切点。

联合选择：参考序列的一个**方向**（正向或反向互补）及其**循环移位**、恰好 `copies` 个覆盖整条读段的
**连续非空分段**，以及每段的一次**全局比对**（Needleman–Wunsch，M/D/I 单碱基代价均为 1），按

1. 总编辑数（`total_edits`）
2. 单段最大编辑数（`max_segment_edits`）

字典序最小化；最优解释按 `(strand, shift, boundaries, CIGAR)` 稳定排序（省略 `strand_mode` 时为 `(shift, boundaries, CIGAR)`）。

- 唯一最优 → `unique`：方向、唯一移位、各段边界、编辑数、CIGAR 及可回放对齐双行串。
- 多个最优 → `ambiguous`：返回稳定排序后的前两份见证及 `more_witnesses`。
- 无可行解释 → HTTP 422 `constraint_failed`，区分
  `segment_length`（结构性长度不可能）与 `per_segment_edit_budget`（预算超限），
  后者给出无预算下最近分段及每个超限分段的位置、所需编辑数与超出量。

## 链方向 `strand_mode`

双链文库的串联条码可能来自参考链或其反向互补链。请求可选 `strand_mode`：

| 取值 | 含义 |
| --- | --- |
| 省略 | 等价于 `forward`，并**完全保持**原有请求、响应与失败体不变（不出现 `strand` 等新字段） |
| `forward` | 仅以提交的参考序列参与解码；显式给出时见证带 `strand: "forward"` |
| `reverse` | 仅以参考序列的**反向互补序列**参与原共同切点解码；`shift` 相对于该反向互补参考 |
| `auto` | 两个方向的**全部可行解释合并**后沿用同一两级目标裁决，不优先保留先计算的方向 |

- `auto` 下仅一个方向可行时，直接判定该方向（`unique`，见证的 `strand` 标明方向），调用方据此
  区分**方向已确定**。
- 两个方向同优时返回 `ambiguous`，见证按 `(strand, shift, boundaries, CIGAR)` 稳定排列
  （`forward` 在前），调用方可识别**方向不确定**。
- 两个方向在预算内均不可行时，仍返回可定位的 `constraint_failed`，表示**读段本身无法解释**；
  `feasible_strands` 给出两个方向各自的可行移位（可能均为空）。
- 非法取值按字段校验错误（422）拒绝。

每份见证的 `shift` 相对于**该方向**的参考序列；`boundaries` 始终使用所提交读段的坐标；
`segments[].reference` 是该方向的定向参考，配合原始读段片段与 CIGAR 可直接回放对齐双行。

## 目录

```
app/solver.py     核心算法：带限 NW + 全最优 CIGAR 枚举、前后缀分段 DP、见证重建、失败诊断（支持正/反链）
app/main.py       FastAPI 服务：POST /api/concatemers/decode 与 /health
tests/            56 个测试，含与穷举参考实现（双向）的一致性校验
scripts/verify.py 一次性校验：等健康 → pytest → 构建自检 → 插/缺/替冒烟 → 正/反链与跨方向歧义冒烟
Dockerfile        python:3.11-slim，内置容器健康检查
docker-compose.yml 可配置宿主机端口；verify 一次性服务（依赖 api 健康后启动）
```

## 用 Docker Compose 启动

```bash
# 默认宿主机端口 8000
docker compose up -d --build

# 自定义宿主机端口
HOST_PORT=9090 docker compose up -d --build
```

API 自带 `/health` 健康检查；compose 也配置了 healthcheck。

## 一次性校验（退出码报告结果）

```bash
docker compose run --build --rm verify
echo $?      # 0 通过，非 0 失败
```

`verify` 服务通过 `depends_on: condition: service_healthy` 等待 API 健康，然后依次：

1. 轮询 `/health`；
2. 在容器内执行全部代码测试（pytest）；
3. 构建自检（依赖版本、应用可导入、路由数）；
4. 对**运行中的 API** 发起解码冒烟，覆盖同一请求中的替换、缺失与插入，
   并额外核对 unique / ambiguous / infeasible 三种响应；
5. 链方向冒烟：反向互补读段在 `forward`/`reverse` 下的不同结论、`auto` 唯一判定方向、
   跨方向同优返回歧义，以及非法 `strand_mode` 被字段拒绝。

## 请求示例

```bash
curl -s -X POST http://localhost:8000/api/concatemers/decode \
  -H 'Content-Type: application/json' \
 -d '{
       "reference": "ACGTACGATC",
       "read": "ATGTACGATCACGTACGATACGTACGATCA",
       "copies": 3,
       "max_edits": 1
     }'
```

上面的读段三段分别含 1 个替换、1 个缺失、1 个插入，响应（节选）：

```json
{
  "status": "unique",
  "objective": {"total_edits": 3, "max_segment_edits": 1},
  "witness": {
    "shift": 0,
    "boundaries": [[0, 10], [10, 19], [19, 30]],
    "segments": [
      {"start": 0, "end": 10, "reference": "ACGTACGATC", "read": "ATGTACGATC",
       "edits": 1, "cigar": "10M",
       "aligned_reference": "ACGTACGATC", "marker": " ^        ", "aligned_read": "ATGTACGATC"},
      {"start": 10, "end": 19, "reference": "ACGTACGATC", "read": "ACGTACGAT",
       "edits": 1, "cigar": "9M1D",
       "aligned_reference": "ACGTACGATC", "marker": "         ^", "aligned_read": "ACGTACGAT-"},
      {"start": 19, "end": 30, "reference": "ACGTACGATC", "read": "ACGTACGATCA",
       "edits": 1, "cigar": "10M1I",
       "aligned_reference": "ACGTACGATC-", "marker": "          ^", "aligned_read": "ACGTACGATCA"}
    ]
  }
}
```

歧义（如同聚物导致移位不可区分）返回：

```json
{
  "status": "ambiguous",
  "objective": {"total_edits": 0, "max_segment_edits": 0},
  "witnesses": [ {…第一份…}, {…第二份…} ],
  "more_witnesses": true
}
```

### 方向判定示例

```bash
curl -s -X POST http://localhost:8000/api/concatemers/decode \
  -H 'Content-Type: application/json' \
 -d '{
       "reference": "ACGTACGATC",
       "read": "GATCGTACGTGATCGTACGTGATCGTACGT",
       "copies": 3,
       "max_edits": 0,
       "strand_mode": "auto"
     }'
```

读段是反向互补参考 `GATCGTACGT` 的三次串联；`auto` 唯一判定方向（节选）：

```json
{
  "status": "unique",
  "objective": {"total_edits": 0, "max_segment_edits": 0},
  "witness": {"strand": "reverse", "shift": 0, "segments": [ {"reference": "GATCGTACGT", …} ]},
  "rotated_reference": "GATCGTACGT"
}
```

若读段同时等距匹配两个方向（例如非周期性的反向互补回文参考），则返回
`"status": "ambiguous"`，两份见证分别带 `"strand": "forward"` 与 `"strand": "reverse"`，
按 `(strand, shift, boundaries, CIGAR)` 稳定排列。

## 本地开发（无 Docker）

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python -m pytest
.venv/bin/uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## 算法说明

- **全局比对**：参考长 ≤20、单段预算 ≤3，使用带宽为 `max_edits` 的带限
  Needleman–Wunsch；随后回溯枚举该距离下的**全部**最优 CIGAR（M/D/I），
  按“操作数优先、再按操作与计数”的确定性顺序排列。结果经 `lru_cache` 复用。
- **分段**：对每个（方向、循环移位）做前缀 DP `(段数, 读段偏移) → (总编辑, 单段最大)`，
  并维护镜像后缀表；只沿“前缀 + 后缀 == 全局最优”的边做 DFS 重建，
  按段长升序、CIGAR 已排序的顺序产出，天然得到稳定排序的见证。
- **方向裁决**：`reverse` 用反向互补参考；`auto` 将两个方向所有 `(方向, 移位)` 的最优值
  汇入同一池后取全局最优，先计算的方向不享有任何优先级。跨方向同优即 `ambiguous`；
  见证排序以方向（`forward` 在前）为最前键。
- **失败诊断**：主流程在预算内无解时，对每个方向与移位计算参考到任意相关子串的
  编辑距离矩阵，再做一次**无单段预算**的分段 DP，找到全局最近分段，
  报告每个超限分段；若连无预算下都无法覆盖（结构性长度约束），直接定位为
  `segment_length`。
