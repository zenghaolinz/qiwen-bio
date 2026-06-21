# 启问 Bio

面向“蛋白结构 - 通路 - 细胞表型”多尺度推理平台的第一阶段 MVP。当前版本先跑通可复现的工程闭环：

`蛋白序列 -> 序列校验 -> 理化特征 -> 演示预测器 -> 证据链 -> Markdown 报告`

## 当前能力

- 校验并规范化标准蛋白质序列（20 种标准氨基酸）
- 计算长度、近似分子量、疏水/带电/芳香残基比例和净电荷代理值
- 通过可替换的 `ProteinPredictor` 接口运行 AMP 演示预测
- 输出带来源、证据类型和置信度的结构化证据链
- 提供 FastAPI、OpenAPI 文档和响应式 Web Demo
- 按 UniProt accession 或人类基因名查询 reviewed UniProt 记录
- 展示功能注释、GO 条目数量及可用的 AlphaFold DB 入口
- 动态获取最新 AlphaFold PDB，统计整体及分档 pLDDT
- 对 `R175H` 形式的突变进行残基一致性校验和位点置信度分析
- 计算非局部 CA 接触图和突变位点 8 Å 空间邻域
- 使用原生 Canvas 展示可拖拽旋转的主链和接触图，无需前端 CDN
- 接入 STRING 高置信互作与过程/通路富集结果
- 输出带来源、互作分数和富集 FDR 的 typed evidence graph
- 使用确定性双环 Canvas 展示蛋白、过程、通路和两类证据边
- 通过 NCBI E-utilities 检索与蛋白和图谱 term 相关的 PubMed 记录
- 提取 PMID、题名、作者、期刊、日期和 DOI，并生成可追踪 Markdown 引用
- 一次请求整合 UniProt、AlphaFold、STRING、PubMed 与序列分析结果
- 输出六层证据完整度评分、缺失层警告和服务端统一 Markdown 报告
- 使用固定 revision 的 ESM-2 8M 模型生成 320 维 mean-pooled embedding
- 按模型、revision、pooling 和序列内容寻址缓存 embedding

> 当前 AMP 预测器是未训练的透明启发式基线，只用于验证系统流程，不可用于科研、临床或实验决策。

## 启动

Windows 可直接双击 `start.bat`，或在终端运行：

```powershell
.\start.bat
```

手动启动方式：

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
uvicorn qiwen_bio.api:app --reload
```

启用真实 ESM-2 embedding：

```powershell
python -m pip install -e ".[model]"
```

模型权重首次使用时从 Hugging Face 下载。若本机设置了不可用的镜像，可临时恢复官方端点：

```powershell
$env:HF_ENDPOINT='https://huggingface.co'
```

浏览器打开 `http://127.0.0.1:8000`，API 文档位于 `http://127.0.0.1:8000/docs`。

## 测试

```powershell
pytest
```

## API 示例

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/analyze `
  -ContentType 'application/json' `
  -Body '{"name":"Demo peptide","sequence":"KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK"}'
```

按基因名或 UniProt ID 查询并分析：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/analyze/uniprot `
  -ContentType 'application/json' `
  -Body '{"identifier":"TP53","organism_id":9606,"mutation":"R175H"}'
```

单独请求 AlphaFold 结构置信度分析：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/structure/alphafold `
  -ContentType 'application/json' `
  -Body '{"accession":"P04637","mutation":"R175H"}'
```

请求 STRING 局部证据图：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/graph/string `
  -ContentType 'application/json' `
  -Body '{"identifier":"TP53","organism_id":9606,"limit":8,"required_score":700}'
```

请求 PubMed 文献记录与 Markdown 引用段落：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/literature/pubmed `
  -ContentType 'application/json' `
  -Body '{"protein":"TP53","context_terms":["Cell cycle","p53 signaling pathway"],"limit":5}'
```

生成统一综合报告：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/report/comprehensive `
  -ContentType 'application/json' `
  -Body '{"identifier":"TP53","organism_id":9606,"mutation":"R175H"}'
```

提取 ESM-2 embedding：

```powershell
Invoke-RestMethod -Method Post -Uri http://127.0.0.1:8000/api/v1/embedding/esm2 `
  -ContentType 'application/json' `
  -Body '{"sequence":"KWKLFKKIEKVGQNIRDGIIKAGPAVAVVGQATQIAK"}'
```

## 后续路线

当前主线是阶段 1C2：整理可再分发的 AMP 数据集、按序列相似性划分数据，并训练经过校准的真实分类器。

模型路线仍需完成：

`qiwen_bio.predictors.ProteinPredictor` 是模型替换边界。真实版本应实现：

1. 离线 ESM-2/ESM-C embedding 缓存和批处理。
2. 使用来源清晰的 AMP 数据集训练分类器。
3. 按序列相似度聚类划分训练、验证、测试集，防止同源泄漏。
4. 保存模型卡、数据版本、阈值、校准与外部测试指标。

结构分析中的接触定义为 CA 距离不超过 8 Å，并排除序列间隔不超过 2 的局部主链邻接。该结果只表示几何邻近，不能直接解释为生化互作、致病性或稳定性变化。

当前通路节点来自 STRING enrichment 返回的 KEGG、Reactome 和 WikiPathways 类别，并不是对这些数据库的直接 API 接入。互作分数和富集 FDR 是数据库证据，不代表因果关系或表型结论。

PubMed 模块只提供上下文检索和书目元数据。检索命中不等于文献支持某个生物结论；使用前仍需阅读摘要或全文并评估研究设计与证据质量。

综合报告的 100 分是证据覆盖度：序列 15、UniProt 注释 20、结构 20、互作 15、过程/通路 15、文献 15。它不表示结论正确率、致病概率、模型置信度或实验成功率。可选数据库失败时报告仍会返回，并明确列出缺失层。

当前 embedding 模型为 `facebook/esm2_t6_8M_UR50D`，固定 revision `c731040fcd8d73dceaa04b0a8e6329b345b0f5df`，最多接收 1022 个残基。缓存位于 `data/embeddings/`，不纳入 Git。当前开发环境安装的是 CPU 版 PyTorch；只有安装 CUDA 版 PyTorch 时 provider 才会自动使用 GPU。Embedding 是特征表示，不是功能或 AMP 分类结论。
