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

> 当前 AMP 预测器是未训练的透明启发式基线，只用于验证系统流程，不可用于科研、临床或实验决策。

## 启动

```powershell
python -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
uvicorn qiwen_bio.api:app --reload
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

## 下一阶段接口

`qiwen_bio.predictors.ProteinPredictor` 是模型替换边界。真实版本应实现：

1. 离线 ESM-2/ESM-C embedding 缓存和批处理。
2. 使用来源清晰的 AMP 数据集训练分类器。
3. 按序列相似度聚类划分训练、验证、测试集，防止同源泄漏。
4. 保存模型卡、数据版本、阈值、校准与外部测试指标。
5. 下载 AlphaFold 结构文件并计算 pLDDT、突变位点与功能区距离。
