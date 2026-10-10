# CCLE表达数据审计

本地完整历史RSEM gene TPM可用，无需下载新表达矩阵；原始文件只读。

- TPM：1019个样本 × 57820个基因；样本列与counts完全一致。
- 全基因TPM列和：最小999994.630，中位1000000.110，最大1000006.030；未发现负值或非有限值。
- 所有57820个完整ENSG.version唯一；去版本重复数0。未合并基因、未筛选变异基因。
- matched GENCODE19 GTF精确版本匹配56238；protein_coding 19836。所有基因保留。
- ACH映射1019，有Broad PatientID 1017，ACH重复行0，新旧映射冲突0。全部未映射样本也保留。
- 本地counts已去除标准GCT的版本和维度前导行，是加工过的tab表，不能仅按扩展名当作完整原始GCT。它只用于header审计，不作为本次REO输入。
- 原TPM和log2(TPM+1)均以float64保存，NPY行为样本、列为基因，对应TSV中0基索引；REO可直接读取原TPM以免依赖转换。输入用round-trip浮点解析，未先转成float32。
- 对数变换秩序检查覆盖全部1019个样本、58917561个排序后相邻比较；原相等52608124对、严格不等6309437对，ties改变0、严格次序改变0。排序后相邻检查足以验证所有基因对的相等和大小关系，无需二次方枚举。此项仅保证本次双精度数据及转换，后续不要再转float32构造REO。
- 仅核验本地TPM的结构、数值完整性、版本和官方来源路径，未重新下载官方TPM来逐字节比对，因此不声称本地与远端字节相同。
- Model.csv固定DepMap24Q4 v1，来自Broad Figshare并验证官方MD5及大小；含PatientID但无ParentModelID。PatientID在官方README定义为源于相同组织样本的模型标识，并非完整患者谱系；需联合Sanger parent/patient信息后冻结分组。

## 输出

- `rebuild_v2/data_processed/expression_tpm.npy`
- `rebuild_v2/data_processed/expression_log2tpm.npy`
- `rebuild_v2/data_processed/expression_samples.tsv`
- `rebuild_v2/data_processed/expression_genes.tsv`
- `rebuild_v2/data_processed/expression_audit.json`

下载来源、MD5/SHA256和输入哈希：`rebuild_v2/data_raw/expression/expression_provenance.json`。

## 官方来源

- [Broad CCLE历史文件目录](https://data.broadinstitute.org/ccle/)：同一20180929 gene TPM、counts及配套GENCODE19 GTF。
- [Ghandi et al., Nature 2019](https://doi.org/10.1038/s41586-019-1186-3)：CCLE资源研究。
- [DepMap24Q4 v1固定记录](https://doi.org/10.25452/figshare.plus.27993248.v1)：Model.csv与README，CC BY4.0。
- 2018表达与2024身份注释属于不同用途的固定版本；未混入2024表达重新量化结果。
