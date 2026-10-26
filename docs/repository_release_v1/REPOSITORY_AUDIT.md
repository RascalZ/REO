# 科研仓库规范检查

日期：2026-10-26。对象：RascalZ/REO，检查起点595161c。结论：**学术代码仓库的基本结构已具备并完成本轮整理；尚不是完整公开复现发布。**

## 当前访问及权限

GitHub连接器实时返回仓库为private，当前账号RascalZ的权限为admin，且有push/maintain/pull/triage。main的protected字段为false；最近五次已有CI均成功；未发现已创建的GitHub Release。见github_access_audit.json。

账号权限与连接器可执行权限不同：分支保护管理端点返回403 “Resource not accessible by integration”。因此不能凭账号admin就声称当前连接器能更改所有设置。本轮通过文件和工作流完成规范化，未改公开状态、协作者权限或分支保护。工作流的默认token权限已收紧为contents: read。实际公开或指定审稿人访问仍需单独落实。

## 本轮已完成

| 项目 | 结果 |
|---|---|
| 首页 | 改为研究问题、结果范围、两条quickstart、访问表、复现/引用/许可入口；原README逐字节保留 |
| 复现导航 | 新增根目录REPRODUCING.md，区分toy/聚合检查、真实模型推理与完整训练；不把现有环境检查称作新环境安装 |
| CI | 新增推理作业，运行公开资产摘要检查和16项手写fixture合同测试；保留原有科学合同检查 |
| 引用 | 保留CFF及真实开发版本，移除未经正式Release支持的date-released；未虚构论文、DOI、个人作者或归档标识 |
| 许可 | 现有MIT及排除上游数据/原稿/legacy的边界保持；新增贡献说明 |
| 入口验证 | 从Git已跟踪文件导出干净代码副本，用现有固定Python3.12环境执行6条命令全部通过；toy预测与手算期望一致 |
| 敏感模式 | 92个可达提交中633个允许文本blob、约16.3MB，13类常见密钥/私钥模式候选0；排除数据、CSV、归档和其他非允许对象 |

历史扫描范围和排除项见history_secret_scan.json。该有界扫描不是完整无敏感信息保证；未读取工作区原始CTRP响应、索引或个体预测。新检查不构成原始科学结果独立复现。Git科学源、冻结配置和结果不变。

## 仍待真实完成的发布事项

- 仓库仍私有；提供实际编辑/审稿人或公众访问方式。
- 64个固定模型、完整衍生身份/折分及预处理状态仍在ignored本地资产中。需要许可核对和真实托管/受控获取，不能用README或SHA清单代替。
- 确认代码贡献者的正式引用信息；当前集体贡献者元数据不冒充已确认的个人代码作者。
- 在真实发布时确定不可变版本及归档标识；目前没有GitHub Release或归档DOI。未在本轮生成假发布事实。

这些事项不妨碍当前代码审阅和材料整理，但不能将当前仓库称为完整开放模型基准。三处历史文稿构建脚本仍有本机Node fallback，均先查PATH；本轮保留其历史字节，不为风格改写冻结/历史代码。

## v4包与清理

v4投稿材料ZIP含26项，已有主稿字节不变；补充材料全部迁入v4，不再依赖旧v2目录。包不是最终作者声明及真实访问安排的替代品，阅读版标记仍在原稿。见submission_package_manifest.json及ZIP内00_READ_BEFORE_SUBMISSION.txt。

v1、v2及旧中文v2的57个非QA文件已逐项SHA和ZIP CRC核验后归档到本地local_archives；187个可再生成QA文件被清理。共删除243个文件，移除38,032,860字节旧散文件；一个4,765,224字节的旧v2 ZIP被其他进程占用而保留。v3/v4、科研数据、结果、模型、配置、Git历史及本轮前已有其他工作树修改未删。实际清理见cleanup_receipt.json。

## 依据

- [BMC Genomics研究论文指南](https://link.springer.com/journal/12864/submission-guidelines/research-article)：资料获取和作者声明需可核实；本轮不冒充已满足尚未完成项。
- [GitHub引用文件说明](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/about-citation-files)：根目录CITATION.cff提供软件引用入口；无论文/DOI时不填虚构的preferred-citation。
- fresh_checkout_checks.json记录真实命令、环境、退出码及边界；GitHub CI的本轮新结果在实际完成后单独记录，不能用之前成功代替。
