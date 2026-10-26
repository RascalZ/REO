# 旧阅读副本的可恢复归档

2026-10-26按用户要求清理旧v1/v2版本。原始投稿历史、科研数据、结果、模型、冻结配置及Git历史不在清理范围。

本地Git忽略的`local_archives/author_review_v1_v2_history_2026-10-26.zip`保留三个旧阅读目录的全部57个非QA文件，包括旧稿MD/Word/PDF、Bib、旧包ZIP、附加文件和清单。压缩包内ARCHIVE_MANIFEST.json记录来源路径、原字节SHA及被丢弃的可再生QA项。当前v4所需补充材料已先逐字节迁入v4，再进行清理。

如需复查历史构建或旧版本，用ZIP恢复到单独工作目录；恢复时保持`manuscript_v1/author_review_v1`、`author_review_v2`和`author_review_v2_zh`相对路径，并核验manifest。不要覆盖当前文件。原历史审计里的路径和SHA描述的是当时状态，不应因为归档而批量改写。QA渲染图可从归档Word/PDF重新生成，不保留为永久重复文件。

v3中英文及当前v4保留。旧v2目录目前仅剩一个被占用而不能删除的旧ZIP，详见cleanup_receipt.json。它已有归档副本；本轮未强行关闭占用进程。

此归档含作者信息，仅保存在本地，不随GitHub或v4投稿包分发。
