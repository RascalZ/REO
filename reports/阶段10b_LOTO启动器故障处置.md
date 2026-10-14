# 阶段10b：LOTO 启动器日志模式故障处置

日期：2026-10-14。v1.1 的正式无拟合预检和独立验证已经通过；第一次启动控制器时，在任何模型worker运行前暴露Windows资源采样CSV字段模式缺陷。

## 观察到的故障

Windows平台的`process_tree_memory()`除CSV表头已声明的资源字段外，还返回`system_commit_bytes`和`system_commit_limit_bytes`。v1.1控制器将完整内存字典写入CSV，但字段名列表漏掉这两列，首个资源采样即触发`ValueError: dict contains fields not in fieldnames`。失败点在`run_controller()`的资源采样阶段，子进程尚未启动。

本地保留的v1.1运行记录显示：`status=failed`、0/1,472个任务完成、0个worker秒、0个检查点、无拟合事件文件、无OOF预测文件；没有模型拟合。该运行目录仍保留为审计证据，不继续或覆盖。

## 版本处置

建立新的LOTO执行修订v1.2和独立run目录，只修复控制器资源CSV schema，加入Windows内存API返回的两列。药物、组织、模型、样本、训练/测试拆分、种子、特征选择、调参网格、评价指标、NA/缺失政策、任务规模和资源阈值均沿用v1.1。v1.2须重新生成不可变锁，重新完成无拟合预检和独立拆分/清单验证，之后才开始模型拟合。

故障执行提交：`1f5b271`（计划/预检文档记录提交）；此前已冻结的v1.1代码及锁仍在Git历史中，原锁SHA-256为`689d7fe85610701b2b91e4fc4d089bdd478dfee545ad24eda9fc926c237460e7`。v1.1模型结果为零个；CTRP未访问。此为软件日志缺陷，不是研究假设或模型结果，不得解释为模型失败。

## v1.2重新冻结与正式预检

修正代码在提交`a286df384fa666250963ecf26d5fcfbf4cae5b3f`中，LOTO v1.2执行锁SHA-256为`78e06d75904fee4023090978d476a5970110d69b16a77ae655a0eb191611c697`。新run目录`runs/gdsc2-cancer-v2-fixed8-tissue-loto-v1_2`完成了无拟合预检；独立validator返回`preflight_validated_no_model_fits`。锁定工作量仍为1,472任务、40,112拟合、39,736外层预测；184个组织×药物拆分均通过训练/测试组隔离和inner-fold完整性检查，内层manifest 109,274行。稀有类警告与v1.1完全相同：176个上下文至少一个组织只有1组，8个上下文有2组，均保留原样。

v1.2正式preflight manifest SHA-256：`0031198974ede4abd7f5bc4741c8455f755930741d134146258a786ac8a10222`；task plan为`9862bb25d06e32a1759792ed248a82ce4ef167d23c95e046608561be62b1a484`；inner manifest为`0ef218bfec2cc9d0e3012f1f4bd3aab61867d2b3c10fa93385ad032ad35e96b8`。开发响应输入仍与P5e冻结来源相同，CTRP响应文件未打开。预检通过，v1.2现在具备启动模型拟合的条件。
