flowchart TB
  %% ===== 预处理 =====
  A([开始]) --> B[输入 Sx,Sy,Sz; Ny×Nx]
  B --> C[FFT: Fx,Fy,Fz; pow_map=|Fx|^2+|Fy|^2+|Fz|^2;<br/>P_tot_RMS = sqrt(sum(pow_map)/Ntot)]
  C -->|P_tot_RMS < thr_para=2.5e-7| Z1[[返回 "Para"]]
  C -->|否则| D[构造 kx,ky（fftshift 顺序）；<br/>定位 π 线索引 ix_pi, iy_pi]

  %% ===== π 线占比与直接 Néel =====
  D --> E[计算：<br/>P_lines = π 行(iy_pi,:) + π 列(:,ix_pi) - 交点<br/>P_tot = sum(pow_map)<br/>frac_lines = sqrt(P_lines/(P_tot+ε))]
  E -->|frac_lines ≥ lines_frac_thresh=0.70 且<br/> (P_inter/P_lines) ≥ 0.90| Z2[[返回 "Néel"]]
  E --> F{frac_lines ≥ 0.70 ?}
  F -->|是| F1[候选 cand = π 行 ∪ π 列；<br/>denom = P_lines]
  F -->|否| F2[候选 cand = 全部网格点；<br/>denom = P_tot]

  %% ===== B3：主模选取（±q 成对直到95%）=====
  F1 --> G
  F2 --> G
  G[按 pow_map 从大到小排序 cand；<br/>若 weak: P_tot_RMS < 5e-4，<br/>则仅保留 pow ≥ 0.15·P_max] --> H[自高到低遍历：<br/>对每个 q 去重（含 ±q），<br/>累积配对功率 P_acc，直到 P_acc/denom ≥ 0.95]
  H --> I{收集的主模数 l}

  %% ===== B4：三模及以上 → 共线性/Strange =====
  I -->|l ≥ 3| J[计算 collinearity_rms(vecs[0])]
  J -->|rms < thr_col = 2.5e-7| Z3[[返回 "Other-collinear"]]
  J -->|否则| Z4[[返回 "Strange"]]

  %% ===== 单模分支 =====
  I -->|l == 1| K[第一主模 (iy1,ix1), v1, q1=(qx1,qy1)]
  K --> L{η = mod(π - (iy1==iy_pi ? qx1 : qy1), 2π)/(2π) < eta_tol = 5e-3 ?}
  L -->|是| Z5[[返回 "Néel"]]
  L -->|否| M[计算：<br/>rms = collinearity_rms(v1)<br/>r_perp = transverse_ratio_after_align(v1)]
  M -->|r_perp < 0.08 或 rms < thr_col| Z6[[返回 "Stripe"]]
  M -->|否则| Z7[[返回 "Spiral"]]

  %% ===== 双模及以上分支：Beat / IBS（CpBS/ClBS）=====
  I -->|l ≥ 2| N[第二主模 (iy2,ix2), q2=(qx2,qy2)]
  N --> O[判断是否在同一轴：<br/>same_x: on_qx_axis(q1)&on_qx_axis(q2)<br/>same_y: on_qy_axis(q1)&on_qy_axis(q2)<br/>distinct_along_axis(q1,q2)]
  O --> P{(same_x 或 same_y) 且 沿轴可区分?}
  P -->|是| Q[尝试“救援”IBS：<br/>若 same_x，则看 qy 轴是否存在<br/>axis_presence_frac ≥ 0.10 的峰并替换第二主模；<br/>若 same_y，则看 qx 轴是否存在相应峰]
  Q --> R{救援成功?}
  R -->|否| Z8[[返回 "Beat"]]
  R -->|是| S[继续]
  P -->|否| S

  %% ===== 正交轴与功率均衡检查 =====
  S --> T{是否正交轴: <br/>(q1 在 qx 轴且 q2 在 qy 轴) 或反之?}
  T -->|否| Z9[[返回 "Strange"]]
  T -->|是| U[计算配对功率：a1=pair_power(q1), a2=pair_power(q2)<br/>balance = |a1-a2|/(a1+a2)<br/>ibs_power_frac=(a1+a2)/denom]
  U --> V{阈值通过？<br/>balance ≤ balance_max，ibs_power_frac ≥ ibs_power_floor}
  V -->|否| Z10[[返回 "Strange"]]
  V -->|是| W[将 v1 对齐 x 轴：得 v1r,v2r；<br/>dot=|<v1r,v2r>|，amp_equal=|‖v2r‖-1|≤amp_match_tol=5e-3；<br/>且 |v1r_y|, |v1r_z| ≤ coplanar_tol=5e-3]
  W --> X{dot ≤ coplanar_tol 且 amp_equal 且<br/>|v1r_y|≤tol 且 |v1r_z|≤tol ?}
  X -->|是| Z11[[返回 "CpBS"]]
  X -->|否| Z12[[返回 "ClBS"]]
