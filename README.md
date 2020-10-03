# randoness
Scripts to quantify the plasticity of instantaneous single-cell growth rate as a function of metabolic state

1. M1_low.m
  - quantify pre-shift (low) and post-shift (high) growth rate for M1 (low-adapted cells)
  - single upshift data
  - two replicates

2. M2_high.m
  - quantify pre-shift (high) and post-shift (low) growth rate for M2
  - single downshift data
  - two replicates

3. Mi_intermediate.m
  - quantify post-shift (high) and post-shift (low) growth rate for Mi
  - steady low nutrient shifted to T=60 min fluctuations between high (30 min) and low (30 min)
  - high and low growth rate pairs quantified for 6 successive nutrient periods
  - three replicates

4. plasticity_across_M.m
  - plot growth rate in high and low nutrient from M1, M2 and Mi