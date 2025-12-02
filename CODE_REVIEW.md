# Monte Carlo Simulation 코드 리뷰

## 📋 개요
세 개의 잠재변수(latent variables)로 구성된 모델에서 piecewise structural relation을 검증하기 위한 Monte Carlo 시뮬레이션 코드를 리뷰했습니다.

## 🔍 발견된 주요 문제점

### 1. **코드 중복 문제**
- **문제**: 설정 부분이 두 번 나타남 (처음에 한 번, 나중에 다시 한 번)
- **위치**: 
  - 첫 번째: 라인 1-60 (초기 설정)
  - 두 번째: 라인 200-230 (중복 설정)
- **영향**: 혼란을 야기하고, 두 번째 설정이 첫 번째를 덮어쓸 수 있음
- **해결**: 중복 제거, 단일 설정 블록으로 통합

### 2. **Psi 계산 불일치**
- **문제**: Psi 계산이 두 번 나타나며 서로 다른 시드 사용
  - 첫 번째: `set.seed(12345)`, `n_large = 500000`
  - 두 번째: `set.seed(123)`, `n_large = 100000`
- **영향**: 실제 시뮬레이션에서 사용되는 psi 값이 불명확
- **해결**: 단일 계산 블록으로 통합, 시드 일관성 유지

### 3. **각 Replication 결과 저장 누락**
- **문제**: 사용자가 요청한 "각 replication마다 결과를 저장" 기능이 없음
- **현재**: 마지막에 한 번만 `mc_simulation_results.csv`로 저장
- **해결**: 각 replication마다 개별 CSV 파일 저장 추가

### 4. **결과 저장 정보 부족**
- **문제**: Posterior mean만 저장하고 있음
- **누락된 정보**:
  - Posterior SD
  - 95% 신뢰구간 (2.5%, 97.5%)
  - Rhat 값 (수렴 진단)
  - True parameter values (비교용)
- **해결**: 모든 통계량 포함하도록 확장

### 5. **에러 처리 개선 필요**
- **문제**: 
  - Stan 에러만 체크하고, NULL 체크 없음
  - Summary 추출 실패 시 처리 없음
- **해결**: 더 포괄적인 에러 처리 및 실패 추적 추가

### 6. **수렴 진단 개선**
- **문제**: Rhat < 1.1 체크만 있음
- **개선점**: 
  - 개별 파라미터별 Rhat 저장
  - 수렴 실패 원인 분석 가능하도록
- **해결**: 각 파라미터별 Rhat 저장

### 7. **진행 상황 추적**
- **문제**: 실패한 replication 수 추적 없음
- **해결**: 성공/실패 카운터 추가

## ✅ 개선된 버전의 주요 변경사항

### 1. **코드 구조 정리**
- 중복 제거
- 논리적 블록으로 재구성:
  - User Input
  - Helper Functions & Setup
  - Measurement Model Setup
  - Stan Model Definition
  - Main Simulation Loop
  - Results Combination

### 2. **각 Replication 결과 저장**
```r
# 각 replication마다 개별 파일로 저장
rep_file <- file.path(output_dir, paste0("rep_", sprintf("%04d", r), ".csv"))
write.csv(res_row, rep_file, row.names = FALSE)
```

### 3. **포괄적인 결과 저장**
- Posterior mean, SD
- 95% CI (lower, upper)
- Rhat 값 (각 파라미터별)
- True parameter values (비교용)
- Convergence status

### 4. **향상된 에러 처리**
- Stan fitting 에러 체크
- NULL 체크
- Summary 추출 실패 체크
- 실패한 replication 추적

### 5. **최종 결과 통합**
- 모든 개별 파일을 읽어서 통합
- Summary statistics 출력
- Convergence rate 계산

## 📊 Stan 모델 검토

### ✅ 잘 구현된 부분
1. **Smooth hinge function**: `log1p_exp`를 사용한 부드러운 hinge 구현
2. **Horseshoe prior**: Regularized horseshoe prior로 delta 파라미터에 대한 sparse prior 구현
3. **Lambda 구조**: Marker variable (첫 번째 indicator = 1.0) 올바르게 구현
4. **Structural model**: 
   - Phantom variable (P) → eta1, eta2
   - eta1 → eta3 (piecewise)
   - eta2 → eta3 (linear)

### ⚠️ 주의사항
1. **Sharpness parameter (k=20)**: 
   - 값이 클수록 sharp hinge에 가까워짐
   - 수렴 문제가 발생할 수 있으니 모니터링 필요
   
2. **Prior 설정**:
   - `beta_p ~ normal(0, 0.5)`: 실제 값이 0.3이므로 적절
   - `tau_p ~ normal(0, 1)`: 실제 값이 0.0이므로 적절
   - `delta_raw_p`: Horseshoe prior로 sparse하게 처리

3. **Adaptation**:
   - `adapt_delta = 0.90`: 기본값보다 낮음. 수렴 문제 시 0.95로 증가 고려
   - `max_treedepth = 12`: 추가 권장

## 🎯 추가 권장사항

### 1. **시뮬레이션 설정 저장**
```r
# 시뮬레이션 설정을 별도 파일로 저장
sim_settings <- list(
  N = N, n_reps = n_reps,
  lambda_val = lambda_val,
  beta_P = beta_P, tau_P = tau_P, delta_P = delta_P,
  beta_13 = beta_13, tau_13 = tau_13, delta_13 = delta_13,
  beta_23 = beta_23,
  psi11 = psi11, psi22 = psi22, psi33 = psi33
)
saveRDS(sim_settings, file.path(output_dir, "simulation_settings.rds"))
```

### 2. **진행 상황 저장**
- 중간 결과를 주기적으로 저장 (예: 50 replication마다)
- 시뮬레이션 중단 시 재개 가능하도록

### 3. **수렴 진단 강화**
- ESS (Effective Sample Size) 체크
- Divergence 체크
- Tree depth 경고 체크

### 4. **병렬 처리 고려**
- `foreach` 패키지를 사용한 replication 병렬 처리
- 단, Stan 모델 컴파일은 한 번만 수행

### 5. **결과 시각화 함수**
- Bias 계산 및 시각화
- Coverage rate 계산
- RMSE 계산

## 📝 사용 방법

### 개선된 버전 사용:
```r
source("monte_carlo_simulation_improved.R")
```

### 결과 확인:
```r
# 개별 replication 결과
read.csv("sim_data_1/rep_0001.csv")

# 통합 결과
results <- read.csv("sim_data_1/mc_simulation_results_combined.csv")

# Bias 계산 예시
bias_beta13 <- mean(results$beta13_mean) - results$beta13_true[1]
```

## 🔧 수정이 필요한 부분 (사용자 확인 필요)

1. **Theta_diag 계산**: 
   - 원본: `rep(c(theta_vec, theta_vec, theta_vec))`
   - 개선: `rep(theta_vec, 3)` (더 명확)

2. **Lambda_mat 생성**: 
   - 현재 구조는 올바르지만, 주석으로 명확히 설명 추가

3. **Stan 모델의 beta30**: 
   - Intercept로 보이는데, 데이터 생성 시 intercept가 0인지 확인 필요

## ✨ 결론

전반적으로 코드 구조는 좋지만, 몇 가지 개선이 필요했습니다. 개선된 버전은:
- ✅ 각 replication 결과를 개별 저장
- ✅ 더 포괄적인 통계량 저장
- ✅ 향상된 에러 처리
- ✅ 코드 중복 제거
- ✅ 결과 통합 및 요약

이제 시뮬레이션을 실행하면 각 replication의 상세한 결과를 확인할 수 있고, 최종적으로 통합된 결과도 얻을 수 있습니다.
