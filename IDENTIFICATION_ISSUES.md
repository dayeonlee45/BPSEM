# Identification Issues in Harring Spline SEM Model

## 발견된 문제점

### 🔴 **문제 1: Phantom Variable과 Intercept의 식별 불가능 (CRITICAL)**

**문제:**
```stan
eta1 = alpha_zero[1] + p_vec + error1
```

이 식에서 `alpha_zero[1]`과 `p_vec`의 **합만 식별**되고, 각각은 식별되지 않습니다.

**예시:**
- Case 1: `alpha_zero[1] = 0.5`, `p_vec = [0.1, 0.2, ...]` → `eta1 = [0.6, 0.7, ...]`
- Case 2: `alpha_zero[1] = 0.0`, `p_vec = [0.6, 0.7, ...]` → `eta1 = [0.6, 0.7, ...]`

두 경우가 **동일한 likelihood**를 생성하므로 무한히 많은 해가 존재합니다.

**해결책:**
- `alpha_zero[1] = 0`으로 고정 (p_vec이 이미 N(0,1)이므로)
- 모델을 `eta1 = p_vec + error1`로 변경

---

### 🟡 **문제 2: 잠재변수 평균의 제약 부족 (MODERATE)**

**문제:**
```stan
matrix[N, M] eta;  // N×M = 300개 파라미터
```

`eta`의 각 열(잠재변수)의 평균이 제약되지 않아, `eta`의 평균이 바뀌면 `alpha_zero`도 함께 바뀌어 같은 likelihood를 생성할 수 있습니다.

**해결책:**
- `transformed parameters`에서 각 잠재변수의 평균을 0으로 제약:
  ```stan
  for (m in 1:M) {
    real mean_eta_m = mean(eta_raw[, m]);
    eta[, m] = eta_raw[, m] - mean_eta_m;
  }
  ```

---

### 🟢 **문제 3: p_vec의 평균 제약 (MINOR)**

**문제:**
`p_vec ~ std_normal()` prior가 있지만, 샘플 평균이 정확히 0이 아닐 수 있습니다.

**해결책:**
- Soft constraint 추가:
  ```stan
  sum(p_vec) ~ normal(0, 0.1 * sqrt(N));
  ```

---

### ⚪ **문제 4: delta=0일 때 tau 식별 불가 (EXPECTED)**

**문제:**
`delta_23 = 0`일 때 `tau23`은 식별되지 않습니다. 이는 **정상적인 동작**입니다 (knot이 효과가 없으므로).

**해결책:**
- 특별한 조치 불필요 (delta=0이면 knot 위치는 의미 없음)

---

## 적용된 수정사항

### 1. **alpha_zero[1] 제거 및 고정**
```stan
// 변경 전
vector[M] alpha_zero;
eta1 = alpha_zero[1] + p_vec + error;

// 변경 후
real alpha_zero_2;  // eta2만
real alpha_zero_3;  // eta3만
eta1 = p_vec + error;  // intercept = 0 고정
```

### 2. **eta 평균 제약**
```stan
// transformed parameters에서
for (m in 1:M) {
  real mean_eta_m = mean(eta_raw[, m]);
  eta[, m] = eta_raw[, m] - mean_eta_m;
}
```

### 3. **p_vec 평균 제약**
```stan
p_vec ~ std_normal();
sum(p_vec) ~ normal(0, 0.1 * sqrt(N));
```

---

## 검증 방법

수정 후 다음을 확인하세요:

1. **R-hat < 1.01**: 모든 파라미터의 R-hat이 1.01 미만
2. **Effective Sample Size**: 충분한 ESS (최소 400 이상)
3. **Trace plots**: 수렴 패턴 확인
4. **Posterior 분산**: 과도하게 큰 SE가 감소했는지 확인

---

## 추가 권장사항

만약 여전히 문제가 있다면:

1. **더 강한 제약**: `sum(p_vec) ~ normal(0, 0.01 * sqrt(N))` (더 엄격)
2. **Hard constraint**: `sum(p_vec) = 0` (transformed parameters에서)
3. **더 많은 iterations**: `iter = 10000, warmup = 5000`
4. **더 큰 샘플 크기**: N=100 → N=200 이상
