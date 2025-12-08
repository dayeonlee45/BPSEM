# Identification Analysis for Harring Spline SEM Model
# 
# 주요 Identification 문제점 분석:

# ==============================================================================
# 문제 1: Phantom Variable (p_vec) Identification
# ==============================================================================
# 
# 모델에서:
#   eta1 = alpha_zero[1] + p_vec + error1
#   eta2 = alpha_zero[2] + spline(p_vec) + error2
#
# 문제: alpha_zero[1]과 p_vec의 합이 식별되지만, 각각은 식별되지 않음!
# 
# 해결책: 
#   - Option 1: alpha_zero[1] = 0으로 고정 (p_vec의 평균이 0이므로)
#   - Option 2: p_vec의 평균을 0으로 제약 (이미 std_normal() prior이지만 명시적 제약 필요)
#   - Option 3: eta1의 intercept를 0으로 고정하고 p_vec만 추정

# ==============================================================================
# 문제 2: Latent Variables (eta) Identification  
# ==============================================================================
#
# eta는 N×M 행렬 (N=100, M=3 → 300개 파라미터)
# Measurement model: Y = eta * Lambda' + error
#
# 식별 조건:
#   - 각 잠재변수에 marker variable (loading=1)이 있음 ✓
#   - 잠재변수들의 scale이 고정되어야 함
#   - 하지만 eta의 평균이 제약되지 않음!
#
# 문제: eta의 평균이 자유롭게 움직일 수 있음
#   - eta의 평균이 바뀌면 alpha_zero도 함께 바뀌어 같은 likelihood 생성 가능
#
# 해결책:
#   - eta의 평균을 0으로 제약 (transformed parameters에서)
#   - 또는 measurement model에서 관측변수의 평균이 0이므로 자동으로 제약됨 (이미 center=TRUE)

# ==============================================================================
# 문제 3: Spline Parameter Identification
# ==============================================================================
#
# Harring spline: alpha2 (average slope), alpha3 (half-difference), tau (knot)
#
# 문제: 
#   - alpha3가 0에 가까우면 tau가 식별되지 않음 (flat line)
#   - alpha2와 alpha3의 조합이 여러 방식으로 표현 가능할 수 있음
#
# 현재 코드에서:
#   - delta_23 = 0이므로 alpha3_23 = 0
#   - 이 경우 tau23이 식별되지 않음 (정상, delta=0이면 knot 의미 없음)

# ==============================================================================
# 문제 4: Scale Identification
# ==============================================================================
#
# 잠재변수들의 scale:
#   - eta1, eta2, eta3: marker variable로 scale 고정 ✓
#   - p_vec: std_normal() prior로 scale 고정 ✓
#
# 하지만 alpha_zero와 eta의 관계에서:
#   - eta1 = alpha_zero[1] + p_vec
#   - alpha_zero[1]이 자유롭게 움직이면 p_vec의 scale이 간접적으로 영향받을 수 있음

# ==============================================================================
# 가장 심각한 문제: Phantom Variable + Intercept Identification
# ==============================================================================
#
# eta1 = alpha_zero[1] + p_vec + error1
#
# 이 식에서:
#   - alpha_zero[1] + p_vec의 합만 식별됨
#   - 각각은 식별되지 않음 (infinite solutions)
#
# 예: alpha_zero[1] = 0.5, p_vec = [0.1, 0.2, ...]
#     vs alpha_zero[1] = 0.0, p_vec = [0.6, 0.7, ...]
#     → 같은 eta1 예측값 생성!

# ==============================================================================
# 해결책
# ==============================================================================

cat("
=== Identification Issues Found ===

1. CRITICAL: Phantom variable (p_vec) and alpha_zero[1] are not separately identified
   - eta1 = alpha_zero[1] + p_vec + error
   - Only the sum is identified, not each component
   - Fix: Set alpha_zero[1] = 0 (since p_vec ~ N(0,1))

2. MODERATE: Latent variable means (eta) may not be fully constrained
   - eta matrix has N×M parameters
   - Measurement model centers Y, but eta means still free
   - Fix: Constrain mean(eta[,m]) = 0 in transformed parameters

3. MINOR: When delta=0, tau is not identified (expected behavior)
   - tau23 when delta_23=0: not identified (OK, knot has no effect)

=== Recommended Fixes ===

Fix 1: Constrain alpha_zero[1] = 0
Fix 2: Add constraint: mean(eta[,m]) = 0 for each m
Fix 3: Consider fixing p_vec mean to 0 explicitly (already has std_normal prior)
")
