# Repeated-Measures Agreement Analysis
This project implements a mixed-effects approach to agreement analysis using repeated PEFR measurements from the classic Bland & Altman dataset:
🔗 https://www-users.york.ac.uk/~mb55/meas/ba.htm
It evaluates:
* Device agreement (between methods)
* Repeatability (within device across occasions)

The approach extends classical Bland–Altman analysis to handle repeated measurements per subject in a single hierarchical mixed model:
value = μ + device + occasion 
        + subject + (subject × device) + (subject × occasion) 
        + residual

This separates:

* Between-subject variation
* Within-subject (measurement) variation
* Device and occasion effects

Two Modelling Approaches are employed:

1️⃣ Homogeneous Residual Variance
e ~ N(0, σ²)

This assumes that all measurements share the same residual variance.
LoA = bias ± 1.96 × √(2σ²)

### Pros
* Simple and stable
* Matches classical Bland–Altman
* Easy to interpret

### Cons
* Assumes equal precision across devices
* Cannot detect heteroscedasticity

2️⃣ Heterogeneous Residual Variance (by Device)
e ~ N(0, σ²_device)

In this model each device has its own variance.
LoA = bias ± 1.96 × √(σ²₁ + σ²₂)

### Pros
* More realistic (matches measurement theory)
* Identifies differences in device precision
* Better representation of agreement

### Cons
* More complex
* Requires more data
* Less stable in bootstrap

Whichever model is used, the code provides a range of agreement and repeatability statistics derived from model variance components:

**Device Agreement**
* Bias
* Limits of Agreement (LoA)
* Repeatability Coefficient (RC)
* Coefficient of Variation (CV)
* Lin’s Concordance Correlation Coefficient (CCC)

**Occasion Repeatability**
* Intraclass Correlation (ICC)
* Bias
* LoA
* RC
* CV

## Comparison of approaches
| Feature                  | Homogeneous Variance | Heterogeneous Variance |
|-------------------------|--------------------|-----------------------|
| Residual variance       | Single σ²          | Device-specific σ²    |
| Simplicity              | ✅ High            | ❌ Lower              |
| Stability               | ✅ High            | ⚠ Moderate            |
| Realism                 | ⚠ Limited          | ✅ Improved           |
| PEFR suitability        | ⚠ Approximate      | ✅ Better             |
| LoA calculation         | Simple             | More complex          |
| Interpretation          | Straightforward    | More nuanced          |
👉 Therefore:

When to use each approach:
| Scenario                          | Recommended Model |
|----------------------------------|------------------|
| Quick exploratory analysis       | Homogeneous      |
| Small sample size                | Homogeneous      |
| Need simple Bland–Altman output  | Homogeneous      |
| Devices differ in precision      | Heterogeneous    |
| Evidence of heteroscedasticity   | Heterogeneous    |
| High-quality inference required  | Heterogeneous    |
| Comparing device performance     | Heterogeneous    |

## Key Takeaway
Homogeneous  → simple, robust, approximate  
Heterogeneous → realistic, detailed, data-dependent 
