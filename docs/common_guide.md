# StatFaRmer Statistical Methods Reference

> **Purpose**: Statistical methods and model selection logic used in StatFaRmer analysis

## Model Selection Logic

### Automatic Model Selection (`chooseModelType()`)

**Factor Cardinality Rules**:
- **≤10 levels**: Classic ANOVA (balanced designs)
- **11-20 levels**: Mixed models (unbalanced designs, blocking factors)
- **>20 levels**: Fallback to ANOVA (computational efficiency)

**Decision Process**:
1. Analyze factor cardinality from data
2. Select appropriate model type based on rules
3. Execute model with error handling
4. Validate results and assumptions
5. Return comprehensive results

### Model Resolution (`resolveModel()`)

**Process Flow**:
1. **Factor Analysis**: Count levels for each factor
2. **Model Selection**: Apply cardinality rules
3. **Execution**: Run selected model with timeout protection
4. **Validation**: Check convergence and assumptions
5. **Fallback**: Automatic fallback for failed models

## Statistical Models

### Classic ANOVA
- **Use Case**: Balanced designs with ≤10 factor levels
- **Assumptions**: Normality, homogeneity of variances, independence
- **Features**: Fixed effects, two-way interactions
- **Implementation**: `aov()` function

### Mixed Effects Models
- **Use Case**: Unbalanced designs, 11-20 factor levels, blocking factors
- **Random Effects**: `(1|dbscan_cluster)`, `(1|timestamp_group)`
- **Fixed Effects**: Treatment, cultivar, interactions
- **Implementation**: `lme4::lmer()` function
- **Validation**: Convergence checking, random effects validation

### Spline Models
- **Use Case**: Time-series analysis, growth curves
- **Temporal Clustering**: DBSCAN clusters as time points
- **Features**: Smooth trends, treatment comparisons, peak analysis
- **Implementation**: `mgcv::gam()` with spline terms

## Random Effects Validation

### Validation Criteria (`validateRandomEffects()`)
- **Sufficient Observations**: Minimum observations per random effect level
- **Factor Structure**: Proper factor encoding and levels
- **Convergence Feasibility**: Check for singular fits
- **Blocking Factor Detection**: Automatic identification of appropriate random effects

### Blocking Factors
- **`dbscan_cluster`**: Temporal clustering from DBSCAN
- **`timestamp_group`**: Time-based grouping
- **`unit`**: Never used as random effect (design constraint)

## Statistical Assumptions

### Normality Testing
- **Shapiro-Wilk Test**: Primary test for small samples
- **Kolmogorov-Smirnov Test**: Alternative for larger samples
- **QQ Plots**: Visual assessment of normality

### Homogeneity of Variances
- **Levene's Test**: Robust to non-normality
- **Bartlett's Test**: Assumes normality
- **Visual Assessment**: Residual plots

### Independence
- **Design-Based**: Ensured by experimental design
- **Temporal Independence**: Addressed by DBSCAN clustering
- **Spatial Independence**: Addressed by proper randomization

## Post-hoc Analysis

### Tukey HSD (`runTukeyHSD()`)
- **Multiple Comparison Correction**: Controls family-wise error rate
- **Significance Letters**: Generate compact letter display
- **Effect Size Calculations**: Standardized differences
- **Implementation**: `TukeyHSD()` function

### Significance Letters (`generateLetters()`)
- **Compact Letter Display**: Standard statistical notation
- **Group Comparisons**: Clear treatment comparisons
- **Plot Integration**: Compatible with visualization functions

## Model Diagnostics

### Residual Analysis
- **Normality**: QQ plots, normality tests
- **Homoscedasticity**: Residual vs fitted plots
- **Independence**: Temporal and spatial patterns
- **Outliers**: Leverage and influence measures

### Convergence Checking
- **Mixed Models**: Check convergence warnings
- **Singular Fits**: Detect overparameterization
- **Boundary Fits**: Random effect variance near zero
- **Fallback Strategy**: Automatic model simplification

## Performance Considerations

### Computational Efficiency
- **Model Selection**: Automatic selection based on data size
- **Timeout Protection**: Prevent infinite computation
- **Memory Management**: Efficient data handling
- **Caching**: Store expensive model fits

### Large Dataset Handling
- **Median Aggregation**: Robust to outliers, faster computation
- **Factor Reduction**: Combine similar categories
- **Data Subsetting**: Filter for specific analyses
- **Progressive Complexity**: Start simple, add complexity

---

**Implementation**: These methods are used in the [Main Application](app_guide.md) for statistical analysis.

**Data Preparation**: Data is prepared using the [Master Wizard](wizard_guide.md) with appropriate preprocessing for these statistical methods.