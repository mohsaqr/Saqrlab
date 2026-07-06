# Simulate Survival Data (Cox Proportional Hazards)

Generate right-censored survival data from a Cox proportional-hazards
generative model. The hazard for subject \\i\\ is \$\$h(t \mid X_i) =
h_0(t) \exp(X_i \beta),\$\$ where \\h_0(t)\\ is a parametric baseline
hazard (Weibull, exponential, or Gompertz). Event times are drawn by
inverse-CDF sampling given each subject's linear predictor, then
administrative/exponential censoring is applied so that roughly
`censoring_rate` of observations are censored. Designed so that
`coxph(Surv(time, status) ~ ., data)` recovers `betas` at large `n`.

## Usage

``` r
simulate_survival(
  n = 300,
  n_covariates = 2,
  betas = NULL,
  baseline = c("weibull", "exponential", "gompertz"),
  lambda = 0.1,
  shape = 1,
  censoring_rate = 0.3,
  covariate_type = c("normal", "binary"),
  seed = NULL
)
```

## Arguments

- n:

  Integer. Number of subjects. Default: 300.

- n_covariates:

  Integer. Number of covariates \\p\\. Default: 2. Ignored when `betas`
  is supplied (its length wins).

- betas:

  Numeric vector or NULL. True log hazard ratios, one per covariate.
  When NULL, a moderate alternating sequence is auto-generated (e.g.
  `c(0.5, -0.3, 0.5, -0.3, ...)`).

- baseline:

  Character. Baseline hazard family: `"weibull"` (default),
  `"exponential"`, or `"gompertz"`.

- lambda:

  Positive numeric. Baseline scale parameter. Default: 0.1.

- shape:

  Positive numeric. Baseline shape parameter. For Weibull this is the
  Weibull shape; for Gompertz the level parameter; ignored for
  exponential. Default: 1.

- censoring_rate:

  Numeric in \[0, 1). Target proportion of censored observations.
  Default: 0.3.

- covariate_type:

  Character. Covariate distribution: `"normal"` (default, standard
  normal) or `"binary"` (Bernoulli(0.5), coded 0/1).

- seed:

  Integer or NULL. Random seed.

## Value

A [`saqr_sim`](https://pak.dynasite.org/Saqrlab/reference/saqr_sim.md)
object with:

- `$data`:

  data.frame with columns `time` (numeric event or censoring time),
  `status` (1 = event, 0 = censored), and `X1`...`Xp` (covariates).

- `$params`:

  list with `betas` (true log hazard ratios), `baseline`, `lambda`,
  `shape`, `covariate_type`, and `realized_censoring_rate`.

## Details

Inverse-CDF event times use the closed-form survival inverse for each
baseline. With \\U \sim \mathrm{Unif}(0,1)\\ and linear predictor \\\eta
= X\beta\\: Weibull \\t = (-\log U / (\lambda
e^{\eta}))^{1/\mathrm{shape}}\\; exponential \\t = -\log U / (\lambda
e^{\eta})\\; Gompertz \\t = (1/\mathrm{shape}) \log(1 -
\mathrm{shape}\log U / (\lambda e^{\eta}))\\.

## Examples

``` r
r <- simulate_survival(n = 400, n_covariates = 2, seed = 1)
print(r)
#> saqr_sim [survival]  400 x 4  (seed=1)
#>   params: betas, baseline, lambda, shape, covariate_type, realized_censoring_rate 
#>   cols:   time, status, X1, X2 
r$params$betas
#>   X1   X2 
#>  0.5 -0.3 
r$params$realized_censoring_rate
#> [1] 0.2975
head(as.data.frame(r))
#>         time status         X1         X2
#> 1 15.3385827      0 -0.6264538  1.0744410
#> 2 19.1493259      0  0.1836433  1.8956548
#> 3  4.2603715      0 -0.8356286 -0.6029973
#> 4  0.3510546      1  1.5952808 -0.3908678
#> 5  1.3033325      1  0.3295078 -0.4162220
#> 6 18.5459318      0 -0.8204684 -0.3756574

# Binary covariates, exponential baseline, explicit effects
r2 <- simulate_survival(n = 500, betas = c(0.8, -0.5),
                        baseline = "exponential",
                        covariate_type = "binary", seed = 42)
summary(r2)
#> Simulation type: survival
#> Seed: 42
#> 
#> --- Data ---
#>       time               status            X1             X2       
#>  Min.   : 0.002791   Min.   :0.000   Min.   :0.00   Min.   :0.000  
#>  1st Qu.: 1.863159   1st Qu.:0.000   1st Qu.:0.00   1st Qu.:0.000  
#>  Median : 4.211101   Median :1.000   Median :0.00   Median :0.000  
#>  Mean   : 6.916516   Mean   :0.674   Mean   :0.47   Mean   :0.476  
#>  3rd Qu.: 9.141489   3rd Qu.:1.000   3rd Qu.:1.00   3rd Qu.:1.000  
#>  Max.   :51.034301   Max.   :1.000   Max.   :1.00   Max.   :1.000  
#> 
#> --- Parameters ---
#> List of 6
#>  $ betas                  : Named num [1:2] 0.8 -0.5
#>   ..- attr(*, "names")= chr [1:2] "X1" "X2"
#>  $ baseline               : chr "exponential"
#>  $ lambda                 : num 0.1
#>  $ shape                  : num 1
#>  $ covariate_type         : chr "binary"
#>  $ realized_censoring_rate: num 0.326
```
