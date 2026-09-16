# Evaluate a prediction model on training and test data

Fits a logistic regression model on `train`, evaluates discrimination on
both `train` and `test`, and compares the two areas under the curve with
DeLong's test for unpaired samples.

## Usage

``` r
evaluate_model(train, test, formula, positive = NULL, ci_level = 0.95)
```

## Arguments

- train, test:

  Data frames with the model variables.

- formula:

  Model formula, for example `outcome ~ marker1 + marker2`.

- positive:

  Value of the outcome that represents the event.

- ci_level:

  Confidence level of the reported AUC intervals.

## Value

A tibble with one row: sample sizes and event counts of both data sets,
the training and test AUC with confidence intervals, and the DeLong
p-value comparing them.

## Examples

``` r
set.seed(1)
train <- data.frame(y = rbinom(80, 1, 0.5), x = rnorm(80))
test <- data.frame(y = rbinom(40, 1, 0.5), x = rnorm(40))
evaluate_model(train, test, y ~ x)
#> # A tibble: 1 × 11
#>   n_train n_event_train auc_train train_ci_lower train_ci_upper n_test
#>     <int>         <int>     <dbl>          <dbl>          <dbl>  <int>
#> 1      80            39     0.515          0.387          0.644     40
#> # ℹ 5 more variables: n_event_test <int>, auc_test <dbl>, test_ci_lower <dbl>,
#> #   test_ci_upper <dbl>, delong_p <dbl>
```
