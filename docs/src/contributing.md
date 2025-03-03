# Contributing

GeochemistryTools uses the default [style guide](https://docs.julialang.org/en/v1/manual/style-guide/)
but makes some deviations.

## Changes from default Julia style guide

- `always_for_in = true`
  - `for_in_replacement = ∈`
- `import_to_using = true`
- `whitespace_typedefs = true`
- `whitespace_ops_in_indices = true`
- `remove_extra_newlines = true`
- `short_to_long_function_def = true`
- `always_use_return = true`
- `format_docstrings = true`
- `conditional_to_if = true`
- `normalize_line_endings = "unix"`
- `conditional_to_if = true`

# Naming

Module and type names use `PascalCase`: `module GeometricStatistics`, `struct OrthogonalPolynomial`

Function names should be `snake_case`.
Modifiers should come after the main concept of a function.

- `geometric_mean` instead of`geomean` or `gmean`
- `geometric_mean_with_zeros`

Names should be descriptive, and avoid abbreviations or shorthand unless they are widely
understood or would be excessively verbose. For example:

- `fit_eivlr` instead of `fit_errors_in_variables_regression`

Non-exported internal functions should be prefixed with an underscore: `_evilr_mahon`

Variable names should be `snake_case`, and where applicable follow mathematical conventions.
Unicode characters are not to be used for any public API as this is exclusionary to people
who use terminals that cannot use Unicode characters. They are fine as internal variables.

Names of constants should be UPPERCASE SNAKE_CASE, with an underscore prefix if they are
internal only.