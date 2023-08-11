# Style Guide

When developing Jems, please take into account the following style pointers.

We use 4 spaces as the indent marker, and use a line length of 120. This is not strictly enforced, but try to keep
overrunning lines to a minimum.

Docstrings are demarked by three double quotes:

```julia
"""
    my_func(a::Number)

This function does something cool with number `a`.
"""
```

Comments are marked with the hashtag and a space:

```julia
# this loop does good stuff!
for i = 1:10
    a += 1
end
```

Inline comments should be separated by at least two spaces:

```julia
c = a + b  # this is high level stuff!
```

## JuliaFormatter

The file `.JuliaFormatter.toml` is be used in conjunction with `JuliaFormatter.jl` to automatically format source files
according to our adopted style:

```julia
JuliaFormatter.format(".")  # formats the whole directory of source files
```

The main function of the formatter is that it will automatically fold long lines, and inserts spaces around operators.
One disadvantance is that it ignores comments and does not yet handle docstrings (even if the `.toml` file explicitly
says to include docstrings, this is a [bug](https://github.com/domluna/JuliaFormatter.jl/issues/649)).
